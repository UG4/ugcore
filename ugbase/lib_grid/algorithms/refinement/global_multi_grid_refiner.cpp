/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#include <cassert>
#include "common/profiler/profiler.h"
#include "global_multi_grid_refiner.h"
#include "lib_grid/algorithms/algorithms.h"
#include "lib_grid/file_io/file_io.h"

//define PROFILE_GLOBAL_MULTI_GRID_REFINER if you want to profile
//the refinement code.
#define PROFILE_GLOBAL_MULTI_GRID_REFINER
#ifdef PROFILE_GLOBAL_MULTI_GRID_REFINER
	#define GMGR_PROFILE_FUNC()	PROFILE_FUNC_GROUP("grid")
	#define GMGR_PROFILE(name)	PROFILE_BEGIN_GROUP(name, "grid")
	#define GMGR_PROFILE_END()	PROFILE_END()
#else
	#define GMGR_PROFILE_FUNC()
	#define GMGR_PROFILE(name)
	#define GMGR_PROFILE_END()
#endif

using namespace std;

namespace ug
{

GlobalMultiGridRefiner::
GlobalMultiGridRefiner(IRefinementCallback* refCallback) :
	IRefiner(refCallback),
	m_pMG(NULL)
{
}

GlobalMultiGridRefiner::
GlobalMultiGridRefiner(MultiGrid& mg, IRefinementCallback* refCallback) :
	IRefiner(refCallback)
{
	m_pMG = NULL;
	assign_grid(mg);
}

GlobalMultiGridRefiner::~GlobalMultiGridRefiner()
{
	if(m_pMG)
		m_pMG->unregister_observer(this);
}

void GlobalMultiGridRefiner::grid_to_be_destroyed(Grid* grid)
{
	m_pMG = NULL;
}

void GlobalMultiGridRefiner::assign_grid(MultiGrid& mg)
{
	assign_grid(&mg);
}

void GlobalMultiGridRefiner::assign_grid(MultiGrid* mg)
{
	if(m_pMG){
		m_pMG->unregister_observer(this);
		m_pMG = NULL;
	}
	
	if(mg){
		m_pMG = mg;
		set_message_hub(mg->message_hub());
		m_pMG->register_observer(this, OT_GRID_OBSERVER);
	}
}


void GlobalMultiGridRefiner::
num_marked_edges_local(std::vector<int>& numMarkedEdgesOut)
{
	num_marked_elems<Edge>(numMarkedEdgesOut);
}

void GlobalMultiGridRefiner::
num_marked_faces_local(std::vector<int>& numMarkedFacesOut)
{
	num_marked_elems<Face>(numMarkedFacesOut);
}

void GlobalMultiGridRefiner::
num_marked_volumes_local(std::vector<int>& numMarkedVolsOut)
{
	num_marked_elems<Volume>(numMarkedVolsOut);
}


template <class TElem>
void GlobalMultiGridRefiner::
num_marked_elems(std::vector<int>& numMarkedElemsOut)
{
	numMarkedElemsOut.clear();
	if(!m_pMG)
		return;
	numMarkedElemsOut.resize(m_pMG->num_levels(), 0);
	if(m_pMG->num_levels() > 0)
		numMarkedElemsOut.back() = m_pMG->num<TElem>(m_pMG->top_level());
}

////////////////////////////////////////////////////////////////////////
void GlobalMultiGridRefiner::perform_refinement()
{
	UG_DLOG(LIB_GRID, 1, "GlobalMultiGridRefiner\n");

	GMGR_PROFILE_FUNC();
		
	assert(m_pMG && "refiner has to be assigned to a multi-grid!");
	if(!m_pMG)
		return;

//	the multi-grid
	MultiGrid& mg = *m_pMG;

//	check if a refinement-callback is set.
//	if not, we'll automatically set one, if a position attachment is available
	bool localRefCallbackSet = false;
	if(!m_refCallback){
		if(mg.has_vertex_attachment(aPosition)){
			localRefCallbackSet = true;
			m_refCallback = new RefinementCallbackLinear<APosition>(mg, aPosition);
		}
		else if(mg.has_vertex_attachment(aPosition2)){
			localRefCallbackSet = true;
			m_refCallback = new RefinementCallbackLinear<APosition2>(mg, aPosition2);
		}
		else if(mg.has_vertex_attachment(aPosition1)){
			localRefCallbackSet = true;
			m_refCallback = new RefinementCallbackLinear<APosition1>(mg, aPosition1);
		}		
	}

//	make sure that the required options are enabled.
	if(mg.num_volumes() > 0){
		if(!mg.option_is_enabled(VOLOPT_AUTOGENERATE_FACES))
		{
			LOG("WARNING in GlobalMultiGridRefiner::refine(): auto-enabling VOLOPT_AUTOGENERATE_FACES.\n");
			mg.enable_options(VOLOPT_AUTOGENERATE_FACES);
		}
	}
	
	if(mg.num_faces() > 0){
		if(!mg.option_is_enabled(FACEOPT_AUTOGENERATE_EDGES))
		{
			LOG("WARNING in GlobalMultiGridRefiner::refine(): auto-enabling FACEOPT_AUTOGENERATE_EDGES.\n");
			mg.enable_options(FACEOPT_AUTOGENERATE_EDGES);
		}
	}

	if(mg.num_levels() == 0)
		return;

//	the old top level
	int oldTopLevel = mg.num_levels() - 1;
//todo: Ghosts have to be ignored in the specified grid-object-collection - otherwise problems may occur
//		e.g. in the AdaptionSurfaceGridFunction during prolongation of piecewise const functions... (M. Breit)
	m_messageHub->post_message(GridMessage_Adaption(GMAT_GLOBAL_REFINEMENT_BEGINS,
													mg.get_grid_objects(oldTopLevel)));

	UG_DLOG(LIB_GRID, 1, "REFINER: reserving memory...");

//	reserve enough memory to speed up the algo
	GMGR_PROFILE(GMGR_Reserve);
	{
		int l = oldTopLevel;

		GMGR_PROFILE(GMGR_ReserveVrtData);
		mg.reserve<Vertex>(mg.num<Vertex>() +
					+ mg.num<Vertex>(l) + mg.num<Edge>(l)
					+ mg.num<Quadrilateral>(l) + mg.num<Hexahedron>(l));
		GMGR_PROFILE_END();

		GMGR_PROFILE(GMGR_ReserveEdgeData);
		mg.reserve<Edge>(mg.num<Edge>()
					+ 2 * mg.num<Edge>(l) + 3 * mg.num<Triangle>(l)
					+ 4 * mg.num<Quadrilateral>(l) + 3 * mg.num<Prism>(l)
					+ mg.num<Tetrahedron>(l)
					+ 4 * mg.num<Pyramid>(l) + 6 * mg.num<Hexahedron>(l));
		GMGR_PROFILE_END();

		GMGR_PROFILE(GMGR_ReserveFaceData);
		mg.reserve<Face>(mg.num<Face>()
					+ 4 * mg.num<Face>(l) + 10 * mg.num<Prism>(l)
					+ 8 * mg.num<Tetrahedron>(l)
					+ 9 * mg.num<Pyramid>(l) + 12 * mg.num<Hexahedron>(l));
		GMGR_PROFILE_END();

		GMGR_PROFILE(GMGR_ReserveVolData);
		mg.reserve<Volume>(mg.num<Volume>()
					+ 8 * mg.num<Tetrahedron>(l) + 8 * mg.num<Prism>(l)
					+ 6 * mg.num<Pyramid>(l) + 8 * mg.num<Hexahedron>(l));
		GMGR_PROFILE_END();
	}
	GMGR_PROFILE_END();
	UG_DLOG(LIB_GRID, 1, " done.\n");

	UG_DLOG(LIB_GRID, 1, " refinement begins.\n");
//	notify derivates that refinement begins
	refinement_step_begins();
	
//	cout << "num marked edges: " << m_selMarks.num<Edge>() << endl;
//	cout << "num marked faces: " << m_selMarks.num<Face>() << endl;

//	we want to add new elements in a new layer.
	bool bHierarchicalInsertionWasEnabled = mg.hierarchical_insertion_enabled();
	if(!bHierarchicalInsertionWasEnabled)
		mg.enable_hierarchical_insertion(true);


//	some buffers
	vector<Vertex*> vVrts;
	vector<Vertex*> vEdgeVrts;
	vector<Vertex*> vFaceVrts;
	vector<Edge*>	vEdges;
	vector<Face*>		vFaces;
	vector<Volume*>		vVols;
	
//	some repeatedly used objects
	EdgeDescriptor ed;
	FaceDescriptor fd;
	VolumeDescriptor vd;

	UG_DLOG(LIB_GRID, 1, "  creating new vertices\n");

//	create new vertices from marked vertices
	for(VertexIterator iter = mg.begin<Vertex>(oldTopLevel);
		iter != mg.end<Vertex>(oldTopLevel); ++iter)
	{
		if(!refinement_is_allowed(*iter))
			continue;
			
		Vertex* v = *iter;

	//	create a new vertex in the next layer.
		//GMGR_PROFILE(GMGR_Refine_CreatingVertices);
		Vertex* nVrt = *mg.create_by_cloning(v, v);

	//	allow refCallback to calculate a new position
		if(m_refCallback)
			m_refCallback->new_vertex(nVrt, v);
		//GMGR_PROFILE_END();
	}


	UG_DLOG(LIB_GRID, 1, "  creating new edges\n");

//	create new vertices and edges from marked edges
	for(EdgeIterator iter = mg.begin<Edge>(oldTopLevel);
		iter != mg.end<Edge>(oldTopLevel); ++iter)
	{
		if(!refinement_is_allowed(*iter))
			continue;

	//	collect_objects_for_refine removed all edges that already were
	//	refined. No need to check that again.
		Edge* e = *iter;

	//	debug: make sure that both vertices may be refined
/*		#ifdef UG_DEBUG
			if(!refinement_is_allowed(e->vertex(0))
				|| !refinement_is_allowed(e->vertex(1)))
			{
				UG_LOG("Can't refine edge between vertices ");
				if(mg.has_vertex_attachment(aPosition)){
					Grid::VertexAttachmentAccessor<APosition> aaPos(mg, aPosition);
					UG_LOG(aaPos[e->vertex(0)] << " and " << aaPos[e->vertex(1)] << endl);
				}
				else if(mg.has_vertex_attachment(aPosition2)){
					Grid::VertexAttachmentAccessor<APosition2> aaPos(mg, aPosition2);
					UG_LOG(aaPos[e->vertex(0)] << " and " << aaPos[e->vertex(1)] << endl);
				}
				else if(mg.has_vertex_attachment(aPosition1)){
					Grid::VertexAttachmentAccessor<APosition1> aaPos(mg, aPosition1);
					UG_LOG(aaPos[e->vertex(0)] << " and " << aaPos[e->vertex(1)] << endl);
				}
			}
		#endif // UG_DEBUG
*/

		assert(refinement_is_allowed(e->vertex(0))
				&& refinement_is_allowed(e->vertex(1)));

		//GMGR_PROFILE(GMGR_Refine_CreatingEdgeVertices);
	//	create two new edges by edge-split
		RegularVertex* nVrt = *mg.create<RegularVertex>(e);

	//	allow refCallback to calculate a new position
		if(m_refCallback)
			m_refCallback->new_vertex(nVrt, e);
		//GMGR_PROFILE_END();

	//	split the edge
		//GMGR_PROFILE(GMGR_Refine_CreatingEdges);
		Vertex* substituteVrts[2];
		substituteVrts[0] = mg.get_child_vertex(e->vertex(0));
		substituteVrts[1] = mg.get_child_vertex(e->vertex(1));

		e->refine(vEdges, nVrt, substituteVrts);
		assert((vEdges.size() == 2) && "RegularEdge refine produced wrong number of edges.");
		mg.register_element(vEdges[0], e);
		mg.register_element(vEdges[1], e);
		//GMGR_PROFILE_END();
	}


	UG_DLOG(LIB_GRID, 1, "  creating new faces\n");

//	create new vertices and faces from marked faces
	for(FaceIterator iter = mg.begin<Face>(oldTopLevel);
		iter != mg.end<Face>(oldTopLevel); ++iter)
	{
		if(!refinement_is_allowed(*iter))
			continue;
			
		Face* f = *iter;
	//	collect child-vertices
		vVrts.clear();
		for(uint j = 0; j < f->num_vertices(); ++j)
			vVrts.push_back(mg.get_child_vertex(f->vertex(j)));

	//	collect the associated edges
		vEdgeVrts.clear();
		//bool bIrregular = false;
		for(uint j = 0; j < f->num_edges(); ++j)
			vEdgeVrts.push_back(mg.get_child_vertex(mg.get_edge(f, j)));

		//GMGR_PROFILE(GMGR_Refine_CreatingFaces);
		Vertex* newVrt;
		if(f->refine(vFaces, &newVrt, &vEdgeVrts.front(), NULL, &vVrts.front())){
		//	if a new vertex was generated, we have to register it
			if(newVrt){
				//GMGR_PROFILE(GMGR_Refine_CreatingVertices);
				mg.register_element(newVrt, f);
			//	allow refCallback to calculate a new position
				if(m_refCallback)
					m_refCallback->new_vertex(newVrt, f);
				//GMGR_PROFILE_END();
			}

		//	register the new faces and assign status
			for(size_t j = 0; j < vFaces.size(); ++j)
				mg.register_element(vFaces[j], f);
		}
		else{
			LOG("  WARNING in Refine: could not refine face.\n");
		}
		//GMGR_PROFILE_END();
	}


	UG_DLOG(LIB_GRID, 1, "  creating new volumes\n");

//	only used for tetrahedron or octahedron refinement
	vector<vector3> corners(6, vector3(0, 0, 0));

//	create new vertices and volumes from marked volumes
	for(VolumeIterator iter = mg.begin<Volume>(oldTopLevel);
		iter != mg.end<Volume>(oldTopLevel); ++iter)
	{
		if(!refinement_is_allowed(*iter))
			continue;

		Volume* v = *iter;
		//GMGR_PROFILE(GMGR_Refining_Volume);

	//	collect child-vertices
		//GMGR_PROFILE(GMGR_CollectingVolumeVertices);
		vVrts.clear();
		for(uint j = 0; j < v->num_vertices(); ++j)
			vVrts.push_back(mg.get_child_vertex(v->vertex(j)));
		//GMGR_PROFILE_END();

	//	collect the associated edges
		vEdgeVrts.clear();
		//GMGR_PROFILE(GMGR_CollectingVolumeEdgeVertices);
		//bool bIrregular = false;
		for(uint j = 0; j < v->num_edges(); ++j)
			vEdgeVrts.push_back(mg.get_child_vertex(mg.get_edge(v, j)));
		//GMGR_PROFILE_END();

	//	collect associated face-vertices
		vFaceVrts.clear();
		//GMGR_PROFILE(GMGR_CollectingVolumeFaceVertices);
		for(uint j = 0; j < v->num_faces(); ++j)
			vFaceVrts.push_back(mg.get_child_vertex(mg.get_face(v, j)));
		//GMGR_PROFILE_END();

	//	if we're performing tetrahedral or octahedral refinement, we have to collect
	//	the corner coordinates, so that the refinement algorithm may choose
	//	the best interior diagonal.
		vector3* pCorners = NULL;
		if((v->num_vertices() == 4) && m_refCallback){
			for(size_t i = 0; i < 4; ++i){
				m_refCallback->current_pos(&corners[i].x(), v->vertex(i), 3);
			}
			pCorners = &corners.front();
		}
		if((v->reference_object_id() == ROID_OCTAHEDRON) && m_refCallback){
			for(size_t i = 0; i < 6; ++i){
				m_refCallback->current_pos(&corners[i].x(), v->vertex(i), 3);
			}
			pCorners = &corners.front();
		}

		Vertex* newVrt;
		if(v->refine(vVols, &newVrt, &vEdgeVrts.front(), &vFaceVrts.front(),
					NULL, RegularVertex(), &vVrts.front(), pCorners)){
		//	if a new vertex was generated, we have to register it
			if(newVrt){
				mg.register_element(newVrt, v);
			//	allow refCallback to calculate a new position
				if(m_refCallback)
					m_refCallback->new_vertex(newVrt, v);
			}

		//	register the new faces and assign status
			for(size_t j = 0; j < vVols.size(); ++j)
				mg.register_element(vVols[j], v);
		}
		else{
			LOG("  WARNING in Refine: could not refine volume.\n");
		}
		//GMGR_PROFILE_END();
	}

//	done - clean up
	if(!bHierarchicalInsertionWasEnabled)
		mg.enable_hierarchical_insertion(false);

//	notify derivates that refinement ends
	refinement_step_ends();
	
//	clear the refinement-callback if necessary
	if(localRefCallbackSet){
		delete m_refCallback;
		m_refCallback = NULL;
	}

	m_messageHub->post_message(GridMessage_Adaption(GMAT_GLOBAL_REFINEMENT_ENDS,
													mg.get_grid_objects(oldTopLevel)));

	UG_DLOG(LIB_GRID, 1, "  refinement done.");
}

bool GlobalMultiGridRefiner::save_marks_to_file(const char* filename)
{
	GMGR_PROFILE(GlobalMultiGridRefiner_save_marks_to_file);
	if(!m_pMG){
		UG_THROW("ERROR in GlobalMultiGridRefiner::save_marks_to_file: No grid assigned!");
	}

	MultiGrid& mg = *m_pMG;
	SubsetHandler sh(mg);

	AssignGridToSubset(mg, sh, 2);
	int lvl = mg.num_levels() - 1;
	sh.assign_subset(mg.begin<Vertex>(lvl), mg.end<Vertex>(lvl), 0);
	sh.assign_subset(mg.begin<Edge>(lvl), mg.end<Edge>(lvl), 0);
	sh.assign_subset(mg.begin<Face>(lvl), mg.end<Face>(lvl), 0);
	sh.assign_subset(mg.begin<Volume>(lvl), mg.end<Volume>(lvl), 0);
	
	sh.subset_info(0).name = "refine";
	sh.subset_info(1).name = "no marks";

	AssignSubsetColors(sh);

	return SaveGridToFile(mg, sh, filename);
}

}//	end of namespace
