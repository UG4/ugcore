/*
 * Copyright (c) 2009-2015:  G-CSC, Goethe University Frankfurt
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

#include <vector>
#include "hanging_node_refiner_2d_irn.h"
#include "lib_grid/algorithms/geom_obj_util/geom_obj_util.h"

using namespace std;

namespace ug
{

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//BEGIN Debug methods
/*
static void PrintEdgeCoords(Edge* e, Grid::VertexAttachmentAccessor<APosition>& aaPos)
{
	LOG("(" << aaPos[e->vertex(0)].x() << ", " << aaPos[e->vertex(0)].y() <<
			"), (" << aaPos[e->vertex(1)].x() << ", " << aaPos[e->vertex(1)].y() << ")");
}

static void PrintConstrainingEdgeInfo(ConstrainingEdge* ce, Grid::VertexAttachmentAccessor<APosition>& aaPos)
{
	LOG("coords of constrained edges:" << endl);
	for(EdgeIterator cdEIter = ce->constrained_edges_begin();
						cdEIter != ce->constrained_edges_end(); ++cdEIter)
	{
		PrintEdgeCoords(*cdEIter, aaPos);
		LOG(endl);
	}

	LOG("coords of constrained vertices:" << endl);
	for(VertexIterator vrtIter = ce->constrained_vertices_begin();
		vrtIter != ce->constrained_vertices_end(); ++vrtIter)
	{
		LOG("pos: (" << aaPos[*vrtIter].x() << ", " << aaPos[*vrtIter].y() << ")");
		HangingVertex* hv = dynamic_cast<HangingVertex*>(*vrtIter);
		if(hv)
		{
			LOG(" local coord: " << hv->get_local_coordinate_1() << endl);
		}
		else
		{
			LOG(" no hanging vertex." << endl);
		}
	}
}
*/
/*
static bool CheckHangingNodeDegree(Grid& grid, uint irregularityRule, Grid::VertexAttachmentAccessor<APosition>& aaPos)
{
	bool retVal = true;
//	check if a constraining edge contains more hanging-vertices than irrgularityRule
	for(ConstrainingEdgeIterator iter = grid.begin<ConstrainingEdge>();
		iter != grid.end<ConstrainingEdge>(); ++iter)
	{
		ConstrainingEdge* ce = *iter;
		if(ce->num_constrained_vertices() > irregularityRule)
		{
			LOG("irregular constraining edge detected:" << endl);
			PrintEdgeCoords(ce, aaPos);
			PrintConstrainingEdgeInfo(ce, aaPos);
			retVal = false;
		}
	}
	return retVal;
}
*/
//END Debug methods

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	temporary helper methods
//TODO: move this method to a better place.
///	returns the local coordinate that v has relative to e
static number GetLocalVertexCoordinate(Edge* e, Vertex* vrt)
{
	if(vrt == e->vertex(0))
		return 0;
	else if(vrt == e->vertex(1))
		return 1.0;
	else if(ConstrainedVertex::type_match(vrt))
	{
		ConstrainedVertex* hv = dynamic_cast<ConstrainedVertex*>(vrt);
		return hv->get_local_coordinate_1();
	}

	LOG("WARNING in GetLocalVertexCoordinate(...): RegularVertex does not lie on edge. Returning -1." << endl);
	return -1.0;
}
/*
static void CalculateCenter(vector3& vOut, const vector3& v1, const vector3& v2)
{
	vector3 v;
	VecAdd(v, v1, v2);
	VecScale(vOut, v, 0.5);
}
*/
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	implementation of HangingNodeRefiner
HangingNodeRefiner2D_IRN::
HangingNodeRefiner2D_IRN(SPRefinementProjector projector) :
	IRefiner(refCallback),
	m_pGrid(nullptr),
	m_aVertex(false),
	m_irregularityRule(1)
{
}

HangingNodeRefiner2D_IRN::
HangingNodeRefiner2D_IRN(Grid& grid, SPRefinementProjector projector) :
	IRefiner(refCallback),
	m_pGrid(nullptr),
	m_aVertex(false),
	m_irregularityRule(1)
{
	assign_grid(grid);
}

HangingNodeRefiner2D_IRN::HangingNodeRefiner2D_IRN(const HangingNodeRefiner2D_IRN& hnr) : m_aVertex(false),
										m_irregularityRule(1)
{
	assert(!"WARNING in HangingNodeRefiner: Copy-Constructor not yet implemented. Please use references!");
	throw(int(0));
}

HangingNodeRefiner2D_IRN::~HangingNodeRefiner2D_IRN()
{
	if(m_pGrid)
	{
		m_pGrid->unregister_observer(this);
		m_pGrid->detach_from_edges(m_aVertex);
		m_pGrid->detach_from_faces(m_aVertex);
	}
}

void HangingNodeRefiner2D_IRN::assign_grid(Grid& grid)
{
	set_grid(&grid);
}

void HangingNodeRefiner2D_IRN::
set_grid(Grid* grid)
{
	if(m_pGrid)
	{
		m_pGrid->unregister_observer(this);
		m_selMarkedElements.assign_grid(nullptr);
		m_selScheduledElements.assign_grid(nullptr);
		m_pGrid->detach_from_edges(m_aVertex);
		m_pGrid->detach_from_faces(m_aVertex);
		m_pGrid = nullptr;
	}
	
	if(grid){
		m_pGrid = grid;
		grid->register_observer(this, ObserverType::OT_GRID_OBSERVER);
		m_selMarkedElements.assign_grid(*grid);
		m_selMarkedElements.enable_autoselection(false);
		m_selMarkedElements.enable_selection_inheritance(false);

		m_selScheduledElements.assign_grid(*grid);
		m_selScheduledElements.enable_autoselection(false);
		m_selScheduledElements.enable_selection_inheritance(false);

		grid->attach_to_edges_dv(m_aVertex, nullptr, false);
		grid->attach_to_faces_dv(m_aVertex, nullptr, false);

		m_aaVertexEDGE.access(*grid, m_aVertex);
		m_aaVertexFACE.access(*grid, m_aVertex);
	}
}

void HangingNodeRefiner2D_IRN::grid_to_be_destroyed(Grid* grid)
{
	if(m_pGrid)
		set_grid(nullptr);
}

void HangingNodeRefiner2D_IRN::clear_marks()
{
	m_selMarkedElements.clear();
}

void HangingNodeRefiner2D_IRN::mark_for_refinement(Edge* e)
{
	assert(m_pGrid && "ERROR in HangingNodeRefiner::mark_for_refinement(...): No grid assigned.");

	mark(e);
}

void HangingNodeRefiner2D_IRN::mark_for_refinement(Face* f)
{
	assert(m_pGrid && "ERROR in HangingNodeRefiner::mark_for_refinement(...): No grid assigned.");

	mark(f);
}

void HangingNodeRefiner2D_IRN::mark_for_refinement(Volume* v)
{
	assert(m_pGrid && "ERROR in HangingNodeRefiner::mark_for_refinement(...): No grid assigned.");

	mark(v);
}

bool HangingNodeRefiner2D_IRN::
set_irregularity_rule(uint irregularityRule)
{
	if(irregularityRule == 0)
		return false;
	m_irregularityRule = irregularityRule;
	return true;
}

uint HangingNodeRefiner2D_IRN::
get_irregularity_rule()
{
	return m_irregularityRule;
}
		
void HangingNodeRefiner2D_IRN::refine()
{
	uint irregularityRule = m_irregularityRule;
	
	assert(m_pGrid && "ERROR in HangingNodeRefiner::refine(...): No grid assigned.");
	assert(irregularityRule > 0 && "ERROR in HangingNodeRefiner::refine(...): irregularityRule has to be > 0.");

	if(irregularityRule < 1)
		return;

	Grid& grid = *m_pGrid;

//	check if a refinement-callback is set.
//	if not, we'll automatically set one, if a position attachment is available
	bool localRefCallbackSet = false;
	if(!m_refCallback){
		if(grid.has_vertex_attachment(aPosition)){
			localRefCallbackSet = true;
			m_refCallback = new RefinementCallbackLinear<APosition>(grid, aPosition);
		}
		else if(grid.has_vertex_attachment(aPosition2)){
			localRefCallbackSet = true;
			m_refCallback = new RefinementCallbackLinear<APosition2>(grid, aPosition2);
		}
	}
	
//	check grid options.
	if(!grid.option_is_enabled(GRIDOPT_FULL_INTERCONNECTION))
	{
		LOG("WARNING in HangingNodeRefiner::refine(): grid option GRIDOPT_FULL_INTERCONNECTION auto-enabled." << endl);
		grid.enable_options(GRIDOPT_FULL_INTERCONNECTION);
	}

//	containers used for temporary results
	EdgeVec 	vEdges;
	FaceVec 	vFaces;
	VolumeVec	vVolumes;

//	this outer loop makes sure, that ConstrainingEdges that have been created during
//	refinement won't violate the irregularityRule
//	we have to enter it at least once (since we have to call collect_objects_for_refine())
	do
	{
	//	fills the queues with the elements that have to be refined.
		collect_objects_for_refine();

	//	call this virtual method to allow derivates to execute some code
	//	before actual refinement takes place.
		refinement_step_begins();
		
	//	some that are being created during the following loop have to be
	//	refined them self. Those edges are pushed to vDelayedRefineEdges
	//	and are refined in the next iteration.
		EdgeVec vDelayedRefineEdges;

	//	iterate while there are edges scheduled for refinement.
		while(!m_selScheduledElements.empty<Edge>())
		{
		//	get the first edge and remove it from the queue.
			Edge* e = *m_selScheduledElements.begin<Edge>();
			m_selScheduledElements.deselect(e);

		//	depending on the type of edge we will proceed now.
			if(ConstrainedEdge::type_match(e))
			{
			//	check if the associated constraining edge is already scheduled
			//	for refinement. If not we have to check if the irregularity-rule
			//	would be hurt by refinement. If so the constraining edge has to
			//	be refined itself. Add it to the queue in this case.
				ConstrainedEdge* constrainedEdge = dynamic_cast<ConstrainedEdge*>(e);
				ConstrainingEdge* constrainingEdge = dynamic_cast<ConstrainingEdge*>(
													constrainedEdge->get_constraining_object());
				if(constrainingEdge != nullptr)
				{
					if(!is_marked(constrainingEdge))
					{
						if(constrainingEdge->num_constrained_vertices() >= irregularityRule)
						{
						//	the edge has to be refined.
						//	since we have to make sure that all marked faces have been created,
						//	before the constrainingEdge will be refined (this is important!!!)
						//	the edge will be scheduled for refinement in the next iteration.
							mark(constrainingEdge);
							vDelayedRefineEdges.push_back(constrainingEdge);
						}
					}
				}

				refine_constrained_edge(constrainedEdge);
			}
			else if(ConstrainingEdge::type_match(e))
			{
				ConstrainingEdge* constrainingEdge = dynamic_cast<ConstrainingEdge*>(e);

			//	if not all associated elements of higher dimension are marked for refinement,
			//	the constraining edge has to be conserved. (if the constraining edge was marked
			//	for refinement because of too many hanging nodes, then all associated elements
			//	will be marked.
				bool bAllMarked = true;
			//	check faces
				{
					CollectFaces(vFaces, grid, e);
					for(uint i = 0; i < vFaces.size(); ++i)
					{
						if(!is_marked(vFaces[i]))
						{
							LOG("found unmarked face: ");
							if(ConstrainingFace::type_match(vFaces[i]))
							{LOG("constraining face");}
							else if(ConstrainedFace::type_match(vFaces[i]))
								{LOG("constrained face");}
							else
								{LOG("normal face");}
							LOG(" with attachment index: " << grid.get_attachment_data_index(vFaces[i]) << endl);
							bAllMarked = false;
							break;
						}
					}
				}

			//	check volumes (only if bAllMarked is not already invalidated).
				if(bAllMarked && grid.num_volumes() > 0)
				{
					CollectVolumes(vVolumes, grid, e);
					for(uint i = 0; i < vVolumes.size(); ++i)
					{
						if(!is_marked(vVolumes[i]))
						{
						LOG("found unmarked volume\n");
							bAllMarked = false;
							break;
						}
					}
				}

			//	if all elements are marked, then we will perform a regular split.
				if(bAllMarked)
				{
				//	In order to correctly refine a constraining edge, we have to create
				//	two new edges connecting the constraining edges end-points with its
				//	central hanging node. The type of the newly created edges depends
				//	on the number of hanging nodes in the upper and lower part of the
				//	constraining edge.

				//	refine the edge
					bool scheduledEdgeReplaced[2];
					Edge* newEdges[2];
					refine_constraining_edge(constrainingEdge, &newEdges[0], &newEdges[1],
											scheduledEdgeReplaced[0], scheduledEdgeReplaced[1]);

				//	check if a newly created edge is a constraining edge with too many hanging nodes.
				//	check as well, whether the new edge replaced an old marked one.
					for(uint i = 0; i < 2; ++i)
					{
						if(scheduledEdgeReplaced[i])
						{
						//	this should be obsolete, since edges are refined in a fixed order - constraining edges are the last.
							assert(!"this piece of code shouldn't be called - at least if the programmer (yeah - thats me...) was correct.");

						//	a marked edge has been replaced. We have to schedule the new edge for refinement.
						//	this should never happen, since such a replaced edge would normally be a constrained-edge.
						//	constrained edges are however refined before constraining edges - due to the storage in
						//	the selector.
							//mark(newEdges[i]);
							//m_selScheduledElements.select(newEdges[i]);
						}
						ConstrainingEdge* nCE = dynamic_cast<ConstrainingEdge*>(newEdges[i]);
						if(nCE)
						{
						//	check how many hanging vertices lie on nCE.
						//	if there are more than allowed by the irregularityRule
						//	then schedule the edge for further refinement.
						//	The edge may not be refined until the adjacent faces have been refined.
						//	it has thus to be pushed to vDelayedRefineEdges.
							if(nCE->num_constrained_vertices() > irregularityRule)
							{
							//	the edge has to be refined - but not until the adjacent faces have been refined.
								mark(nCE);
								vDelayedRefineEdges.push_back(nCE);
							}
						}
					}
				}
				else
				{
					LOG("NOT ALL NEIGHBOURS MARKED... NO PROBLEM?!?\n");

				//	if not all neighbours are marked, the edge will not be splitted.
				//	find the center and store it with the edge for later use.
					Vertex* centerVrt = nullptr;
					for(VertexIterator vrtIter = constrainingEdge->constrained_vertices_begin();
						vrtIter != constrainingEdge->constrained_vertices_end(); ++vrtIter)
					{
						ConstrainedVertex* hv = dynamic_cast<ConstrainedVertex*>(*vrtIter);
						if(hv)
						{
							if(hv->get_local_coordinate_1() == 0.5)
							{
								centerVrt = hv;
								break;
							}
						}
					}

					if(centerVrt)
					{
						set_center_vertex(constrainingEdge, centerVrt);
					}
					else
					{
						assert(!"Program shouldn't reach this point!");
					}
				}
			}
			else
			{
			//	e is considered to be a standard edge.
			//	Depending on the neighborhood of e, we have to decide whether
			//	e shall be splitted into two normal edges, or if e has to be
			//	splitted by the insertion of a hanging vertex - thus transforming
			//	e to a constraining edge and creating two new constrained edges.
			//	check if all adjacent faces and volumes are selected.
				bool bAllSelected = true;

				CollectFaces(vFaces, grid, e);
				for(uint i = 0; i < vFaces.size(); ++i)
				{
					if(!is_marked(vFaces[i]))
					{
						bAllSelected = false;
						break;
					}
				}

			//TODO:	Volumes.

				if(bAllSelected)
				{
					refine_edge_with_normal_vertex(e);
				}
				else
				{
					refine_edge_with_hanging_vertex(e);
				}
			}
		}

	//	Edges are refined now.

	//	Refine Faces next.
	//	if there are no volumes, the face can be refined as normal
		if(grid.num_volumes() == 0)
		{
			while(!m_selScheduledElements.empty<Face>())
			{
				Face* f = *m_selScheduledElements.begin<Face>();
				m_selScheduledElements.deselect(f);

				refine_face_with_normal_vertex(f);
			}
		}
		else
		{
			while(!m_selScheduledElements.empty<Face>())
			{
				Face* f = *m_selScheduledElements.begin<Face>();
				m_selScheduledElements.deselect(f);

				//TODO:	distinguish face-types.
				if(ConstrainedFace::type_match(f))
				{
				}
				else if(ConstrainingFace::type_match(f))
				{
					ConstrainingFace* cf = dynamic_cast<ConstrainingFace*>(f);

				//	if not all associated elements of higher dimension are marked for refinement,
				//	the constraining face has to be conserved. (if the constraining face was marked
				//	for refinement because of too many hanging nodes, then all associated elements
				//	will be marked).
					bool bAllMarked = true;

				//	check volumes
					CollectVolumes(vVolumes, grid, f);
					for(uint i = 0; i < vVolumes.size(); ++i)
					{
						if(!is_marked(vVolumes[i]))
						{
						LOG("found unmarked volume\n");
							bAllMarked = false;
							break;
						}
					}

					if(bAllMarked)
					{
					//
					}
					else
					{
					//	the central hanging vertex has to be set as center-vertex for volume-refinement
						if(cf->num_constrained_vertices() > 0)
						{
						//TODO:	find the central vertex and call set_center_vertex with cf and the central vertex.
						//	don't delete cf. It has to remain since there are still volumes that reference cf.
							assert(!"quadrilaterals are not yet supported.");
						}
					}
				}
				else
				{
				//	f is considered to be a standard face.
				//	if all adjacent volumes are selected, then the face will be refined as normal
				//	if not, we have to create a constraining and some constrained faces.
					bool bAllSelected = true;
					CollectVolumes(vVolumes, grid, f);

					for(uint i = 0; i < vVolumes.size(); ++i)
					{
						if(!is_marked(vVolumes[i]))
						{
							bAllSelected = false;
							break;
						}
					}

					if(bAllSelected)
						refine_face_with_normal_vertex(f);
					else
					{
						refine_face_with_hanging_vertex(f);
					}
				}
			}
		}

	//	Faces are refined now.
	//	Refine Volumes next.
		while(!m_selScheduledElements.empty<Volume>())
		{
			Volume* v = *m_selScheduledElements.begin<Volume>();
			m_selScheduledElements.deselect(v);

			refine_volume_with_normal_vertex(v);
		}

	//	erase faces that are no longer needed.
		if(grid.num_volumes() > 0)
		{
			FaceIterator iter = m_selMarkedElements.begin<Face>();
			while(iter != m_selMarkedElements.end<Face>())
			{
				Face* f = *iter;
				++iter;
				CollectVolumes(vVolumes, grid, f);
				if(vVolumes.size() == 0)
				{
				//	erase
					grid.erase(f);
					LOG("erased face\n");
				}
			}
		}

	//	erase edges that are no longer needed.
		{
			EdgeIterator iter = m_selMarkedElements.begin<Edge>();
			while(iter != m_selMarkedElements.end<Edge>())
			{
				Edge* e = *iter;
				++iter;
				CollectFaces(vFaces, grid, e);
				if(vFaces.size() == 0)
				{
				//	erase
					grid.erase(e);
				}
			}
		}

	//	clear all marks
		clear_marks();

	//	mark the elements in the delayed-refinement containers
		for(uint i = 0; i < vDelayedRefineEdges.size(); ++i)
			mark_for_refinement(vDelayedRefineEdges[i]);

		vDelayedRefineEdges.clear();
		
	//	call this virtual method to allow derivates to execute some code
	//	after refinement took place.
		refinement_step_ends();
		
	}while(!m_selMarkedElements.empty());
	
//	clear the refinement-callback if necessary
	if(localRefCallbackSet){
		delete m_refCallback;
		m_refCallback = nullptr;
	}
}

void HangingNodeRefiner2D_IRN::collect_objects_for_refine()
{
//	Algorithm Layout:
//		The marks have to be adjusted. All edges on which a vertex has to be generated
//		should be marked.

//	make sure the queues are empty.
	m_selScheduledElements.clear();
	Grid& grid = *m_pGrid;

//	containers used for temporary results
	vector<Edge*> vEdges;
	vector<Face*> vFaces;
	vector<Volume*> vVolumes;


//	assert correct selection. see HangingVertexRefiner description.
//	select all faces and volumes adjacent to selected edges.
	{
		for(EdgeIterator iter = m_selMarkedElements.begin<Edge>();
			iter != m_selMarkedElements.end<Edge>(); ++iter)
		{
			m_selScheduledElements.select(*iter);
			CollectFaces(vFaces, grid, *iter);
			for(uint i = 0; i < vFaces.size(); ++i)
				mark_for_refinement(vFaces[i]);

			if(grid.num_volumes() > 0)
			{
				CollectVolumes(vVolumes, grid, *iter);
				for(uint i = 0; i < vVolumes.size(); ++i)
					mark_for_refinement(vVolumes[i]);
			}
		}
	}

//	now we'll select all edges adjacent to selected elements (faces and volumes).
//	these edges have to be refined later on.
	{
		for(FaceIterator iter = m_selMarkedElements.begin<Face>();
			iter != m_selMarkedElements.end<Face>(); ++iter)
		{
			m_selScheduledElements.select(*iter);
			CollectEdges(vEdges, grid, *iter);
			for(uint i = 0; i < vEdges.size(); ++i)
			{
				Edge* e = vEdges[i];
				if(!is_marked(e))
				{
					mark(e);
					m_selScheduledElements.select(e);
				}
			}

			if(grid.num_volumes() > 0)
			{
				CollectVolumes(vVolumes, grid, *iter);
				for(uint i = 0; i < vVolumes.size(); ++i)
					mark_for_refinement(vVolumes[i]);
			}
		}
	}

//	iterate over all volumes and mark adjacent edges and faces.
	{
		for(VolumeIterator iter = m_selMarkedElements.begin<Volume>();
			iter != m_selMarkedElements.end<Volume>(); ++iter)
		{
			m_selScheduledElements.select(*iter);
			CollectEdges(vEdges, grid, *iter);
			for(uint i = 0; i < vEdges.size(); ++i)
			{
				Edge* e = vEdges[i];
				if(!is_marked(e))
				{
					mark(e);
					m_selScheduledElements.select(e);
				}
			}

			CollectFaces(vFaces, grid, *iter);
			for(uint i = 0; i < vFaces.size(); ++i)
			{
				Face* f = vFaces[i];
				if(!is_marked(f))
				{
					mark(f);
					m_selScheduledElements.select(f);
				}
			}
		}
	}
/*
	LOG("num marked edges: " << m_selScheduledElements.num_selected<Edge>() << endl);
	LOG("num marked faces: " << m_selScheduledElements.num_selected<Face>() << endl);
	LOG("num marked volumes: " << m_selScheduledElements.num_selected<Volume>() << endl);
*/
//	we have to make sure that adjacent elements of selected constraining-edges
//	are all marked for refinement, too.
/*
	{
	//	first push all constraining edges to a queue
		queue<ConstrainingEdge*> qConstrainingEdges;
		for(ConstrainingEdgeIterator iter = m_selMarkedElements.begin<ConstrainingEdge>();
			iter != m_selMarkedElements.end<ConstrainingEdge>(); ++iter)
			qConstrainingEdges.push(*iter);

	//	do the same for constraining faces
		queue<ConstrainingFace*> qConstrainingFaces;
		for(ConstrainingFaceIterator iter = m_selMarkedElements.begin<ConstrainingFace>();
			iter != m_selMarkedElements.end<ConstrainingFace>(); ++iter)
			qConstrainingFaces.push(*iter);

	//	select connected faces and volumes together with their edges.
		while(!qConstrainingEdges.empty())
		{
			ConstrainingEdge* ce = qConstrainingEdges.front();
			qConstrainingEdges.pop();

			CollectFaces(vFaces, grid, ce);
			for(uint i = 0; i < vFaces.size(); ++i)
			{
				Face* f = vFaces[i];
				if(!is_marked(f))
				{
					mark(f);
					m_selScheduledElements.select(f);

				//	select the faces edges
					CollectEdges(vEdges, grid, f);
					for(uint j = 0; j < vEdges.size(); ++j)
					{
						Edge* e = vEdges[j];
						if(!is_marked(e))
						{
							mark(e);
							m_selScheduledElements.select(e);
							if(ConstrainingEdge::type_match(e))
								qConstrainingEdges.push((ConstrainingEdge*)e);
						}
					}
				//TODO: check whether f is a constraining face and push it to the queue if so.
				}
			}

			if(grid.num_volumes() > 0)
			{
				CollectVolumes(vVolumes, grid, ce);
				for(uint i = 0; i < vVolumes.size(); ++i)
				{
					Volume* v = vVolumes[i];
					mark(v);
					m_selScheduledElements.select(v);

				//	select the volumes edges
					CollectEdges(vEdges, grid, v);
					for(uint j = 0; j < vEdges.size(); ++j)
					{
						Edge* e = vEdges[j];
						if(!is_marked(e))
						{
							mark(e);
							m_selScheduledElements.select(e);
							if(ConstrainingEdge::type_match(e))
								qConstrainingEdges.push((ConstrainingEdge*)e);
						}
					}

				//	select the volumes faces
					CollectFaces(vFaces, grid, v);
					for(uint j = 0; j < vFaces.size(); ++j)
					{
						Face* f = vFaces[j];
						if(!is_marked(f))
						{
							mark(f);
							m_selScheduledElements.select(f);
						//TODO: if f is a constraining face, it has to be pushed to qConstrainingFaces.
						}
					}
				}
			}
		}
	}
*/
}

////////////////////////////////////////////////////////////////////////
//	implementation of refine-methods.
void HangingNodeRefiner2D_IRN::refine_constrained_edge(ConstrainedEdge* constrainedEdge)
{
	Grid& grid = *m_pGrid;

//	get the constraining object.
//TODO: add face support.
	ConstrainingEdge* constrainingEdge = dynamic_cast<ConstrainingEdge*>(
										constrainedEdge->get_constraining_object());
	if(constrainingEdge != nullptr)
	{
	//	a constrained edge is splitted into two constrained edges.
	//	the constraining element of ce will constrain the new edges, too.
	//	create the new vertex first (a hanging vertex)
		ConstrainedVertex* hv = *grid.create<ConstrainedVertex>(constrainedEdge);
	
	//	allow refCallback to calculate a new position
		if(m_refCallback)
			m_refCallback->new_vertex(hv, constrainedEdge);

		set_center_vertex(constrainedEdge, hv);

	//	calculate the local coordinates of hv.
		number localCoord = 0.5 * (GetLocalVertexCoordinate(constrainingEdge, constrainedEdge->vertex(0))
								+ GetLocalVertexCoordinate(constrainingEdge, constrainedEdge->vertex(1)));

		hv->set_local_coordinate_1(localCoord);

		constrainingEdge->add_constrained_object(hv);
		hv->set_constraining_object(constrainingEdge);

	//	unlink cde and ce
		constrainingEdge->unconstrain_object(constrainedEdge);
		constrainedEdge->set_constraining_object(nullptr);

	//	split the constrained edge and register new edges at the grid.
	//	constrainedEdge may not be deleted, since associated faces would
	//	be deleted too. Furthermore, hv is associated with it.
		vector<ConstrainedEdge*> vConstrainedEdges;
		if(constrainedEdge->refine(vConstrainedEdges, hv))
		{
			for(uint i = 0; i < vConstrainedEdges.size(); ++i)
			{
				ConstrainedEdge* nCE = vConstrainedEdges[i];
				grid.register_element(nCE, constrainedEdge);
				constrainingEdge->add_constrained_object(nCE);
				nCE->set_constraining_object(constrainingEdge);
			}
		}
		else
		{
		//	if this happens then there's a strange problem.
			assert(!"ERROR in RefineConstrainedEdge(...) - ConstrainedEdge refine failed due to unknown problem.");
			LOG("ERROR in RefineConstrainedEdge(...) - ConstrainedEdge refine failed due to unknown problem." << endl);
		}
	}
	else
	{
		assert(!"ERROR in RefineConstrainedEdge(...) - Constraining object of ConstrainedEdge should be a constraining edge!");
		LOG("ERROR in RefineConstrainedEdge(...) - Constraining object of ConstrainedEdge should be a constraining edge!" << endl);
	}
}

void HangingNodeRefiner2D_IRN::refine_constraining_edge(ConstrainingEdge* constrainingEdge,
												Edge** ppEdge1Out, Edge** ppEdge2Out,
												bool& scheduledEdgeReplaced1Out, bool& scheduledEdgeReplaced2Out)
{
	Grid& grid = *m_pGrid;
	scheduledEdgeReplaced1Out = scheduledEdgeReplaced2Out = false;

//	if a constraining edge is refined, its central hanging-vertex
//	will transform to a normal vertex.
//	Two new edges will be generated. Depending on the number and
//	distribution of hanging vertices on the constraining edge
//	this can be either normal edges or constraining edges again.
	ConstrainingEdge* ce = constrainingEdge;
	bool bHasHVrts[2] = {false, false};// indicates whether hanging nodes lie on the lower and/or the upper part.
	ConstrainedVertex* centerVrt = nullptr;

//	iterate through the constrained vertices of the edge and gather vertex-information
	for(VertexIterator hVrtIter = ce->constrained_vertices_begin();
		hVrtIter != ce->constrained_vertices_end(); ++hVrtIter)
	{
		if(ConstrainedVertex::type_match(*hVrtIter))
		{
			ConstrainedVertex* hv = dynamic_cast<ConstrainedVertex*>(*hVrtIter);
			if(hv->get_local_coordinate_1() < 0.5 - SMALL)
				bHasHVrts[0] = true;
			else if(hv->get_local_coordinate_1() > 0.5 + SMALL)
				bHasHVrts[1] = true;
			else
				centerVrt = hv;
		}
	}

//	we will only proceed if a centerVrt has been found
	if(centerVrt)
	{
	//	store resulting edges temporarily in this array
		Edge* newEdges[2] = {nullptr, nullptr};
		bool scheduledEdgeReplaced[2] = {false, false};

	//	centerVrt has to be transformed into a normal vertex.
	//	unlink it first from the constraining edge
		ce->unconstrain_object(centerVrt);
		vector3 vPos;

		RegularVertex* nCenterVrt = *grid.create_and_replace<RegularVertex>(centerVrt);

	//	store this vertex for later use on ce
		set_center_vertex(ce, nCenterVrt);

	//	create the new edges. The types of the new edges depend on bHasHVrts
	//	Handle both parts in the following loop.
		for(uint i = 0; i < 2; ++i)
		{
		//	set up information that depends on i
			number lowBorder = 0.5 * (number)i;
			number highBorder = lowBorder + 0.5;
			Vertex* vrt1, *vrt2;
			switch(i)
			{
				case 0:	vrt1 = ce->vertex(0);	vrt2 = nCenterVrt;	break;
				case 1:	vrt1 = nCenterVrt;	vrt2 = ce->vertex(1);	break;
				default: vrt1 = vrt2 = nullptr; break;//this will never happen!
			}

			if(bHasHVrts[i])
			{
			//	create a constraining edge
			//	check if an old edge already exits that connects those vertices
			//	if so, replace it.
				Edge* tEdge = grid.get_edge(vrt1, vrt2);
				ConstrainingEdge* nCE = nullptr;
				if(tEdge)
				{
				//	store the inserted vertex
					Vertex* tVrt = get_center_vertex(tEdge);
				//	if tEdge is constrained by ce, we have to unconstrain it first.
					if(ce->is_constrained_object(tEdge))
						ce->unconstrain_object(tEdge);
				//	create the new edge.
					nCE = *grid.create_and_replace<ConstrainingEdge>(tEdge);
					set_center_vertex(nCE, tVrt);
				}
				else
				{
					nCE = *grid.create<ConstrainingEdge>(EdgeDescriptor(vrt1, vrt2), ce);
				}

			//	store new edge
				newEdges[i] = nCE;

			//	copy hanging nodes
			//	only consider the ones that lie between lowBorder and highBorder.
				for(VertexIterator hVrtIter = ce->constrained_vertices_begin();
					hVrtIter != ce->constrained_vertices_end(); ++hVrtIter)
				{
					if(ConstrainedVertex::type_match(*hVrtIter))
					{
						ConstrainedVertex* hv = dynamic_cast<ConstrainedVertex*>(*hVrtIter);
						if((hv->get_local_coordinate_1() > lowBorder + SMALL) &&
							(hv->get_local_coordinate_1() < highBorder - SMALL))
						{
							hv->set_constraining_object(nCE);
							nCE->add_constrained_object(hv);
						}
					}
				}

			//	copy constrained edges
				EdgeIterator cdEIter = ce->constrained_edges_begin();
				while(cdEIter != ce->constrained_edges_end())
				{
					Edge* tEdge = *cdEIter;
					++cdEIter;
					if(ConstrainedEdge::type_match(tEdge))
					{
						ConstrainedEdge* cde = dynamic_cast<ConstrainedEdge*>(tEdge);
						bool bCopy = false;
						for(uint j = 0; j < 2; ++j)
						{
							Vertex* tVrt = cde->vertex(j);
						//	tVrt could be ce->vertex(0), ce->vertex(1), nCenterVrt or
						//	an arbitrary hanging node on ce.
						//	nCenterVrt won't help us to decide on which side cde lies...
							if(tVrt != nCenterVrt)
							{
								if((tVrt == vrt1) || (tVrt == vrt2))
								{
								//	since tVrt != nCenterVrt, a match means that cde is connected to
								//	ce's end-point of the correct side. (vrt1 or vrt2 was set to nCenterVrt)
									bCopy = true;
									break;
								}
								else if(ConstrainedVertex::type_match(tVrt) &&
										tVrt != ce->vertex(0) && tVrt != ce->vertex(1))
								{
								//	it was important to check whether tVrt is an end-point of ce, since those
								//	end-points can possibly be HangingVertices with local-coords again.
								//	in regard to ce however, their local coordinates should not be regarded.
									ConstrainedVertex* tmpHVrt = dynamic_cast<ConstrainedVertex*>(tVrt);
									if((tmpHVrt->get_local_coordinate_1() > lowBorder + SMALL) &&
										(tmpHVrt->get_local_coordinate_1() < highBorder - SMALL))
									{
										bCopy = true;
										break;
									}
								}
							}
						}

					//	if the edge belongs to this part of ce then copy it.
						if(bCopy)
						{
							cde->set_constraining_object(nCE);
							nCE->add_constrained_object(cde);
						}
					}
				}
			}
			else
			{
			//	no hanging nodes in the lower part of the edge.
			//	create a normal edge and replace the constrained edge of this part of ce
			//	find it first.
				Edge* replaceMe = nullptr;
				for(EdgeIterator cdEIter = ce->constrained_edges_begin();
					cdEIter != ce->constrained_edges_end(); ++cdEIter)
				{
					if(EdgeContains(*cdEIter, vrt1, vrt2))
					{
						replaceMe = *cdEIter;
						break;
					}
				}

				if(replaceMe)
				{
				//	unconstrain the edge
					ce->unconstrain_object(replaceMe);
				//	store the vertex that was placed on refineMe (if there was any)
					Vertex* tmpVrt = get_center_vertex(replaceMe);
				//	if the edge was scheduled we have to output this
					scheduledEdgeReplaced[i] = m_selScheduledElements.is_selected(replaceMe);
				//	create and replace
					RegularEdge* nEdge = *grid.create_and_replace<RegularEdge>(replaceMe);
				//	assign the stored vertex
					set_center_vertex(nEdge, tmpVrt);
				//	store the new edge
					newEdges[i] = nEdge;
				}
				else
				{
				//	if the program gets to this section, something went wrong. There should be an edge replaceMe,
				//	unless someone externally created a constraining edge without appropriate constrained edges.
					LOG("WARNING in PerformHangingNodeEdgeRefinement(...): Program shouldn't reach this point! Expect undefined behaviour." << endl);
					LOG("constraining refined " << ce << ": ");
					//PrintEdgeCoords(ce, m_aaPos);
					LOG(endl);

					//PrintConstrainingEdgeInfo(ce, m_aaPos);

					assert(!"WARNING in PerformHangingNodeEdgeRefinement(...): Program shouldn't reach this point! Expect undefined behaviour.");
					grid.create<RegularEdge>(EdgeDescriptor(vrt1, vrt2));
				}
			}
		}

	//	adjust local coordinates of hanging nodes
	//	do this here since problems will occur if done earlier.
		for(VertexIterator hVrtIter = ce->constrained_vertices_begin();
			hVrtIter != ce->constrained_vertices_end(); ++hVrtIter)
		{
			if(ConstrainedVertex::type_match(*hVrtIter))
			{
				ConstrainedVertex* hv = dynamic_cast<ConstrainedVertex*>(*hVrtIter);
				if(hv->get_local_coordinate_1() < 0.5 - SMALL)
					hv->set_local_coordinate_1(hv->get_local_coordinate_1() * (number)2);
				else if(hv->get_local_coordinate_1() > 0.5 + SMALL)
					hv->set_local_coordinate_1((hv->get_local_coordinate_1() - 0.5) * (number)2);
			}
		}

	//	clear all constrained edges.
		ce->clear_constrained_objects();

	//	prepare return values
		*ppEdge1Out = newEdges[0];
		*ppEdge2Out = newEdges[1];
		scheduledEdgeReplaced1Out = scheduledEdgeReplaced[0];
		scheduledEdgeReplaced2Out = scheduledEdgeReplaced[1];
	}
	else
	{
		LOG("ERROR in PerformHangingNodeEdgeRefinement(...): Program shouldn't reach this point. No central hanging-vertex was found on ConstrainingEdge." << endl);
		//LOG("dumping grid to DEBUG_GRID_SAVE.obj\n");
		//SaveGridToOBJ(grid, "DEBUG_GRID_SAVE.obj");
	//	save a debug grid. output the coordinates of the edges end-points.
		/*
		if(m_aaPos.valid())
		{
			LOG("coordinates of problematic edge: ");
			PrintEdgeCoords(ce, m_aaPos);
			LOG(endl);
		}
		*/
		LOG("RegularEdge attachment index: " << grid.get_attachment_data_index(ce) << endl);
		LOG("num_constrained_vertices: " << ce->num_constrained_vertices() << endl);
		LOG("coords of constrained vertices:" << endl);
		for(VertexIterator vrtIter = ce->constrained_vertices_begin();
			vrtIter != ce->constrained_vertices_end(); ++vrtIter)
		{
			//LOG("pos: (" << m_aaPos[*vrtIter].x() << ", " << m_aaPos[*vrtIter].y() << ")");
			ConstrainedVertex* hv = dynamic_cast<ConstrainedVertex*>(*vrtIter);
			if(hv)
			{
				LOG(" local coord: " << hv->get_local_coordinate_1() << endl);
			}
			else
			{
				LOG(" no hanging vertex." << endl);
			}
		}

		assert(!"ERROR in PerformHangingNodeEdgeRefinement(...): Program shouldn't reach this point. No central hanging-vertex was found on ConstrainingEdge.");
	}
}

void HangingNodeRefiner2D_IRN::refine_edge_with_normal_vertex(Edge* e)
{
	Grid& grid = *m_pGrid;

	RegularVertex* nVrt = *grid.create<RegularVertex>(e);
	set_center_vertex(e, nVrt);

//	allow refCallback to calculate a new position
	if(m_refCallback)
		m_refCallback->new_vertex(nVrt, e);

//	split the edge
	vector<Edge*> vEdges(2);
	e->refine(vEdges, nVrt);
	assert((vEdges.size() == 2) && "ERROR in RefineEdgeWithNormalVertex: Edge::refine - produced wrong number of edges.");
	grid.register_element(vEdges[0], e);
	grid.register_element(vEdges[1], e);

//	e may not be deleted at this point, since its associated triangles are still needed.
}

void HangingNodeRefiner2D_IRN::refine_edge_with_hanging_vertex(Edge* e)
{
	Grid& grid = *m_pGrid;
//	we have to insert a hanging node.
//	e has to be transformed to a constraining edge at the same time.
	assert(!ConstrainingEdge::type_match(e) && "ERROR in RefineEdgeWithHangingVertex(...): invalid operation. e is a constraining edge.");
	assert(!ConstrainedEdge::type_match(e) && "ERROR in RefineEdgeWithHangingVertex(...): invalid operation. e is a constrained edge.");

	ConstrainingEdge* ce = *grid.create_and_replace<ConstrainingEdge>(e);

	ConstrainedVertex* hv = *grid.create<ConstrainedVertex>(ce);

//	allow refCallback to calculate a new position
	if(m_refCallback)
		m_refCallback->new_vertex(hv, ce);

	set_center_vertex(ce, hv);
	hv->set_constraining_object(ce);
	ce->add_constrained_object(hv);

	hv->set_local_coordinate_1(0.5);

//	two constrained edges have to be created.
	ConstrainedEdge* nEdge[2];
	nEdge[0] = *grid.create<ConstrainedEdge>(EdgeDescriptor(ce->vertex(0), hv), ce);
	nEdge[1] = *grid.create<ConstrainedEdge>(EdgeDescriptor(hv, ce->vertex(1)), ce);

	for(uint i = 0; i < 2; ++i)
	{
		ce->add_constrained_object(nEdge[i]);
		nEdge[i]->set_constraining_object(ce);
	}
}

void HangingNodeRefiner2D_IRN::refine_face_with_normal_vertex(Face* f)
{
//LOG("refining face with normal vertex\n");
	Grid& grid = *m_pGrid;

	vector<Edge*> 	vEdges(f->num_edges());
	vector<Vertex*> vNewEdgeVertices(f->num_edges());
	vector<Face*>		vFaces(f->num_vertices());// heuristic
//	collect all associated edges.
	CollectEdges(vEdges, grid, f);
	size_t numEdges = vEdges.size();

	assert(numEdges == f->num_edges() && "ERROR in RefineFaceWithNormalVertex(...): associated edges missing.");

//	each should have an associated vertex. sort them into vNewEdgeVertices.
	for(size_t i = 0; i < numEdges; ++i)
	{
		Edge* e = vEdges[i];
		int edgeIndex = GetEdgeIndex(f, e);

		assert((edgeIndex >= 0) && (edgeIndex < (int)vEdges.size()) && "ERROR in RefineFaceWithNormalVertex(...): unknown problem in CollectEdges / GetEdgeIndex.");
		assert((get_center_vertex(e) != nullptr) && "ERROR in RefineFaceWithNormalVertex(...): no new vertex on refined edge.");
		vNewEdgeVertices[edgeIndex] = get_center_vertex(e);
	}

//	we'll perform a regular refine
	Vertex* vNewVrt = nullptr;
	f->refine_regular(vFaces, &vNewVrt, vNewEdgeVertices, nullptr, RegularVertex(), nullptr);

//	if a new vertex has been created during refine, then register it at the grid.
	if(vNewVrt)
	{
		grid.register_element(vNewVrt, f);

	//	allow refCallback to calculate a new position
		if(m_refCallback)
			m_refCallback->new_vertex(vNewVrt, f);
	}

	for(uint i = 0; i < vFaces.size(); ++i)
	{
		grid.register_element(vFaces[i], f);
	}

//	erase the old face if no selected volumes are in the grid.
//	otherwise it will be erased later on.
	if(grid.num_volumes() == 0)
		grid.erase(f);
}

void HangingNodeRefiner2D_IRN::refine_face_with_hanging_vertex(Face* f)
{
//LOG("refining face with hanging vertex\n");
	Grid& grid = *m_pGrid;

	vector<Edge*> 	vEdges(f->num_edges());
	vector<Vertex*> vNewEdgeVertices(f->num_edges());
	vector<Face*>		vFaces(f->num_vertices());// heuristic
//	collect all associated edges.
	CollectEdges(vEdges, grid, f);
	size_t numEdges = vEdges.size();

	assert(numEdges == f->num_edges() && "ERROR in RefineFaceWithNormalVertex(...): associated edges missing.");

//	each should have an associated vertex. sort them into vNewEdgeVertices.
	for(size_t i = 0; i < numEdges; ++i)
	{
		Edge* e = vEdges[i];
		int edgeIndex = GetEdgeIndex(f, e);

		assert((edgeIndex >= 0) && (edgeIndex < (int)vEdges.size()) && "ERROR in RefineFaceWithNormalVertex(...): unknown problem in CollectEdges / GetEdgeIndex.");
		assert((get_center_vertex(e) != nullptr) && "ERROR in RefineFaceWithNormalVertex(...): no new vertex on refined edge.");
		vNewEdgeVertices[edgeIndex] = get_center_vertex(e);
	}

//	check if a new vertex has to be created inside the face.
//	if so do it
	Vertex* vNewVrt = nullptr;
	if(f->num_vertices() > 3)
	{
		vNewVrt = *grid.create<ConstrainedVertex>(f);

	//	allow refCallback to calculate a new position
		if(m_refCallback)
			m_refCallback->new_vertex(vNewVrt, f);

	//TODO:	assign a local coordinate and store the vertex with the face
		assert(!"code missing: assign local coordinate.");
	}

//	the face has to be replaced by a constraining face.
//	we'll perform a switch here
	ConstrainingFace* constrainingFace = nullptr;
	uint faceType = f->reference_object_id();
	switch(faceType)
	{
		case ROID_TRIANGLE:
			{
				assert(f->num_vertices() == 3 && "bad number of vertices in face with ROID_TRIANGLE");
				//LOG(grid.get_attachment_data_index(f->vertex(0)) << endl);
			//	create the constraining triangle and replace the old face.
				constrainingFace = *grid.create_and_replace<ConstrainingTriangle>(f);
				//LOG(grid.get_attachment_data_index(constrainingFace->vertex(0)) << endl);

			//	create the constrained faces.
			//	the following triangle will not be registered at the grid. Just a temporary one.
				ConstrainedTriangle constrainedTri(constrainingFace->vertex(0),
													constrainingFace->vertex(1),
													constrainingFace->vertex(2));

			//	refine the constrainedTri and register new ones.
				Vertex* tmpVrt;
				constrainedTri.refine_regular(vFaces, &tmpVrt, vNewEdgeVertices, nullptr, RegularVertex(), nullptr);
			}
			break;
//TODO: add support for quadrilaterals
		default:
			assert(!"quadrilateral support missing.");
			break;
	}

	if(constrainingFace)
	{
	//	register the new faces
		for(uint i = 0; i < vFaces.size(); ++i)
		{
			ConstrainedFace* cdf = dynamic_cast<ConstrainedFace*>(vFaces[i]);
			assert(cdf && "constrained face refine did produce faces which are not constrained.");
			if(cdf)
			{
				grid.register_element(cdf, constrainingFace);
				cdf->set_constraining_object(constrainingFace);
				constrainingFace->add_constrained_object(cdf);
			}
		}

	//	we have to link the new constrained edges which have been auto-generated between the constrained faces.
	//	Since only edges that lie inside the constraining face are newly created, and since only those
	//	have to be linked with the constraining face, the following algorithm will be ok for
	//	triangles and quadrilaterals.
	//	Check for each new edge-vertex, if an edge exists with the new center vertex or with it's next neighbor.
		uint numNewEdgeVertices = vNewEdgeVertices.size();
		for(uint i = 0; i < numNewEdgeVertices; ++i)
		{
			if(vNewVrt)
			{
				ConstrainedEdge* e = dynamic_cast<ConstrainedEdge*>(grid.get_edge(vNewEdgeVertices[i], vNewVrt));
				if(e)
				{
				//	link e with the constraining face
					e->set_constraining_object(constrainingFace);
					constrainingFace->add_constrained_object(e);
				}
			}

		//	check if a constrained edge exists between the vertex and its next neighbor
			Vertex* vNext = vNewEdgeVertices[(i + 1) % numNewEdgeVertices];
			ConstrainedEdge* e = dynamic_cast<ConstrainedEdge*>(grid.get_edge(vNewEdgeVertices[i], vNext));
			if(e)
			{
			//	link e with the constraining face
				e->set_constraining_object(constrainingFace);
				constrainingFace->add_constrained_object(e);
			}
		}
	}
}

void HangingNodeRefiner2D_IRN::refine_constraining_face(std::vector<Face*>& vNewFacesOut,
												ConstrainingFace* constrainingFace)
{
//	associated edges of a constraining face are always constraining edges.
//	new faces are created.
//	we have to distribute the constrained faces to the new faces.
//	afterward we have to determine for each new face if it will
//	be a normal face or again by a constraining face.

//	should constrained faces have a coordinate by themselves?
//	this could be the coordinate of the center-vertex.
//	unsure how this could be computed.

//	other ways to compute the distribution: ???
}

void HangingNodeRefiner2D_IRN::refine_volume_with_normal_vertex(Volume* v)
{
	Grid& grid = *m_pGrid;

	vector<Edge*> 	vEdges(v->num_edges());
	vector<Vertex*> vNewEdgeVertices(v->num_edges());
	vector<Face*>		vFaces(v->num_faces());
	vector<Vertex*>	vNewFaceVertices(v->num_faces());
	vector<Volume*>		vVolumes(8);// heuristic
//	collect all associated edges.
	CollectEdges(vEdges, grid, v);
	uint numEdges = vEdges.size();

	assert(numEdges == v->num_edges() && "ERROR in RefineVolumeWithNormalVertex(...): associated edges missing.");

//	each should have an associated vertex. sort them into vNewEdgeVertices.
	for(uint i = 0; i < numEdges; ++i)
	{
		Edge* e = vEdges[i];
		int edgeIndex = GetEdgeIndex(v, e);

		assert((edgeIndex >= 0) && (edgeIndex < (int)vEdges.size()) && "ERROR in RefineVolumeWithNormalVertex(...): unknown problem in CollectEdges / GetEdgeIndex.");
		assert((get_center_vertex(e) != nullptr) && "ERROR in RefineVolumeWithNormalVertex(...): no new vertex on refined edge.");
		vNewEdgeVertices[edgeIndex] = get_center_vertex(e);
	}

//	get the vertices contained by adjacent faces
	CollectFaces(vFaces, grid, v);

	for(uint i = 0; i < vFaces.size(); ++i)
	{
	//TODO: assign vertices
		vNewFaceVertices[i] = nullptr;
	}

//TODO:	check if a new vertex has to be created inside the volume.
//	if so do it
/*
	Vertex* vNewVrt = nullptr;
	if(f->num_vertices() > 3)
	{
		vNewVrt = *grid.create<RegularVertex>(f);
	//	assign a new position
		if(m_aaPos.valid())
			m_aaPos[vNewVrt] = CalculateCenter(f, m_aaPos);
	}
*/
//	refine the volume and register new volumes at the grid.
	Vertex* createdVrt = nullptr;
	v->refine(vVolumes, &createdVrt, &vNewEdgeVertices.front(),
			  &vNewFaceVertices.front(), nullptr, RegularVertex(), nullptr);


	for(uint i = 0; i < vVolumes.size(); ++i)
	{
		grid.register_element(vVolumes[i], v);
	}

//	erase the old volume.
	grid.erase(v);
}

}//	end of namespace
