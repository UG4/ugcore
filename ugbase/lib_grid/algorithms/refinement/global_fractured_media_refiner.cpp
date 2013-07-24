//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y12 m04 d25

#include <cassert>
#include "common/profiler/profiler.h"
#include "global_fractured_media_refiner.h"
#include "lib_grid/algorithms/algorithms.h"
#include "lib_grid/common_attachments.h"

//define PROFILE_GLOBAL_FRACTURED_DOMAIN_REFINER if you want to profile
//the refinement code.
//#define PROFILE_GLOBAL_FRACTURED_DOMAIN_REFINER
#ifdef PROFILE_GLOBAL_FRACTURED_DOMAIN_REFINER
	#define GFDR_PROFILE_FUNC()	PROFILE_FUNC_GROUP("grid")
	#define GFDR_PROFILE(name)	PROFILE_BEGIN_GROUP(name, "grid")
	#define GFDR_PROFILE_END()	PROFILE_END()
#else
	#define GFDR_PROFILE_FUNC()
	#define GFDR_PROFILE(name)
	#define GFDR_PROFILE_END()
#endif

using namespace std;

namespace ug
{

//template <class TAPosition>
//GlobalFracturedMediaRefiner<TAPosition>::
GlobalFracturedMediaRefiner::
GlobalFracturedMediaRefiner(IRefinementCallback* refCallback) :
	IRefiner(refCallback),
	m_pMG(NULL),
	m_pSH(NULL)
{
//	initialize m_aPos to the default position attachment
	//m_aPos = GetDefaultPositionAttachment<TAPosition>();
}


//template <class TAPosition>
//GlobalFracturedMediaRefiner<TAPosition>::
GlobalFracturedMediaRefiner::
GlobalFracturedMediaRefiner(MultiGrid& mg, IRefinementCallback* refCallback) :
	IRefiner(refCallback),
	m_pSH(NULL)
{
	m_pMG = NULL;
	assign_grid(mg);
//	initialize m_aPos to the default position attachment
	//m_aPos = GetDefaultPositionAttachment<TAPosition>();
}


//template <class TAPosition>
//GlobalFracturedMediaRefiner<TAPosition>::
GlobalFracturedMediaRefiner::
~GlobalFracturedMediaRefiner()
{
	if(m_pMG)
		m_pMG->unregister_observer(this);
}


//template <class TAPosition>
//void GlobalFracturedMediaRefiner<TAPosition>::
void GlobalFracturedMediaRefiner::
grid_to_be_destroyed(Grid* grid)
{
	m_pMG = NULL;
}


//template <class TAPosition>
//void GlobalFracturedMediaRefiner<TAPosition>::
void GlobalFracturedMediaRefiner::
assign_grid(MultiGrid& mg)
{
	assign_grid(&mg);
}


//template <class TAPosition>
//void GlobalFracturedMediaRefiner<TAPosition>::
void GlobalFracturedMediaRefiner::
assign_grid(MultiGrid* mg)
{
	if(m_pMG){
		m_pMG->unregister_observer(this);
		m_marker.assign_grid(NULL);
		m_pMG = NULL;
	}
	
	if(mg){
		m_pMG = mg;
		set_message_hub(mg->message_hub());
		m_pMG->register_observer(this, OT_GRID_OBSERVER);
		m_marker.assign_grid(mg);
	}
}

//template <class TAPosition>
//void GlobalFracturedMediaRefiner<TAPosition>::
void GlobalFracturedMediaRefiner::
mark_as_fracture(int subInd, bool isFracture)
{
	if(subInd < 0)
		return;

	while(subInd > (int)m_subsetIsFracture.size())
		m_subsetIsFracture.push_back(false);
	if(subInd == (int)m_subsetIsFracture.size())
		m_subsetIsFracture.push_back(isFracture);
	else
		m_subsetIsFracture[subInd] = isFracture;
}

//template <class TAPosition>
//bool GlobalFracturedMediaRefiner<TAPosition>::
bool GlobalFracturedMediaRefiner::
is_fracture(int subInd)
{
	if((subInd < 0) || (subInd >= (int)m_subsetIsFracture.size()))
		return false;
	return m_subsetIsFracture[subInd];
}

//template <class TAPosition>
//void GlobalFracturedMediaRefiner<TAPosition>::
void GlobalFracturedMediaRefiner::
perform_refinement()
{
	GFDR_PROFILE_FUNC();
	UG_DLOG(LIB_GRID, 1, "GlobalFracturedMediaRefiner\n");

//	get the grid
	if(!m_pMG)
		UG_THROW("No grid assigned!");

	MultiGrid& mg = *m_pMG;

//	adjust marks for refinement
	adjust_marks();

//	if a debug file was specified, we'll now save the marks to that file
	if(!m_adjustedMarksDebugFilename.empty())
		save_marks_to_file(m_adjustedMarksDebugFilename.c_str());

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
			LOG("WARNING in GlobalFracturedMediaRefiner::refine(): auto-enabling VOLOPT_AUTOGENERATE_FACES.\n");
			mg.enable_options(VOLOPT_AUTOGENERATE_FACES);
		}
	}
	
	if(mg.num_faces() > 0){
		if(!mg.option_is_enabled(FACEOPT_AUTOGENERATE_EDGES))
		{
			LOG("WARNING in GlobalFracturedMediaRefiner::refine(): auto-enabling FACEOPT_AUTOGENERATE_EDGES.\n");
			mg.enable_options(FACEOPT_AUTOGENERATE_EDGES);
		}
	}

//	the old top level
	int oldTopLevel = mg.num_levels() - 1;
	UG_DLOG(LIB_GRID, 1, "REFINER: reserving memory...");

//	reserve enough memory to speed up the algo
//	todo: give better estimate based on fracture-subsets...
	GFDR_PROFILE(GFDR_Reserve);
	{
		int l = oldTopLevel;

		GFDR_PROFILE(GFDR_ReserveVrtData);
		mg.reserve<VertexBase>(mg.num<VertexBase>() +
					+ mg.num<VertexBase>(l) + mg.num<EdgeBase>(l)
					+ mg.num<Quadrilateral>(l) + mg.num<Hexahedron>(l));
		GFDR_PROFILE_END();

		GFDR_PROFILE(GFDR_ReserveEdgeData);
		mg.reserve<EdgeBase>(mg.num<EdgeBase>()
					+ 2 * mg.num<EdgeBase>(l) + 3 * mg.num<Triangle>(l)
					+ 4 * mg.num<Quadrilateral>(l) + 3 * mg.num<Prism>(l)
					+ mg.num<Tetrahedron>(l)
					+ 4 * mg.num<Pyramid>(l) + 6 * mg.num<Hexahedron>(l));
		GFDR_PROFILE_END();

		GFDR_PROFILE(GFDR_ReserveFaceData);
		mg.reserve<Face>(mg.num<Face>()
					+ 4 * mg.num<Face>(l) + 10 * mg.num<Prism>(l)
					+ 8 * mg.num<Tetrahedron>(l)
					+ 9 * mg.num<Pyramid>(l) + 12 * mg.num<Hexahedron>(l));
		GFDR_PROFILE_END();

		GFDR_PROFILE(GFDR_ReserveVolData);
		mg.reserve<Volume>(mg.num<Volume>()
					+ 8 * mg.num<Tetrahedron>(l) + 8 * mg.num<Prism>(l)
					+ 6 * mg.num<Pyramid>(l) + 8 * mg.num<Hexahedron>(l));
		GFDR_PROFILE_END();
	}
	GFDR_PROFILE_END();
	UG_DLOG(LIB_GRID, 1, " done.\n");

	UG_DLOG(LIB_GRID, 1, " refinement begins.\n");
//	notify derivates that refinement begins
	refinement_step_begins();
	
//	cout << "num marked edges: " << m_selMarks.num<EdgeBase>() << endl;
//	cout << "num marked faces: " << m_selMarks.num<Face>() << endl;

//	we want to add new elements in a new layer.
	bool bHierarchicalInsertionWasEnabled = mg.hierarchical_insertion_enabled();
	if(!bHierarchicalInsertionWasEnabled)
		mg.enable_hierarchical_insertion(true);


//	some buffers
	vector<VertexBase*> vVrts;
	vector<VertexBase*> vEdgeVrts;
	vector<VertexBase*> vFaceVrts;
	vector<EdgeBase*>	vEdges;
	vector<Face*>		vFaces;
	vector<Volume*>		vVols;
	
//	some repeatedly used objects
	EdgeDescriptor ed;
	FaceDescriptor fd;
	VolumeDescriptor vd;

	UG_DLOG(LIB_GRID, 1, "  creating new vertices\n");

//	create new vertices from marked vertices
	for(VertexBaseIterator iter = mg.begin<VertexBase>(oldTopLevel);
		iter != mg.end<VertexBase>(oldTopLevel); ++iter)
	{
		if(!refinement_is_allowed(*iter) || !m_marker.is_marked(*iter))
			continue;
			
		VertexBase* v = *iter;

	//	create a new vertex in the next layer.
		//GFDR_PROFILE(GFDR_Refine_CreatingVertices);
		VertexBase* nVrt = *mg.create_by_cloning(v, v);

	//	allow refCallback to calculate a new position
		if(m_refCallback)
			m_refCallback->new_vertex(nVrt, v);
		//GFDR_PROFILE_END();
	}


	UG_DLOG(LIB_GRID, 1, "  creating new edges\n");

//	create new vertices and edges from marked edges
	for(EdgeBaseIterator iter = mg.begin<EdgeBase>(oldTopLevel);
		iter != mg.end<EdgeBase>(oldTopLevel); ++iter)
	{
		if(!refinement_is_allowed(*iter))
			continue;

		EdgeBase* e = *iter;

		if(!(refinement_is_allowed(e->vertex(0)) && refinement_is_allowed(e->vertex(1))))
		{
			UG_THROW("Refinement has to be allowed for both corners of an edge, for"
					" which refinement is allowed.");
		}

	//	associated vertices on next level.
		VertexBase* substituteVrts[2];
		substituteVrts[0] = mg.get_child_vertex(e->vertex(0));
		substituteVrts[1] = mg.get_child_vertex(e->vertex(1));

	//	if the face is not marked, we'll simply clone it to the next level
		if(!m_marker.is_marked(e)){
			ed.set_vertices(substituteVrts[0], substituteVrts[1]);
			mg.create_by_cloning(e, ed, e);
			continue;
		}

		//GFDR_PROFILE(GFDR_Refine_CreatingEdgeVertices);
	//	create two new edges by edge-split
		Vertex* nVrt = *mg.create<Vertex>(e);

	//	allow refCallback to calculate a new position
		if(m_refCallback)
			m_refCallback->new_vertex(nVrt, e);
		//GFDR_PROFILE_END();

	//	split the edge
		//GFDR_PROFILE(GFDR_Refine_CreatingEdges);
		e->refine(vEdges, nVrt, substituteVrts);
		assert((vEdges.size() == 2) && "Edge refine produced wrong number of edges.");
		mg.register_element(vEdges[0], e);
		mg.register_element(vEdges[1], e);
		//GFDR_PROFILE_END();
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

	//	if the face is not marked, we'll simply clone it to the next level
		if(!m_marker.is_marked(f)){
			fd.set_num_vertices(vVrts.size());
			for(size_t i = 0; i < vVrts.size(); ++i)
				fd.set_vertex(i, vVrts[i]);
			mg.create_by_cloning(f, fd, f);
			continue;
		}

	//	collect the child vertices of associated edges
		vEdgeVrts.clear();
		for(uint j = 0; j < f->num_edges(); ++j){
			vEdgeVrts.push_back(mg.get_child_vertex(mg.get_edge(f, j)));
		}

		//GFDR_PROFILE(GFDR_Refine_CreatingFaces);
		VertexBase* newVrt;
		if(f->refine(vFaces, &newVrt, &vEdgeVrts.front(), NULL, &vVrts.front())){
		//	if a new vertex was generated, we have to register it
			if(newVrt){
				//GFDR_PROFILE(GFDR_Refine_CreatingVertices);
				mg.register_element(newVrt, f);
			//	allow refCallback to calculate a new position
				if(m_refCallback)
					m_refCallback->new_vertex(newVrt, f);
				//GFDR_PROFILE_END();
			}

		//	register the new faces and assign status
			for(size_t j = 0; j < vFaces.size(); ++j)
				mg.register_element(vFaces[j], f);
		}
		else{
			LOG("  WARNING in Refine: could not refine face.\n");
		}
		//GFDR_PROFILE_END();
	}


	UG_DLOG(LIB_GRID, 1, "  creating new volumes\n");

//	only used for tetrahedron refinement
	vector<vector3> corners(4, vector3(0, 0, 0));

//	create new vertices and volumes from marked volumes
	for(VolumeIterator iter = mg.begin<Volume>(oldTopLevel);
		iter != mg.end<Volume>(oldTopLevel); ++iter)
	{
		if(!refinement_is_allowed(*iter) || !m_marker.is_marked(*iter))
			continue;

		Volume* v = *iter;
		//GFDR_PROFILE(GFDR_Refining_Volume);

	//	collect child-vertices
		//GFDR_PROFILE(GFDR_CollectingVolumeVertices);
		vVrts.clear();
		for(uint j = 0; j < v->num_vertices(); ++j)
			vVrts.push_back(mg.get_child_vertex(v->vertex(j)));
		//GFDR_PROFILE_END();

	//	collect the associated edges
		vEdgeVrts.clear();
		//GFDR_PROFILE(GFDR_CollectingVolumeEdgeVertices);
		//bool bIrregular = false;
		for(uint j = 0; j < v->num_edges(); ++j)
			vEdgeVrts.push_back(mg.get_child_vertex(mg.get_edge(v, j)));
		//GFDR_PROFILE_END();

	//	collect associated face-vertices
		vFaceVrts.clear();
		//GFDR_PROFILE(GFDR_CollectingVolumeFaceVertices);
		for(uint j = 0; j < v->num_faces(); ++j)
			vFaceVrts.push_back(mg.get_child_vertex(mg.get_face(v, j)));
		//GFDR_PROFILE_END();

	//	if we're performing tetrahedral refinement, we have to collect
	//	the corner coordinates, so that the refinement algorithm may choose
	//	the best interior diagonal.
		vector3* pCorners = NULL;
		if((v->num_vertices() == 4) && m_refCallback){
			for(size_t i = 0; i < 4; ++i){
				m_refCallback->current_pos(&corners[i].x(), v->vertex(i), 3);
			}
			pCorners = &corners.front();
		}

		VertexBase* newVrt;
		if(v->refine(vVols, &newVrt, &vEdgeVrts.front(), &vFaceVrts.front(),
					NULL, Vertex(), &vVrts.front(), pCorners)){
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
		//GFDR_PROFILE_END();
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

	UG_DLOG(LIB_GRID, 1, "  refinement done.");
}

/*
//template <class TAPosition>
template <class TElem>
//void GlobalFracturedMediaRefiner<TAPosition>::
void GlobalFracturedMediaRefiner::
assign_elem_and_side_marks()
{
	typedef typename MultiGrid::traits<TElem>::iterator	ElemIter;
	typedef typename TElem::side	Side;

	if(!m_pMG)
		UG_THROW("No grid assigned!");

	if(!m_pSH)
		UG_THROW("No subset handler assigned");

	if(m_pMG->num_levels() == 0)
		UG_THROW("At least one level has to exist in the associated multi-grid.");

//	we might need an attachment accessor for this one.
//	Grid::VertexAttachmentAccessor<TAPosition> aaPos;
//	if(!aaPos.access(*m_pMG, m_aPos))
//		UG_THROW("You have to specify a valid position attachment through 'set_position_attachment'");

	MultiGrid& mg = *m_pMG;

//	iterate over all elements in the top level and adjust their marks
	int topLvl = mg.num_levels() - 1;

//	we only clear side marks and element marks here, since only those will be set again.
	m_marker.clear();

//	collect associated sides and elements in this vectors
	std::vector<Side*>	sides, sidesCon;
	std::vector<TElem*>	elems;
	std::queue<TElem*>	unclassified;

	for(ElemIter iter = mg.begin<TElem>(topLvl);
		iter != mg.end<TElem>(topLvl); ++iter)
	{
		TElem* e = *iter;

		m_marker.mark(e);

	//	if the element is not a fracture element, we'll mark all sides for refinement
	//	then we're done
		if(!is_fracture_element(e)){
			CollectAssociated(sides, mg, e);
			for(size_t i = 0; i < sides.size(); ++i)
				m_marker.mark(sides[i]);
			continue;
		}


	//	find the side which is connected to a non-fracture element
		CollectAssociated(sides, mg, e);
		Side* outerSide = NULL;
		for(size_t i_side = 0; i_side < sides.size(); ++i_side){
			Side* s = sides[i_side];
			CollectAssociated(elems, mg, s);
			bool gotOne = false;
			for(size_t i_elem = 0; i_elem < elems.size(); ++i_elem){
				if(!is_fracture_element(elems[i_elem])){
					gotOne = true;
					break;
				}
			}

			if(gotOne){
				outerSide = s;
				break;
			}
		}

		if(outerSide){
		//	the side is an interface between the fracture and the surrounding medium.
			m_marker.mark(outerSide);

		//	we now have to also mark the opposing side in e. Get the descriptor
		//	of that side
			typename geometry_traits<Side>::GeneralDescriptor desc;
			if(e->get_opposing_side(outerSide, desc)){
				Side* opSide = mg.get_element(desc);
				if(opSide)
					m_marker.mark(opSide);
			}
			else{
			//	this should only happen if we're on a cap-element
			//	we'll push the element to the list of unclassified elements
				unclassified.push(e);
			}
		}
		else{
			UG_THROW("Invalid fracture structure. Each fracture element has to"
					" be connected to non-fracture element.");
		}
	}

//	we'll count the number of failed tries in a row, to avoid an infinite loop
	size_t numFailedTries = 0;
	while(!unclassified.empty()){
		TElem* e = unclassified.front();
		unclassified.pop();

	//	first check if two sides are already marked (e.g. during an earlier iteration)
	//	if this is the case, we'll ignore the element
		CollectAssociated(sides, mg, e);
		size_t numMarked = num_marked(sides);

		if(numMarked == 2){
		//	the additional side has been marked during iteration over an other element
			continue;
		}
		else if(numMarked != 1){
			UG_THROW("Unknown fracture topology encountered. "
					<< numMarked << " sides marked, but only 1 should be marked."
					" In element at " << GetGeometricObjectCenter(mg, e));
		}

	//	iterate over neighbors which lie in the fracture. If exactly one
	//	neighbor has only one marked side (!=markedSide), then the side
	//	connecting both will be marked.
	//	If no neighbor is found with exactly one marked side, then an error is raised.
	//	If more then one neighbors are found with exactly one marked side, then
	//	the element is pushed back to queue.
		Side* connectingSide = NULL;
		int elemsFound = 0;// set to true, if an element with exactly one marked side is found.
		for(size_t i_side = 0; i_side < sides.size(); ++i_side){
			Side* s = sides[i_side];
			if(m_marker.is_marked(s))
				continue;

		//	get the connected element of the side
			CollectAssociated(elems, mg, s);
			if(elems.size() != 2){
				UG_THROW("Invalid fracture topology detected. Aborting.");
			}

			TElem* eCon;
			if(elems[0] == e)	eCon = elems[1];
			else				eCon = elems[0];

		//	count the number of marked edges in eCon
			CollectAssociated(sidesCon, mg, eCon);

			if(num_marked(sidesCon) == 1){
				++elemsFound;
				connectingSide = s;
			}
		}

		if(elemsFound == 0){
			UG_THROW("Couldn't decide which inner edge to refine."
					" Please check your fractures topology at "
					<< GetGeometricObjectCenter(mg, e));
		}
		else if(elemsFound == 1){
		//	we found the side which we have to mark.
			m_marker.mark(connectingSide);
			numFailedTries = 0;
		}
		else{
		//	we have to push the element back to the queue
			unclassified.push(e);
			++numFailedTries;
		//	prevent an infinite loop
			if(numFailedTries >= unclassified.size()){
				UG_THROW("Couldn't identify fracture structure.");
			}
		}
	}

//	mark sides of marked side elements
	mark_sides_of_marked_top_level_elements<Side>();
}
*/


template <class TElem>
void GlobalFracturedMediaRefiner::
assign_elem_and_side_marks()
{
	typedef typename MultiGrid::traits<TElem>::iterator	ElemIter;
	typedef typename TElem::side	Side;
	typedef typename Side::side		SideOfSide;

	if(!m_pMG)
		UG_THROW("No grid assigned!");

	if(!m_pSH)
		UG_THROW("No subset handler assigned");

	if(m_pMG->num_levels() == 0)
		UG_THROW("At least one level has to exist in the associated multi-grid.");

	MultiGrid& mg = *m_pMG;

//	iterate over all elements in the top level and adjust their marks
	int topLvl = mg.num_levels() - 1;

//	we only clear side marks and element marks here, since only those will be set again.
	m_marker.clear();

//	we'll also use a temporary marked to mark fixed sides.
	BoolMarker	fixedMarker(mg);

//	collect associated sides and elements in this vectors

	std::vector<Side*>	sides, sidesCon;
	std::vector<SideOfSide*> sidesOfSides;	//only used in 3d for fixed marks
	std::vector<TElem*>	elems;
	std::vector<TElem*>	caps;// elements at fracture caps (rim-boundary)

//	first we'll mark all elements which do not lie in a fracture. We'll also mark
//	their sides.
	for(ElemIter iter = mg.begin<TElem>(topLvl);
		iter != mg.end<TElem>(topLvl); ++iter)
	{
		TElem* e = *iter;

		if(!refinement_is_allowed(e))
			continue;

		if(!is_fracture_element(e)){
			m_marker.mark(e);
			CollectAssociated(sides, mg, e);
			for(size_t i = 0; i < sides.size(); ++i){
				m_marker.mark(sides[i]);
			}
		}
	}


//	call a callback, which allows to communicate marks in a parallel environment.
//	the default implementation does nothing
	communicate_marks(m_marker);

//	now iterate over all fracture elements and mark interior sides accordingly
	for(ElemIter iter = mg.begin<TElem>(topLvl);
		iter != mg.end<TElem>(topLvl); ++iter)
	{
		TElem* e = *iter;

		if(!refinement_is_allowed(e))
			continue;

		if(!is_fracture_element(e))
			continue;

		m_marker.mark(e);

	//	find the side which is already marked. If more than one is marked,
	//	we'll ignore the element. If no side is marked, we'll throw an error,
	//	since the fracture does not feature the required topology in this case.
		CollectAssociated(sides, mg, e);
		Side* markedSide = NULL;
		int numMarked = 0;
		for(size_t i_side = 0; i_side < sides.size(); ++i_side){
			Side* s = sides[i_side];
			if(m_marker.is_marked(s)){
				markedSide = s;
				++numMarked;
			}
		}

		if((numMarked == 0) || (numMarked > 2)){
			UG_THROW("Bad fracture topology encountered in element with center "
					 << GetGeometricObjectCenter(mg, e) << "! Aborting.");
		}

		typename geometry_traits<Side>::GeneralDescriptor desc;
		if(e->get_opposing_side(markedSide, desc)){
			Side* opSide = mg.get_element(desc);
			if(opSide){
				if((numMarked == 2) && (!m_marker.is_marked(opSide))){
					UG_THROW("Bad fracture topology encountered in element with center "
							 << GetGeometricObjectCenter(mg, e) << "! Aborting.");
				}
				m_marker.mark(opSide);
			}
			else{
				UG_THROW("Opposing side could not be found in grid. Check grid options.");
			}

		//	we'll iterate over the sides of the element again and mark currently
		//	unmarked sides as fixed sides
			for(size_t i_side = 0; i_side < sides.size(); ++i_side){
				Side* s = sides[i_side];
				if(!m_marker.is_marked(s))
				//	in 3d we'll also have to mark sides of sides, to guarantee,
				//	that fixed marks are communicated correctly.
					fixedMarker.mark(s);

					if(Side::dim == 2){
						CollectAssociated(sidesOfSides, mg, s);
						for(size_t i = 0; i < sidesOfSides.size(); ++i)
							fixedMarker.mark(sidesOfSides[i]);
					}
			}

		}
		else{
		//	this should only happen if we're on a cap-element. Store the element for later use.
			caps.push_back(e);
		}
	}

//	communicate fixed marks. Default implementation does nothing.
//	This can be used by a derived class to communicate marks to other processes.
	communicate_marks(fixedMarker);

//	now adjust marks at cap elements
	for(typename std::vector<TElem*>::iterator iter = caps.begin();
		iter != caps.end(); ++iter)
	{
		TElem* e = *iter;
		m_marker.mark(e);

	//	sides which are marked as fixed or have a fixed side them selfes (3d only)
	//	may not be refined. All others may be refined.
		CollectAssociated(sides, mg, e);
		for(size_t i_side = 0; i_side < sides.size(); ++i_side){
			bool refineSide = true;
			Side* s = sides[i_side];

			if(fixedMarker.is_marked(s)){
				refineSide = false;
				break;
			}
			else{
				if(Side::dim == 2){
					CollectAssociated(sidesOfSides, mg, s);
					for(size_t i = 0; i < sidesOfSides.size(); ++i){
						if(fixedMarker.is_marked(sidesOfSides[i])){
							refineSide = false;
							break;
						}
					}

					if(!refineSide)
						break;
				}
			}

		//	if the side shall be refined, we'll mark it
			if(refineSide)
				m_marker.mark(s);
		}

	}

//	mark sides of marked side elements
	if(Side::HAS_SIDES)
		mark_sides_of_marked_top_level_elements<Side>();
}


template <class TElem>
void GlobalFracturedMediaRefiner::
mark_sides_of_marked_top_level_elements()
{
	typedef typename MultiGrid::traits<TElem>::iterator	ElemIter;
	typedef typename TElem::side	Side;
	if(!m_pMG)
		UG_THROW("No grid assigned!");

	if(m_pMG->num_levels() == 0)
		UG_THROW("At least one level has to exist in the associated multi-grid.");

	MultiGrid& mg = *m_pMG;

//	iterate over all elements in the top level and adjust their marks
	int topLvl = mg.num_levels() - 1;

//	collect sides in this container
	vector<Side*> sides;

	for(ElemIter iter = mg.begin<TElem>(topLvl);
		iter != mg.end<TElem>(topLvl); ++iter)
	{
		TElem* e = *iter;
		if(m_marker.is_marked(e)){
			CollectAssociated(sides, mg, e);
			for(size_t i = 0; i < sides.size(); ++i){
				m_marker.mark(sides[i]);
			}
		}
	}

	if(Side::HAS_SIDES)
		mark_sides_of_marked_top_level_elements<Side>();
}


//template <class TAPosition>
//void GlobalFracturedMediaRefiner<TAPosition>::
void GlobalFracturedMediaRefiner::
adjust_marks()
{
	if(!m_pMG){
		UG_THROW("No grid assigned!");
	}

//	call the actual assign method
	if(m_pMG->num<Volume>() > 0)
		assign_elem_and_side_marks<Volume>();
	else if(m_pMG->num<Face>() > 0)
		assign_elem_and_side_marks<Face>();
	else{
	//	simply mark everything
		m_marker.mark(m_pMG->vertices_begin(), m_pMG->vertices_end());
		m_marker.mark(m_pMG->edges_begin(), m_pMG->edges_end());
	}
}

//template <class TAPosition>
//bool GlobalFracturedMediaRefiner<TAPosition>::
bool GlobalFracturedMediaRefiner::
save_marks_to_file(const char* filename)
{
	if(!m_pMG){
		UG_THROW("ERROR in GlobalFracturedMediaRefiner::save_marks_to_file: No grid assigned!");
	}

	MultiGrid& mg = *m_pMG;
	SubsetHandler sh(mg);

	AssignGridToSubset(mg, sh, 2);
	int lvl = mg.num_levels() - 1;

	for(VolumeIterator iter = mg.begin<Volume>(lvl); iter != mg.end<Volume>(lvl); ++iter){
		if(m_marker.is_marked(*iter)){
			if(is_fracture_element(*iter))
				sh.assign_subset(*iter, 1);
			else
				sh.assign_subset(*iter, 0);
		}
	}

	vector<EdgeBase*> edges;
	for(FaceIterator iter = mg.begin<Face>(lvl); iter != mg.end<Face>(lvl); ++iter){
		Face* f = *iter;
		CollectAssociated(edges, mg, f);
		if(num_marked(edges) != f->num_sides())
			sh.assign_subset(f, 1);
		else
			sh.assign_subset(f, 0);
	}

	for(EdgeBaseIterator iter = mg.begin<EdgeBase>(lvl); iter != mg.end<EdgeBase>(lvl); ++iter){
		if(m_marker.is_marked(*iter))
			sh.assign_subset(*iter, 0);
	}

	for(VertexBaseIterator iter = mg.begin<VertexBase>(lvl); iter != mg.end<VertexBase>(lvl); ++iter){
		if(m_marker.is_marked(*iter))
			sh.assign_subset(*iter, 0);
	}

	sh.subset_info(0).name = "refine regular";
	sh.subset_info(1).name = "refine anisotropic";
	sh.subset_info(2).name = "no marks";

	AssignSubsetColors(sh);

	return SaveGridToFile(mg, sh, filename);
}

template <class TElem>
size_t GlobalFracturedMediaRefiner::
num_marked(const std::vector<TElem*>& elems) const
{
	size_t num = 0;
	for(size_t i = 0; i < elems.size(); ++i){
		if(m_marker.is_marked(elems[i]))
			++num;
	}
	return num;
}

}//	end of namespace
