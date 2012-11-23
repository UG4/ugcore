// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y11 m01 d11

#include <vector>
#include "hanging_node_refiner_base.h"
#include "lib_grid/algorithms/geom_obj_util/geom_obj_util.h"
#include "lib_grid/algorithms/debug_util.h"
#include "lib_grid/algorithms/subset_util.h"
#include "lib_grid/file_io/file_io.h"
#include "ref_mark_adjusters/std_hnode_adjuster.h"

//define PROFILE_HANGING_NODE_REFINER if you want to profile
//the refinement code.
//#define PROFILE_HANGING_NODE_REFINER
#ifdef PROFILE_HANGING_NODE_REFINER
	#define HNODE_PROFILE_FUNC()	PROFILE_FUNC()
	#define HNODE_PROFILE_BEGIN(name)	PROFILE_BEGIN(name)
	#define HNODE_PROFILE_END()	PROFILE_END()
#else
	#define HNODE_PROFILE_FUNC()
	#define HNODE_PROFILE_BEGIN(name)
	#define HNODE_PROFILE_END()
#endif


using namespace std;

namespace ug{

HangingNodeRefinerBase::
HangingNodeRefinerBase(IRefinementCallback* refCallback) :
	IRefiner(refCallback),
	m_pGrid(NULL),
	m_nodeDependencyOrder1(true),
	m_automarkHigherDimensionalObjects(false),
	m_msgIdAdaption(-1)
{
	add_ref_mark_adjuster(StdHNodeAdjuster::create());
}

HangingNodeRefinerBase::
HangingNodeRefinerBase(const HangingNodeRefinerBase&)
{
	throw(UGError("no copy construction allowed."));
}

HangingNodeRefinerBase::
~HangingNodeRefinerBase()
{
	if(m_pGrid)
	{
		m_pGrid->unregister_observer(this);
	}
}

void HangingNodeRefinerBase::
set_grid(Grid* grid)
{
	if(m_pGrid)
	{
		m_pGrid->unregister_observer(this);
		m_selMarkedElements.assign_grid(NULL);
		m_pGrid = NULL;
	}

	if(grid){
		m_pGrid = grid;
		grid->register_observer(this, OT_GRID_OBSERVER);
		m_selMarkedElements.assign_grid(*grid);
		m_selMarkedElements.enable_autoselection(false);
		m_selMarkedElements.enable_selection_inheritance(false);
		m_msgIdAdaption = GridMessageId_Adaption(grid->message_hub());
		set_message_hub(grid->message_hub());
	}
}

void HangingNodeRefinerBase::grid_to_be_destroyed(Grid* grid)
{
	m_pGrid = NULL;
}

void HangingNodeRefinerBase::clear_marks()
{
	m_selMarkedElements.clear();
	m_newlyMarkedRefVrts.clear();
	m_newlyMarkedRefEdges.clear();
	m_newlyMarkedRefFaces.clear();
	m_newlyMarkedRefVols.clear();
}

bool HangingNodeRefinerBase::mark(VertexBase* v, RefinementMark refMark)
{
	assert(m_pGrid && "ERROR in HangingNodeRefinerBase::mark_for_refinement(...): No grid assigned.");
	if(get_mark(v) != refMark){
		if(refinement_is_allowed(v)){
			m_selMarkedElements.select(v, refMark);
			if(refMark & (RM_REGULAR | RM_ANISOTROPIC))
				m_newlyMarkedRefVrts.push_back(v);
			return true;
		}
	}
	return false;
}

bool HangingNodeRefinerBase::mark(EdgeBase* e, RefinementMark refMark)
{
	assert(m_pGrid && "ERROR in HangingNodeRefinerBase::mark_for_refinement(...): No grid assigned.");
	if(get_mark(e) != refMark){
		if(refinement_is_allowed(e)){
			m_selMarkedElements.select(e, refMark);
			if(refMark & (RM_REGULAR | RM_ANISOTROPIC))
				m_newlyMarkedRefEdges.push_back(e);
			return true;
		}
	}
	return false;
}

bool HangingNodeRefinerBase::mark(Face* f, RefinementMark refMark)
{
	assert(m_pGrid && "ERROR in HangingNodeRefinerBase::mark_for_refinement(...): No grid assigned.");

	if(get_mark(f) != refMark){
		if(refinement_is_allowed(f)){
			m_selMarkedElements.select(f, refMark);
			if(refMark & (RM_REGULAR | RM_ANISOTROPIC))
				m_newlyMarkedRefFaces.push_back(f);
			return true;
		}
	}
	return false;
}

bool HangingNodeRefinerBase::mark(Volume* v, RefinementMark refMark)
{
	assert(m_pGrid && "ERROR in HangingNodeRefinerBase::mark_for_refinement(...): No grid assigned.");
	if(get_mark(v) != refMark){
		if(refinement_is_allowed(v)){
			m_selMarkedElements.select(v, refMark);
			if(refMark & (RM_REGULAR | RM_ANISOTROPIC))
				m_newlyMarkedRefVols.push_back(v);
			return true;
		}
	}
	return false;
}


RefinementMark HangingNodeRefinerBase::
get_mark(VertexBase* v)
{
	return (RefinementMark)(m_selMarkedElements.get_selection_status(v)
							& ~HNRM_CONSTRAINED);
}

RefinementMark HangingNodeRefinerBase::
get_mark(EdgeBase* e)
{
	return (RefinementMark)(m_selMarkedElements.get_selection_status(e)
							& ~HNRM_CONSTRAINED);
}

RefinementMark HangingNodeRefinerBase::
get_mark(Face* f)
{
	return (RefinementMark)(m_selMarkedElements.get_selection_status(f)
							& ~HNRM_CONSTRAINED);
}

RefinementMark HangingNodeRefinerBase::
get_mark(Volume* v)
{
	return (RefinementMark)(m_selMarkedElements.get_selection_status(v)
							& ~HNRM_CONSTRAINED);
}


bool HangingNodeRefinerBase::save_marks_to_file(const char* filename)
{
	UG_DLOG(LIB_GRID, 1, "  saving marks to file...\n");
	if(!m_pGrid){
		UG_THROW("ERROR in HangingNodeRefinerBase::save_marks_to_file: No grid assigned!");
	}

	Grid& g = *m_pGrid;
	SubsetHandler sh(g);

	AssignGridToSubset(g, sh, 3);

	Selector& sel = get_refmark_selector();

	for(VertexBaseIterator iter = g.vertices_begin(); iter != g.vertices_end(); ++iter){
		Selector::status_t status = sel.get_selection_status(*iter);
		switch(status){
			case RM_REGULAR: sh.assign_subset(*iter, 0); break;
			case RM_ANISOTROPIC: sh.assign_subset(*iter, 1); break;
			case RM_COARSEN: sh.assign_subset(*iter, 2); break;
		}
	}

	for(EdgeBaseIterator iter = g.edges_begin(); iter != g.edges_end(); ++iter){
		Selector::status_t status = sel.get_selection_status(*iter);
		switch(status){
			case RM_REGULAR: sh.assign_subset(*iter, 0); break;
			case RM_ANISOTROPIC: sh.assign_subset(*iter, 1); break;
			case RM_COARSEN: sh.assign_subset(*iter, 2); break;
		}
	}

	for(FaceIterator iter = g.faces_begin(); iter != g.faces_end(); ++iter){
		Selector::status_t status = sel.get_selection_status(*iter);
		switch(status){
			case RM_REGULAR: sh.assign_subset(*iter, 0); break;
			case RM_ANISOTROPIC: sh.assign_subset(*iter, 1); break;
			case RM_COARSEN: sh.assign_subset(*iter, 2); break;
		}
	}

	for(VolumeIterator iter = g.volumes_begin(); iter != g.volumes_end(); ++iter){
		Selector::status_t status = sel.get_selection_status(*iter);
		switch(status){
			case RM_REGULAR: sh.assign_subset(*iter, 0); break;
			case RM_ANISOTROPIC: sh.assign_subset(*iter, 1); break;
			case RM_COARSEN: sh.assign_subset(*iter, 2); break;
		}
	}

	sh.subset_info(0).name = "refine regular";
	sh.subset_info(1).name = "refine anisotropic";
	sh.subset_info(2).name = "coarsen";
	sh.subset_info(3).name = "no marks";

	AssignSubsetColors(sh);

	return SaveGridToFile(g, sh, filename);
}

void HangingNodeRefinerBase::
enable_node_dependency_order_1(bool bEnable)
{
	m_nodeDependencyOrder1 = bEnable;
	for(size_t i = 0; i < m_refMarkAdjusters.size(); ++i){
		m_refMarkAdjusters[i]->enable_node_dependency_order_1(bEnable);
	}
}

void HangingNodeRefinerBase::perform_refinement()
{
	HNODE_PROFILE_FUNC();
	UG_DLOG(LIB_GRID, 1, "performing hanging-node-refine:\n");

	if(!m_pGrid)
		throw(UGError("ERROR in HangingNodeRefinerBase::refine(...): No grid assigned."));

	if(m_selMarkedElements.grid() != m_pGrid)
		throw(UGError("selector not initialized properly. Use HangingNodeRefinerBase::set_grid."));

	Grid& grid = *m_pGrid;

//	notify the grid's message hub that refinement begins
	grid.message_hub()->post_message(m_msgIdAdaption,
							GridMessage_Adaption(GMAT_HNODE_REFINEMENT_BEGINS));

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
		else if(grid.has_vertex_attachment(aPosition1)){
			localRefCallbackSet = true;
			m_refCallback = new RefinementCallbackLinear<APosition1>(grid, aPosition1);
		}
	}

//	check grid options.
	if(!grid.option_is_enabled(GRIDOPT_AUTOGENERATE_SIDES))
	{
		LOG("WARNING in HangingNodeRefiner::refine(): grid option GRIDOPT_AUTOGENERATE_SIDES auto-enabled." << endl);
		grid.enable_options(GRIDOPT_AUTOGENERATE_SIDES);
	}

//	containers used for temporary results
	vector<EdgeBase*> 	vEdges;
	vector<Face*>	 	vFaces;
	vector<Volume*>		vVols;

//	fills the queues with the elements that have to be refined.
	collect_objects_for_refine();
//	assigns hnode marks
	assign_hnode_marks();

//	if a debug file was specified, we'll now save the marks to that file
	if(!m_adjustedMarksDebugFilename.empty())
		save_marks_to_file(m_adjustedMarksDebugFilename.c_str());

//	call pre_refine to allow derived classes to perform some actions
	pre_refine();

////////////////////////////////
//	ConstrainedEdges
//	constrained edges can be ignored, since they will be refined together with
//	their constraining edges.

	UG_DLOG(LIB_GRID, 1, "  constraining edges.\n");
////////////////////////////////
//	ConstrainingEdges
//	iterate through scheduled cg-edges
	{
		ConstrainingEdgeIterator iter = m_selMarkedElements.begin<ConstrainingEdge>();
		while(iter != m_selMarkedElements.end<ConstrainingEdge>()){
			ConstrainingEdge* cge = *iter;
			++iter;

			if(!marked_for_hnode_refinement(cge) && refinement_is_allowed(cge))
				refine_constraining_edge(cge);
		}
	}

	UG_DLOG(LIB_GRID, 1, "  normal edges.\n");
////////////////////////////////
//	Normal Edges
//	iterate through scheduled edges
	{
		EdgeIterator iter = m_selMarkedElements.begin<Edge>();
		while(iter != m_selMarkedElements.end<Edge>())
		{
			Edge* e = *iter;
			++iter;

			if(!refinement_is_allowed(e))
				continue;
/*
		//	check whether all connected elements are marked.
		//	If so, we can perform a normal refine.
		//	If not, we'll create a new hanging vertex and new constrained edges.
		//	The edge itself will be transformed to a constraining edge.
			bool bAllMarked = true;

			if(grid.num_faces() > 0){
				CollectFaces(vFaces, grid, e);
				for(size_t i = 0; i < vFaces.size(); ++i){
					if(!is_marked(vFaces[i])){
						bAllMarked = false;
						break;
					}
				}
			}

		//	volumes only have to be checked if all faces were marked.
			if(bAllMarked && (grid.num_volumes() > 0)){
				CollectVolumes(vVols, grid, e);
				for(size_t i = 0; i < vVols.size(); ++i){
					if(!is_marked(vVols[i])){
						bAllMarked = false;
						break;
					}
				}
			}

			if(bAllMarked)
				refine_edge_with_normal_vertex(e);
			else
				refine_edge_with_hanging_vertex(e);
*/
			if(marked_for_hnode_refinement(e)){
				refine_edge_with_hanging_vertex(e);
			}
			else{
				refine_edge_with_normal_vertex(e);
			}
		}
	}


////////////////////////////////
//	Faces
//	constrained faces are ignored, since they are automatically refined during
//	constraining-face-refine.

	UG_DLOG(LIB_GRID, 1, "  constraining triangles.\n");
//	constraining triangles
	{
		ConstrainingTriangleIterator iter = m_selMarkedElements.begin<ConstrainingTriangle>();
		while(iter != m_selMarkedElements.end<ConstrainingTriangle>()){
			ConstrainingTriangle* cgf = *iter;
			++iter;
			if(!marked_for_hnode_refinement(cgf) && refinement_is_allowed(cgf))
				refine_constraining_face(cgf);
		}
	}

	UG_DLOG(LIB_GRID, 1, "  constraining quadrilaterals.\n");
//	constraining quadrilaterals
	{
		ConstrainingQuadrilateralIterator iter = m_selMarkedElements.begin<ConstrainingQuadrilateral>();
		while(iter != m_selMarkedElements.end<ConstrainingQuadrilateral>()){
			ConstrainingQuadrilateral* cgf = *iter;
			++iter;
			if(!marked_for_hnode_refinement(cgf) && refinement_is_allowed(cgf))
				refine_constraining_face(cgf);
		}
	}

	UG_DLOG(LIB_GRID, 1, "  normal triangles.\n");
//	normal triangles
	{
		TriangleIterator iter = m_selMarkedElements.begin<Triangle>();
		while(iter != m_selMarkedElements.end<Triangle>()){
			Face* f = *iter;
			++iter;

			if(!refinement_is_allowed(f))
				continue;
/*
		//	check whether all associated volumes are marked
			bool bAllMarked = true;
			if(grid.num_volumes() > 0){
				CollectVolumes(vVols, grid, f);
				for(size_t i = 0; i < vVols.size(); ++i){
					if(!is_marked(vVols[i])){
						bAllMarked = false;
						break;
					}
				}
			}

		//	if all are marked, the face can be normally refined.
		//	If not, we'll create constrained faces.
			if(bAllMarked)
				refine_face_with_normal_vertex(f);
			else{
				refine_face_with_hanging_vertex(f);
			}
*/
			if(marked_for_hnode_refinement(f))
				refine_face_with_hanging_vertex(f);
			else
				refine_face_with_normal_vertex(f);
		}
	}

	UG_DLOG(LIB_GRID, 1, "  normal quadrilaterals.\n");
//	normal quadrilaterals
	{
		QuadrilateralIterator iter = m_selMarkedElements.begin<Quadrilateral>();
		while(iter != m_selMarkedElements.end<Quadrilateral>()){
			Face* f = *iter;
			++iter;

			if(!refinement_is_allowed(f))
				continue;
/*
		//	check whether all associated volumes are marked
			bool bAllMarked = true;
			if(grid.num_volumes() > 0){
				CollectVolumes(vVols, grid, f);
				for(size_t i = 0; i < vVols.size(); ++i){
					if(!is_marked(vVols[i])){
						bAllMarked = false;
						break;
					}
				}
			}

		//	if all are marked, the face can be normally refined.
		//	If not, we'll create constrained faces.
			if(bAllMarked)
				refine_face_with_normal_vertex(f);
			else{
				refine_face_with_hanging_vertex(f);
			}
*/
			if(marked_for_hnode_refinement(f))
				refine_face_with_hanging_vertex(f);
			else
				refine_face_with_normal_vertex(f);

		}
	}


////////////////////////////////
//	Volumes
	UG_DLOG(LIB_GRID, 1, "  volumes.\n");
	{
		VolumeIterator iter = m_selMarkedElements.begin<Volume>();
		while(iter != m_selMarkedElements.end<Volume>()){
			Volume* vol = *iter;
			++iter;
			if(refinement_is_allowed(vol))
				refine_volume_with_normal_vertex(vol);
		}
	}

	UG_DLOG(LIB_GRID, 1, "  refinement done.\n");


////////////////////////////////
//	call post_refine to allow derived classes to perform some actions
	post_refine();


////////////////////////////////
//	Clean up
//	clear the refinement-callback if necessary
	if(localRefCallbackSet){
		delete m_refCallback;
		m_refCallback = NULL;
	}

	clear_marks();
	UG_DLOG(LIB_GRID, 1, "  done.\n");

//	notify the grid's message hub that refinement ends
	grid.message_hub()->post_message(m_msgIdAdaption,
							GridMessage_Adaption(GMAT_HNODE_REFINEMENT_ENDS));
}

template <class TElem>
bool HangingNodeRefinerBase::remove_coarsen_marks()
{
	typedef typename geometry_traits<TElem>::iterator	ElemIter;
	bool removedCoarsenMark = false;
	for(ElemIter iter = m_selMarkedElements.begin<TElem>();
		iter != m_selMarkedElements.end<TElem>();)
	{
		TElem* e = *iter;
		++iter;
		if(get_mark(e) == RM_COARSEN){
			m_selMarkedElements.deselect(e);
			removedCoarsenMark = true;
		}
	}

	return removedCoarsenMark;
}

void HangingNodeRefinerBase::collect_objects_for_refine()
{
	HNODE_PROFILE_FUNC();
	UG_DLOG(LIB_GRID, 1, "hnode_ref-start: collect_objects_for_refine\n");

//	build correct selection. see HangingVertexRefiner description.
//	unmark all elements which are marked for coarsening

	bool removedCoarseMarks = remove_coarsen_marks<VertexBase>();
	removedCoarseMarks |= remove_coarsen_marks<EdgeBase>();
	removedCoarseMarks |= remove_coarsen_marks<Face>();
	removedCoarseMarks |= remove_coarsen_marks<Volume>();

	if(removedCoarseMarks){
		UG_LOG("WARNING in HangingNodeRefinerBase::collect_objects_for_refine: "
				"Removed coarsen marks.\n");
	}

	std::vector<VertexBase*>	newlyMarkedVrts;
	std::vector<EdgeBase*>		newlyMarkedEdges;
	std::vector<Face*>			newlyMarkedFaces;
	std::vector<Volume*>		newlyMarkedVols;

	bool continueAdjustment = true;
	while(continueAdjustment){
	//	we don't simply pass m_newlyMarkedXXX to the adjusters, since we want to
	//	record newly marked elements during adjustment. Those newly marked elems
	//	are then used for the next calls, etc.
		newlyMarkedVrts.swap(m_newlyMarkedRefVrts);
		m_newlyMarkedRefVrts.clear();
		newlyMarkedEdges.swap(m_newlyMarkedRefEdges);
		m_newlyMarkedRefEdges.clear();
		newlyMarkedFaces.swap(m_newlyMarkedRefFaces);
		m_newlyMarkedRefFaces.clear();
		newlyMarkedVols.swap(m_newlyMarkedRefVols);
		m_newlyMarkedRefVols.clear();

	//	call the adjusters
		for(size_t i = 0; i < m_refMarkAdjusters.size(); ++i){
			m_refMarkAdjusters[i]->ref_marks_changed(*this, newlyMarkedVrts,
														newlyMarkedEdges,
														newlyMarkedFaces,
														newlyMarkedVols);
		}

		bool continueRequired =	(!m_newlyMarkedRefVrts.empty())
							 || (!m_newlyMarkedRefEdges.empty())
							 || (!m_newlyMarkedRefFaces.empty())
							 || (!m_newlyMarkedRefVols.empty());

		continueAdjustment = continue_collect_objects_for_refine(continueRequired);
	}

	UG_DLOG(LIB_GRID, 1, "hnode_ref-stop: collect_objects_for_refine\n");
}

void HangingNodeRefinerBase::
assign_hnode_marks()
{
	UG_DLOG(LIB_GRID, 1, "  assigning hnode marks...\n");
//	iterate over all faces and volumes. If the element is not marked, but
//	a side is marked, the side has to be marked for hnode refinement.
//	Note that we won't mark any new elements here - we only adjust the marks.
//	Note also that we won't remove any marks during this algorithm (neither normal
//	nor hnode marks).
	vector<Face*> faces;
	vector<Volume*> vols;

//	the grid
	UG_ASSERT(m_pGrid, "A grid is required to perform this operation!");
	Grid& grid = *m_pGrid;

	if(grid.num<Volume>() > 0){
		for(FaceIterator iter = m_selMarkedElements.begin<Face>();
			iter != m_selMarkedElements.end<Face>(); ++iter)
		{
			Face* f = *iter;
			CollectAssociated(vols, grid, f);
			for(size_t i = 0; i < vols.size(); ++i){
				if(refinement_is_allowed(vols[i])
				   && (!m_selMarkedElements.is_selected(vols[i])))
				{
					mark_for_hnode_refinement(f, true);
					break;
				}
			}
		}
	}
	
	if(grid.num<Face>() > 0){
		for(EdgeBaseIterator iter = m_selMarkedElements.begin<EdgeBase>();
			iter != m_selMarkedElements.end<EdgeBase>(); ++iter)
		{
			EdgeBase* e = *iter;
			CollectAssociated(faces, grid, e);
			for(size_t i = 0; i < faces.size(); ++i){
				if(marked_for_hnode_refinement(faces[i])){
					mark_for_hnode_refinement(e, true);
					break;
				}
				else if(refinement_is_allowed(faces[i])
						&& (!m_selMarkedElements.is_selected(faces[i])))
				{
					mark_for_hnode_refinement(e, true);
					break;
				}
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////
//	implementation of refine-methods.
void HangingNodeRefinerBase::
refine_constraining_edge(ConstrainingEdge* cge)
{
	HNODE_PROFILE_FUNC();

//	make sure that there is one hanging vertex and two constrained edges.
	assert(cge->num_constrained_vertices() == 1 && "bad number of constrained vertices. There has to be exactly 1.");
	assert(cge->num_constrained_edges() == 2 && "bad number of constrained edges. There have to be exactly 2.");

//	the grid
	Grid& grid = *m_pGrid;

//	the central hanging vertex has to be transformed into a normal vertex
	ConstrainedVertex* centralHV = NULL;
	if(cge->num_constrained_vertices() > 0)
		centralHV = dynamic_cast<ConstrainedVertex*>(cge->constrained_vertex(0));

	if(!centralHV){
		UG_LOG("The central hanging vertex of a constraining edge is missing. ignoring edge.\n");
		return;
	}

//	replace the central vertex with a normal vertex
	Vertex* centerVrt = *grid.create_and_replace<Vertex>(centralHV);

//	Iterate over the constrained edges.
//	Unmarked constrained edges will be replaced by a normal edge.
//	Marked ones will be replaced by a ConstrainingEdge. Additionally
//	associated constrained edges will be created together with the
//	new central vertex.
	for(size_t i = 0; i < cge->num_constrained_edges(); ++i){
		EdgeBase* cde = cge->constrained_edge(i);
		if(is_marked(cde)){
			refine_edge_with_hanging_vertex(cde);
		}
		else{
		//	the constrained-edge can be transformed to a normal edge
			grid.create_and_replace<Edge>(cde);
		}
	}

	cge->clear_constrained_objects();
	set_center_vertex(cge, centerVrt);
}

void HangingNodeRefinerBase::
refine_edge_with_normal_vertex(EdgeBase* e, VertexBase** newCornerVrts)
{
	HNODE_PROFILE_FUNC();

//	the grid
	Grid& grid = *m_pGrid;

	Vertex* nVrt = *grid.create<Vertex>(e);
	set_center_vertex(e, nVrt);

//	allow refCallback to calculate a new position
	if(m_refCallback)
		m_refCallback->new_vertex(nVrt, e);

//	split the edge
	vector<EdgeBase*> vEdges(2);
	e->refine(vEdges, nVrt, newCornerVrts);
	assert((vEdges.size() == 2) && "EdgeBase::refine - produced wrong number of edges.");
	grid.register_element(vEdges[0], e);
	grid.register_element(vEdges[1], e);
}

void HangingNodeRefinerBase::
refine_edge_with_hanging_vertex(EdgeBase* e, VertexBase** newCornerVrts)
{
	HNODE_PROFILE_FUNC();

	Grid& grid = *m_pGrid;
//	we have to insert a hanging node.
//	e has to be transformed to a constraining edge at the same time.
	assert(!ConstrainingEdge::type_match(e) && "invalid operation. e is a constraining edge.");
	//assert(!ConstrainedEdge::type_match(e) && "invalid operation. e is a constrained edge.");

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
	if(newCornerVrts){
		nEdge[0] = *grid.create<ConstrainedEdge>(EdgeDescriptor(newCornerVrts[0], hv), ce);
		nEdge[1] = *grid.create<ConstrainedEdge>(EdgeDescriptor(hv, newCornerVrts[1]), ce);
	}
	else{
		nEdge[0] = *grid.create<ConstrainedEdge>(EdgeDescriptor(ce->vertex(0), hv), ce);
		nEdge[1] = *grid.create<ConstrainedEdge>(EdgeDescriptor(hv, ce->vertex(1)), ce);
	}

	for(uint i = 0; i < 2; ++i)
	{
		ce->add_constrained_object(nEdge[i]);
		nEdge[i]->set_constraining_object(ce);
	}
}

void HangingNodeRefinerBase::
refine_face_with_normal_vertex(Face* f, VertexBase** newCornerVrts)
{
	HNODE_PROFILE_FUNC();

//UG_LOG("refine_face_with_normal_vertex\n");
	Grid& grid = *m_pGrid;

	VertexBase* vNewEdgeVertices[MAX_FACE_VERTICES];
	vector<Face*>		vFaces(f->num_vertices());// heuristic

	size_t numEdges = f->num_edges();
	for(size_t i = 0; i < numEdges; ++i){
		EdgeBase* e = grid.get_edge(f, i);

	//	if the face is refined with a regular rule, then every edge has to have
	//	an associated center vertex
		assert((get_mark(f) == RM_ANISOTROPIC) ||
				(get_mark(f) == RM_REGULAR && get_center_vertex(e) != NULL));

	//	assign the center vertex
		vNewEdgeVertices[i] = get_center_vertex(e);
	}

//	we'll perform a regular refine
	VertexBase* nVrt = NULL;
	/*f->refine_regular(vFaces, &nVrt, vNewEdgeVertices, NULL,
					  Vertex(), newCornerVrts);*/
	f->refine(vFaces, &nVrt, vNewEdgeVertices, NULL, newCornerVrts);

//	if a new vertex has been created during refine, then register it at the grid.
	if(nVrt)
	{
		grid.register_element(nVrt, f);

	//	allow refCallback to calculate a new position
		if(m_refCallback)
			m_refCallback->new_vertex(nVrt, f);
	}

	for(uint i = 0; i < vFaces.size(); ++i)
	{
		grid.register_element(vFaces[i], f);
	}

	set_center_vertex(f, nVrt);
}

void HangingNodeRefinerBase::
refine_face_with_hanging_vertex(Face* f, VertexBase** newCornerVrts)
{
	HNODE_PROFILE_FUNC();

//UG_LOG("refine_face_with_hanging_vertex\n");
	Grid& grid = *m_pGrid;

	size_t numVrts = f->num_vertices();
/*
	vector<EdgeBase*> 	vEdges(f->num_edges());
	vector<VertexBase*> vNewEdgeVertices(f->num_edges());
	vector<Face*>		vFaces(numVrts);// heuristic

//todo: iterate over edges directly
//	collect all associated edges.
	CollectEdges(vEdges, grid, f);
	size_t numEdges = vEdges.size();

	assert(numEdges == f->num_edges() && "ERROR in RefineFaceWithNormalVertex(...): associated edges missing.");

//	each should have an associated vertex. sort them into vNewEdgeVertices.
	for(size_t i = 0; i < numEdges; ++i)
	{
		EdgeBase* e = vEdges[i];
		int edgeIndex = GetEdgeIndex(f, e);

		assert((edgeIndex >= 0) && (edgeIndex < (int)vEdges.size()) && "ERROR in RefineFaceWithNormalVertex(...): unknown problem in CollectEdges / GetEdgeIndex.");
		//assert((get_center_vertex(e) != NULL) && "ERROR in RefineFaceWithNormalVertex(...): no new vertex on refined edge.");
		vNewEdgeVertices[edgeIndex] = get_center_vertex(e);
	}
*/
	VertexBase* vNewEdgeVertices[MAX_FACE_VERTICES];
	vector<Face*>		vFaces(f->num_vertices());// heuristic

	size_t numEdges = f->num_edges();
	for(size_t i = 0; i < numEdges; ++i){
		EdgeBase* e = grid.get_edge(f, i);

	//	if the face is refined with a regular rule, then every edge has to have
	//	an associated center vertex
		assert((get_mark(f) == RM_ANISOTROPIC) ||
				((get_mark(f) == RM_REGULAR) && (get_center_vertex(e) != NULL)));

	//	assign the center vertex
		vNewEdgeVertices[i] = get_center_vertex(e);
	}

	ConstrainingFace* cgf = NULL;
	ConstrainedVertex* hv = NULL;

//	the face has to be replaced by a constraining face.
//	we'll perform a switch here depending on the number of vertices
	switch(numVrts)
	{
		case 3:
			{
			//	create the constraining triangle and replace the old face.
				cgf = *grid.create_and_replace<ConstrainingTriangle>(f);

			//	create the constrained faces.
			//	the following triangle will not be registered at the grid. Just a temporary one.
				ConstrainedTriangle constrainedTri(cgf->vertex(0),
													cgf->vertex(1),
													cgf->vertex(2));

			//	refine the constrained tri
				VertexBase* tmpVrt;
				constrainedTri.refine(vFaces, &tmpVrt, vNewEdgeVertices,
									  NULL, newCornerVrts);
			}
			break;
		case 4:
			{
				cgf = *grid.create_and_replace<ConstrainingQuadrilateral>(f);

			//	a central hanging vertex is required
				hv = *grid.create<ConstrainedVertex>(cgf);

			//	allow refCallback to calculate a new position
				if(m_refCallback)
					m_refCallback->new_vertex(hv, cgf);

				set_center_vertex(cgf, hv);
				hv->set_constraining_object(cgf);
				cgf->add_constrained_object(hv);
				hv->set_local_coordinates(0.5, 0.5);

			//	create the constrained faces.
			//	the following quadrilateral will not be registered at the grid. Just a temporary one.
				ConstrainedQuadrilateral cdf(cgf->vertex(0), cgf->vertex(1),
											 cgf->vertex(2), cgf->vertex(3));

			//	refine the constrained quad
				VertexBase* tmpVrt;
				cdf.refine(vFaces, &tmpVrt, vNewEdgeVertices, hv, newCornerVrts);
			}
			break;
		default:
			assert(!"unsupported element type.");
			break;
	}

	if(cgf)
	{
	//	register the new faces
		for(size_t i = 0; i < vFaces.size(); ++i)
		{
			ConstrainedFace* cdf = dynamic_cast<ConstrainedFace*>(vFaces[i]);
			assert(cdf && "constrained face refine did produce faces which are not constrained.");
			if(cdf)
			{
				grid.register_element(cdf, cgf);
				cdf->set_constraining_object(cgf);
				cgf->add_constrained_object(cdf);
			}
		}

	//	we have to link the new constrained edges which have been auto-generated between the constrained faces.
	//	Since only edges that lie inside of the constraining face are newly created, and since only those
	//	have to be linked with the constraining face, the following algorithm will be ok for
	//	triangles and quadrilaterals.
	//	Check for each new edge-vertex, if an edge exists with the new center vertex or with it's next neighbor.
		for(size_t i = 0; i < numEdges; ++i)
		{
			if(hv)
			{
				ConstrainedEdge* e = dynamic_cast<ConstrainedEdge*>(grid.get_edge(vNewEdgeVertices[i], hv));
				if(e)
				{
				//	link e with the constraining face
					e->set_constraining_object(cgf);
					cgf->add_constrained_object(e);
				}
			}
			else{
			//	check if a constrained edge exists between the vertex and its next neighbor
				VertexBase* vNext = vNewEdgeVertices[(i + 1) % numEdges];
				ConstrainedEdge* e = dynamic_cast<ConstrainedEdge*>(grid.get_edge(vNewEdgeVertices[i], vNext));
				if(e)
				{
				//	link e with the constraining face
					e->set_constraining_object(cgf);
					cgf->add_constrained_object(e);
				}
			}
		}
	}
}

void HangingNodeRefinerBase::
refine_constraining_face(ConstrainingFace* cgf)
{
	HNODE_PROFILE_FUNC();

	size_t numVrts = cgf->num_vertices();

//	the grid
	Grid& grid = *m_pGrid;

//	make sure that there is one hanging vertex and two constrained edges.
	UG_ASSERT(cgf->num_constrained_edges() == numVrts,
			 "bad number of constrained edges: " << cgf->num_constrained_edges()
			 << ". There have to be as many as vertices: " << numVrts << "."
			 << "At face with center " << GetGeometricObjectCenter(grid, cgf));
	UG_ASSERT(cgf->num_constrained_faces() == 4,
			  "bad number of constrained faces. There have to be exactly 4. "
			  << "At face with center " << GetGeometricObjectCenter(grid, cgf));

	ConstrainedVertex* centralHV = NULL;
	Vertex* centerVrt = NULL;

	if(numVrts == 4){
	//	the central hanging vertex has to be transformed into a normal vertex
		centralHV = NULL;
		if(cgf->num_constrained_vertices() > 0)
			centralHV = dynamic_cast<ConstrainedVertex*>(cgf->constrained_vertex(0));
/*
		if(!centralHV){
			UG_LOG("The central hanging vertex of a constraining face is missing. ignoring face.\n");
			return;
		}
*/
	//	replace the central vertex with a normal vertex
		if(centralHV)
			centerVrt = *grid.create_and_replace<Vertex>(centralHV);
	}

//	Iterate over the constrained edges.
//	Unmarked constrained edges will be replaced by a normal edge.
//	Marked ones will be replaced by a ConstrainingEdge. Additionally
//	associated constrained edges will be created together with the
//	new central vertex.
	for(size_t i = 0; i < cgf->num_constrained_edges(); ++i){
		EdgeBase* cde = cgf->constrained_edge(i);
		if(is_marked(cde)){
			refine_edge_with_hanging_vertex(cde);
		}
		else{
		//	the constrained-edge can be transformed to a normal edge
			grid.create_and_replace<Edge>(cde);
		}
	}

//	iterate over the constrained faces.
//	If it is marked, we'll replace it by a constraining face and create
//	associated constrained faces.
//	if not, it will simply be transformed to a normal face.
//	To ease implementation we will transform it anyway and if required we will
//	call refine_face_with_hanging_vertex(...).
	for(size_t i = 0; i < cgf->num_constrained_faces(); ++i){
		Face* f = cgf->constrained_face(i);
		if(is_marked(f)){
		//	refine it using hanging_node_refinement.
			refine_face_with_hanging_vertex(f);
		}
		else{
		//	replace it by a normal face
			if(f->num_vertices() == 3)
				f = *grid.create_and_replace<Triangle>(f);
			else
				f = *grid.create_and_replace<Quadrilateral>(f);
		}
	}

//	cgf->clear_constrained_objects();
//	cgf itself now has to be transformed to a normal face
	Face* nFace;
	if(cgf->num_vertices() == 3)
		nFace = *grid.create_and_replace<Triangle>(cgf);
	else
		nFace = *grid.create_and_replace<Quadrilateral>(cgf);

	if(centerVrt)
		set_center_vertex(nFace, centerVrt);
}

void HangingNodeRefinerBase::
refine_volume_with_normal_vertex(Volume* v, VertexBase** newCornerVrts)
{
	HNODE_PROFILE_FUNC();

	Grid& grid = *m_pGrid;

	//vector<EdgeBase*> 	vEdges(v->num_edges());
	vector<VertexBase*> vNewEdgeVertices(v->num_edges());
	//vector<Face*>		vFaces(v->num_faces());
	vector<VertexBase*>	vNewFaceVertices(v->num_faces());
	vector<Volume*>		vVolumes(8);// heuristic
//	collect all associated edges.

	size_t numEdges = v->num_edges();
	for(size_t i = 0; i < numEdges; ++i){
		EdgeBase* e = grid.get_edge(v, i);
		vNewEdgeVertices[i] = get_center_vertex(e);
		UG_ASSERT(marked_anisotropic(v) || vNewEdgeVertices[i],
					"In order to fully refine a volume, all edges have"
					"to contain a new vertex!");
	}

	size_t numFaces = v->num_faces();
	for(size_t i = 0; i < numFaces; ++i){
		Face* f = grid.get_face(v, i);

		/*if(!VolumeContains(v, f))
		{
			UG_LOG("Grid::get_face(vol, ind) returned bad face.");
			MultiGrid* pmg = dynamic_cast<MultiGrid*>(m_pGrid);
			UG_LOG("Vol in level " << pmg->get_level(v));
			UG_LOG(", face in level " << pmg->get_level(f) << endl);
			UG_LOG("positions of volume vertices:");
			for(size_t i_c = 0; i_c < v->num_vertices(); ++i_c){
				UG_LOG(" " << GetGeometricObjectCenter(grid, v->vertex(i_c)));
			}
			UG_LOG("\npositions of face vertices:");
			for(size_t i_c = 0; i_c < f->num_vertices(); ++i_c){
				UG_LOG(" " << GetGeometricObjectCenter(grid, f->vertex(i_c)));
			}
			UG_LOG(endl);
		}*/

		if(f->num_vertices() == 3)
			vNewFaceVertices[i] = NULL;
		else{
			vNewFaceVertices[i] = get_center_vertex(f);
			UG_ASSERT(marked_anisotropic(v) || vNewFaceVertices[i],
					"In order to fully refine a volume, all quadrilateral sides have"
					" to contain a new vertex!"
					<< ", vol-center " << GetGeometricObjectCenter(grid, v)
					<< ", face-center " << GetGeometricObjectCenter(grid, f));
		//todo:remove this!!!
			/*if(!vNewFaceVertices[i]){
				UG_LOG("missing face-vertex: " << GetGeometricObjectCenter(grid, f) << endl);
				UG_LOG("corners of face:");
				for(size_t i_c = 0; i_c < f->num_vertices(); ++i_c){
					UG_LOG(" " << GetGeometricObjectCenter(grid, f->vertex(i_c)));
				}
				UG_LOG(endl);
				UG_LOG("during refinement of volume: " << GetGeometricObjectCenter(grid, v) << endl);
			}*/
		}
	}

//	if we're performing tetrahedral refinement, we have to collect
//	the corner coordinates, so that the refinement algorithm may choose
//	the best interior diagonal.
	vector3 corners[4];
	vector3* pCorners = NULL;
	if((v->num_vertices() == 4) && m_refCallback){
		for(size_t i = 0; i < 4; ++i){
			m_refCallback->current_pos(&corners[i].x, v->vertex(i), 3);
		}
		pCorners = corners;
	}

//	refine the volume and register new volumes at the grid.
	VertexBase* createdVrt = NULL;
	v->refine(vVolumes, &createdVrt, &vNewEdgeVertices.front(),
			  &vNewFaceVertices.front(), NULL, Vertex(), newCornerVrts, pCorners);

	if(createdVrt){
	//	register the new vertex
		grid.register_element(createdVrt, v);

	//	allow refCallback to calculate a new position
		if(m_refCallback)
			m_refCallback->new_vertex(createdVrt, v);
	}

	for(uint i = 0; i < vVolumes.size(); ++i)
		grid.register_element(vVolumes[i], v);
}

}// end of namespace
