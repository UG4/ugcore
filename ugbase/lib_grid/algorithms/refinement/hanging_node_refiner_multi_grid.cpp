// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 11.01.2011 (m,d,y)

#include <sstream>
#include "common/assert.h"
#include "lib_grid/algorithms/multi_grid_util.h"
#include "lib_grid/algorithms/selection_util.h"
#include "hanging_node_refiner_multi_grid.h"
#include "ref_mark_adjusters/mg_hnode_adjuster.h"

//	Only required for debug saves
#include "lib_grid/file_io/file_io.h"

#ifdef UG_PARALLEL
#include "pcl/pcl.h"
#include "lib_grid/parallelization/distributed_grid.h"
#include "lib_grid/parallelization/parallelization_util.h"
#endif

//define PROFILE_HANGING_NODE_REFINER_MG if you want to profile
//the refinement code.
//#define PROFILE_HANGING_NODE_REFINER_MG
#ifdef PROFILE_HANGING_NODE_REFINER_MG
	#define MGHNODE_PROFILE_FUNC()	PROFILE_FUNC()
	#define MGHNODE_PROFILE_BEGIN(name)	PROFILE_BEGIN(name)
	#define MGHNODE_PROFILE_END()	PROFILE_END()
#else
	#define MGHNODE_PROFILE_FUNC()
	#define MGHNODE_PROFILE_BEGIN(name)
	#define MGHNODE_PROFILE_END()
#endif

using namespace std;

namespace ug{

HangingNodeRefiner_MultiGrid::
HangingNodeRefiner_MultiGrid(IRefinementCallback* refCallback) :
	BaseClass(refCallback),
	m_pMG(NULL)
{
	add_ref_mark_adjuster(MGHNodeAdjuster::create());
}

HangingNodeRefiner_MultiGrid::
HangingNodeRefiner_MultiGrid(MultiGrid& mg,
						IRefinementCallback* refCallback) :
	BaseClass(refCallback),
	m_pMG(NULL)
{
	add_ref_mark_adjuster(MGHNodeAdjuster::create());
	set_grid(&mg);
}

HangingNodeRefiner_MultiGrid::
~HangingNodeRefiner_MultiGrid()
{
	set_grid(NULL);
}

void HangingNodeRefiner_MultiGrid::
grid_to_be_destroyed(Grid* grid)
{
	set_grid(NULL);
}

void HangingNodeRefiner_MultiGrid::
assign_grid(MultiGrid& mg)
{
	set_grid(&mg);
}

void HangingNodeRefiner_MultiGrid::
set_grid(MultiGrid* mg)
{
//	call base class implementation
	BaseClass::set_grid(mg);

	m_pMG = mg;
}

bool
HangingNodeRefiner_MultiGrid::
refinement_is_allowed(VertexBase* elem)
{
	return !m_pMG->has_children(elem);
}

bool
HangingNodeRefiner_MultiGrid::
refinement_is_allowed(EdgeBase* elem)
{
	return (!m_pMG->has_children(elem))
			|| elem->is_constraining();
}

bool
HangingNodeRefiner_MultiGrid::
refinement_is_allowed(Face* elem)
{
	return (!m_pMG->has_children(elem))
			|| elem->is_constraining();
}

bool
HangingNodeRefiner_MultiGrid::
refinement_is_allowed(Volume* elem)
{
	return !m_pMG->has_children(elem);
}

void HangingNodeRefiner_MultiGrid::
pre_refine()
{
	if(!m_pMG)
		throw(UGError("HangingNodeRefiner_MultiGrid::refine: No grid assigned."));

	MultiGrid& mg = *m_pMG;

//	Resize the attachment containers

	{
		selector_t& sel = get_refmark_selector();

		MGHNODE_PROFILE_BEGIN(MGHNode_ReserveAttachmentMemory);

		MGHNODE_PROFILE_BEGIN(MGHNODE_ReserveVrtData);
		mg.reserve<VertexBase>(mg.num<VertexBase>() +
					+ sel.num<VertexBase>()
					+ sel.num<EdgeBase>() - sel.num<ConstrainingEdge>()
					+ sel.num<Quadrilateral>()
					+ sel.num<ConstrainedQuadrilateral>()
					+ sel.num<Hexahedron>());
		MGHNODE_PROFILE_END();

		MGHNODE_PROFILE_BEGIN(MGHNODE_ReserveEdgeData);
		mg.reserve<EdgeBase>(mg.num<EdgeBase>()
					+ 2 * (sel.num<EdgeBase>() - sel.num<ConstrainingEdge>())
					+ 3 * (sel.num<Triangle>() + sel.num<ConstrainedTriangle>())
					+ 4 * (sel.num<Quadrilateral>() + sel.num<ConstrainedQuadrilateral>())
					+ 3 * sel.num<Prism>() + sel.num<Tetrahedron>()
					+ 4 * sel.num<Pyramid>() + 6 * sel.num<Hexahedron>());
		MGHNODE_PROFILE_END();

		MGHNODE_PROFILE_BEGIN(MGHNODE_ReserveFaceData);
		mg.reserve<Face>(mg.num<Face>()
					+ 4 * (sel.num<Face>() - sel.num<ConstrainingTriangle>()
							- sel.num<ConstrainingQuadrilateral>())
					+ 10 * sel.num<Prism>()
					+ 8 * sel.num<Tetrahedron>()
					+ 9 * sel.num<Pyramid>() + 12 * sel.num<Hexahedron>());
		MGHNODE_PROFILE_END();

		MGHNODE_PROFILE_BEGIN(MGHNODE_ReserveVolData);
		mg.reserve<Volume>(mg.num<Volume>()
					+ 8 * sel.num<Tetrahedron>() + 8 * sel.num<Prism>()
					+ 6 * sel.num<Pyramid>() + 8 * sel.num<Hexahedron>());
		MGHNODE_PROFILE_END();

		MGHNODE_PROFILE_END();
	}


//	create a child vertex for each marked vertex
	for(sel_vrt_iter iter = m_selMarkedElements.begin<VertexBase>();
		iter != m_selMarkedElements.end<VertexBase>(); ++iter)
	{
		if(marked_refine(*iter) && refinement_is_allowed(*iter)){
			VertexBase* vrt = *mg.create<Vertex>(*iter);
			if(m_refCallback)
				m_refCallback->new_vertex(vrt, *iter);
		}
	}
}

void HangingNodeRefiner_MultiGrid::
process_constraining_edge(ConstrainingEdge* cge)
{
//	call base implementation to perform refinement
	BaseClass::process_constraining_edge(cge);

//	the constrained edge is now a normal edge
	UG_ASSERT(marked_to_normal(cge), "A constraining has to be converted to a"
									 " normal edge when it is refined.");
	VertexBase* centerVrt = get_center_vertex(cge);
	Edge* e = *m_pMG->create_and_replace<Edge>(cge);
	if(centerVrt)
		set_center_vertex(e, centerVrt);
}

void HangingNodeRefiner_MultiGrid::
refine_edge_with_normal_vertex(EdgeBase* e, VertexBase** newCornerVrts)
{
//	collect child corners
	std::vector<VertexBase*> childCorners;
	collect_child_corners(childCorners, e);

//	call original implementation
	BaseClass::refine_edge_with_normal_vertex(e, &childCorners.front());
}

void HangingNodeRefiner_MultiGrid::
refine_edge_with_hanging_vertex(EdgeBase* e, VertexBase** newCornerVrts)
{
//	collect child corners
	std::vector<VertexBase*> childCorners;
	collect_child_corners(childCorners, e);

//	call original implementation
	BaseClass::refine_edge_with_hanging_vertex(e, &childCorners.front());
}

void HangingNodeRefiner_MultiGrid::
refine_face_with_normal_vertex(Face* f, VertexBase** newCornerVrts)
{
//	collect child corners
	std::vector<VertexBase*> childCorners;
	collect_child_corners(childCorners, f);

//	call original implementation
	BaseClass::refine_face_with_normal_vertex(f, &childCorners.front());
}

void HangingNodeRefiner_MultiGrid::
refine_face_with_hanging_vertex(Face* f, VertexBase** newCornerVrts)
{
//	collect child corners
	std::vector<VertexBase*> childCorners;
	collect_child_corners(childCorners, f);

//	call original implementation
	BaseClass::refine_face_with_hanging_vertex(f, &childCorners.front());
}

void HangingNodeRefiner_MultiGrid::
refine_volume_with_normal_vertex(Volume* v, VertexBase** newCornerVrts)
{
//	collect child corners
	std::vector<VertexBase*> childCorners;
	collect_child_corners(childCorners, v);

//	call original implementation
	BaseClass::refine_volume_with_normal_vertex(v, &childCorners.front());
}

VertexBase* HangingNodeRefiner_MultiGrid::
get_center_vertex(EdgeBase* e)
{
	return m_pMG->get_child_vertex(e);
}


void HangingNodeRefiner_MultiGrid::
set_center_vertex(EdgeBase* e, VertexBase* v)
{
//	the child vertex should already have been associated since the
//	multigrid is an observer itself.
	UG_ASSERT(m_pMG->get_child_vertex(e) == v, "child/vertex mismatch.");
}


VertexBase* HangingNodeRefiner_MultiGrid::
get_center_vertex(Face* f)
{
	return m_pMG->get_child_vertex(f);
}


void HangingNodeRefiner_MultiGrid::
set_center_vertex(Face* f, VertexBase* v)
{
//	the child vertex should already have been associated since the
//	multigrid is an observer itself.
	UG_ASSERT(m_pMG->get_child_vertex(f) == v, "child/vertex mismatch.");
}



void HangingNodeRefiner_MultiGrid::
restrict_selection_to_surface_coarsen_elements()
{
	restrict_selection_to_surface_coarsen_elements<VertexBase>();
	restrict_selection_to_surface_coarsen_elements<EdgeBase>();
	restrict_selection_to_surface_coarsen_elements<Face>();
	restrict_selection_to_surface_coarsen_elements<Volume>();
}

template <class TElem>
void HangingNodeRefiner_MultiGrid::
restrict_selection_to_surface_coarsen_elements()
{
	MultiGrid& mg = *m_pMG;
	selector_t& sel = get_refmark_selector();

	for(typename selector_t::template traits<TElem>::iterator iter = sel.template begin<TElem>();
		iter != sel.template end<TElem>();)
	{
		TElem* e = *iter;
		++iter;

	//	make sure that only coarsen-marks are applied
		if(get_mark(e) != RM_COARSEN){
			sel.deselect(e);
			continue;
		}

	//	make sure that the element is a surface element
		if(mg.has_children(e))
			sel.deselect(e);
	}
}


void HangingNodeRefiner_MultiGrid::
restrict_selection_to_coarsen_families()
{
	restrict_selection_to_coarsen_families<VertexBase>();
	restrict_selection_to_coarsen_families<EdgeBase>();
	restrict_selection_to_coarsen_families<Face>();
	restrict_selection_to_coarsen_families<Volume>();
}


template <class TElem>
void HangingNodeRefiner_MultiGrid::
restrict_selection_to_coarsen_families()
{
	typedef typename TElem::geometric_base_object	TBaseElem;
	MultiGrid& mg = *m_pMG;
	selector_t& sel = get_refmark_selector();

//	coarsening is only performed for complete families. If some siblings aren't
//	selected for coarsening, the whole family may not be coarsened. All siblings
//	have to be deselected in this case.

//	mark parents, which have been processed and mark siblings, which are valid
//	candidates. If thus a surface element is encountered, which has a marked
//	parent, but is not marked itself, then it is clear that it is not a valid
//	candidate.
	mg.begin_marking();
	for(typename selector_t::traits<TElem>::iterator iter = sel.begin<TElem>();
		iter != sel.end<TElem>();)
	{
		TElem* e = *iter;
		++iter;

	//	make sure that only surface elements are selected
		UG_ASSERT(mg.num_children<TBaseElem>(e) == 0,
				  "Only surface elements may be passed to this method.");

	//	make sure that only RM_COARSEN marks are used
		UG_ASSERT(get_mark(e) == RM_COARSEN,
				  "Only RM_COARSEN marks may be used in this method.");

	//	if the element is marked, we'll continue (the family is complete)
		if(mg.is_marked(e))
			continue;

	//	get the parent
		TBaseElem* parent = dynamic_cast<TBaseElem*>(mg.get_parent(e));
		if(parent){
			if(mg.is_marked(parent)){
			//	the parent is marked and e is not. We thus have to deselect e
				sel.deselect(e);
				continue;
			}

			mg.mark(parent);

		//	check whether all children of e of type TBaseElem are marked
			bool allMarked = true;
			size_t numChildren = mg.num_children<TBaseElem>(parent);
			for(size_t i = 0; i < numChildren; ++i){
				if(get_mark(mg.get_child<TBaseElem>(parent, i)) != RM_COARSEN){
					allMarked = false;
					break;
				}
			}

			if(allMarked){
			//	mark all children of parent, so that all it is clear that they
			//	belong to a complete family
				for(size_t i = 0; i < numChildren; ++i){
					mg.mark(mg.get_child<TBaseElem>(parent, i));
				}
			}
			else
				sel.deselect(e);
		}
		else{
		//	the parent of this element has a different type or it doesn't exist
		//	at all. The element is thus not regarded as a valid candidate.
			sel.deselect(e);
		}
	}
	mg.end_marking();
}

static void DeselectFamily(ISelector& sel, MultiGrid& mg, VertexBase* elem)
{
	sel.deselect(elem);
	size_t numVertexChildren = mg.num_children<VertexBase>(elem);

	for(size_t i = 0; i < numVertexChildren; ++i)
		DeselectFamily(sel, mg, mg.get_child<VertexBase>(elem, i));
}

static void DeselectFamily(ISelector& sel, MultiGrid& mg, EdgeBase* elem)
{
	sel.deselect(elem);
	size_t numVertexChildren = mg.num_children<VertexBase>(elem);
	size_t numEdgeChildren = mg.num_children<EdgeBase>(elem);

	for(size_t i = 0; i < numVertexChildren; ++i)
		DeselectFamily(sel, mg, mg.get_child<VertexBase>(elem, i));

	for(size_t i = 0; i < numEdgeChildren; ++i)
		DeselectFamily(sel, mg, mg.get_child<EdgeBase>(elem, i));
}

static void DeselectFamily(ISelector& sel, MultiGrid& mg, Face* elem)
{
	sel.deselect(elem);
	size_t numVertexChildren = mg.num_children<VertexBase>(elem);
	size_t numEdgeChildren = mg.num_children<EdgeBase>(elem);
	size_t numFaceChildren = mg.num_children<Face>(elem);

	for(size_t i = 0; i < numVertexChildren; ++i)
		DeselectFamily(sel, mg, mg.get_child<VertexBase>(elem, i));

	for(size_t i = 0; i < numEdgeChildren; ++i)
		DeselectFamily(sel, mg, mg.get_child<EdgeBase>(elem, i));

	for(size_t i = 0; i < numFaceChildren; ++i)
		DeselectFamily(sel, mg, mg.get_child<Face>(elem, i));
}

static void DeselectFamily(ISelector& sel, MultiGrid& mg, Volume* elem)
{
	sel.deselect(elem);
	size_t numVertexChildren = mg.num_children<VertexBase>(elem);
	size_t numEdgeChildren = mg.num_children<EdgeBase>(elem);
	size_t numFaceChildren = mg.num_children<Face>(elem);
	size_t numVolumeChildren = mg.num_children<Volume>(elem);

	for(size_t i = 0; i < numVertexChildren; ++i)
		DeselectFamily(sel, mg, mg.get_child<VertexBase>(elem, i));

	for(size_t i = 0; i < numEdgeChildren; ++i)
		DeselectFamily(sel, mg, mg.get_child<EdgeBase>(elem, i));

	for(size_t i = 0; i < numFaceChildren; ++i)
		DeselectFamily(sel, mg, mg.get_child<Face>(elem, i));

	for(size_t i = 0; i < numVolumeChildren; ++i)
		DeselectFamily(sel, mg, mg.get_child<Volume>(elem, i));
}

static void DeselectFamily(ISelector& sel, MultiGrid& mg, GeometricObject* elem)
{
	switch(elem->base_object_id()){
		case VERTEX:	return DeselectFamily(sel, mg, static_cast<VertexBase*>(elem));
		case EDGE:		return DeselectFamily(sel, mg, static_cast<EdgeBase*>(elem));
		case FACE:		return DeselectFamily(sel, mg, static_cast<Face*>(elem));
		case VOLUME:	return DeselectFamily(sel, mg, static_cast<Volume*>(elem));
	}
}


static void SaveCoarsenMarksToFile(ISelector& sel, const char* filename)
{
	assert(sel.grid());
	Grid& g = *sel.grid();

	SubsetHandler sh(g);

	AssignGridToSubset(g, sh, 8);
	for(VertexBaseIterator iter = g.vertices_begin(); iter != g.vertices_end(); ++iter){
		ISelector::status_t status = sel.get_selection_status(*iter);
		switch(status){
			case RM_COARSEN: sh.assign_subset(*iter, 0); break;
			case HangingNodeRefiner_MultiGrid::HNCM_ALL: sh.assign_subset(*iter, 1); break;
			case HangingNodeRefiner_MultiGrid::HNCM_PARTIAL: sh.assign_subset(*iter, 2); break;
			case HangingNodeRefiner_MultiGrid::HNCM_NONE: sh.assign_subset(*iter, 3); break;
			case HangingNodeRefiner_MultiGrid::HNCM_UNKNOWN: sh.assign_subset(*iter, 4); break;
			case HangingNodeRefiner_MultiGrid::HNCM_INVALID: sh.assign_subset(*iter, 5); break;
			case HangingNodeRefiner_MultiGrid::HNCM_REPLACE: sh.assign_subset(*iter, 6); break;
			case HangingNodeRefiner_MultiGrid::HNCM_NO_NBRS: sh.assign_subset(*iter, 7); break;
		}
	}

	for(EdgeBaseIterator iter = g.edges_begin(); iter != g.edges_end(); ++iter){
		ISelector::status_t status = sel.get_selection_status(*iter);
		switch(status){
			case RM_COARSEN: sh.assign_subset(*iter, 0); break;
			case HangingNodeRefiner_MultiGrid::HNCM_ALL: sh.assign_subset(*iter, 1); break;
			case HangingNodeRefiner_MultiGrid::HNCM_PARTIAL: sh.assign_subset(*iter, 2); break;
			case HangingNodeRefiner_MultiGrid::HNCM_NONE: sh.assign_subset(*iter, 3); break;
			case HangingNodeRefiner_MultiGrid::HNCM_UNKNOWN: sh.assign_subset(*iter, 4); break;
			case HangingNodeRefiner_MultiGrid::HNCM_INVALID: sh.assign_subset(*iter, 5); break;
			case HangingNodeRefiner_MultiGrid::HNCM_REPLACE: sh.assign_subset(*iter, 6); break;
			case HangingNodeRefiner_MultiGrid::HNCM_NO_NBRS: sh.assign_subset(*iter, 7); break;
		}
	}

	for(FaceIterator iter = g.faces_begin(); iter != g.faces_end(); ++iter){
		ISelector::status_t status = sel.get_selection_status(*iter);
		switch(status){
			case RM_COARSEN: sh.assign_subset(*iter, 0); break;
			case HangingNodeRefiner_MultiGrid::HNCM_ALL: sh.assign_subset(*iter, 1); break;
			case HangingNodeRefiner_MultiGrid::HNCM_PARTIAL: sh.assign_subset(*iter, 2); break;
			case HangingNodeRefiner_MultiGrid::HNCM_NONE: sh.assign_subset(*iter, 3); break;
			case HangingNodeRefiner_MultiGrid::HNCM_UNKNOWN: sh.assign_subset(*iter, 4); break;
			case HangingNodeRefiner_MultiGrid::HNCM_INVALID: sh.assign_subset(*iter, 5); break;
			case HangingNodeRefiner_MultiGrid::HNCM_REPLACE: sh.assign_subset(*iter, 6); break;
			case HangingNodeRefiner_MultiGrid::HNCM_NO_NBRS: sh.assign_subset(*iter, 7); break;
		}
	}

	for(VolumeIterator iter = g.volumes_begin(); iter != g.volumes_end(); ++iter){
		ISelector::status_t status = sel.get_selection_status(*iter);
		switch(status){
			case RM_COARSEN: sh.assign_subset(*iter, 0); break;
			case HangingNodeRefiner_MultiGrid::HNCM_ALL: sh.assign_subset(*iter, 1); break;
			case HangingNodeRefiner_MultiGrid::HNCM_PARTIAL: sh.assign_subset(*iter, 2); break;
			case HangingNodeRefiner_MultiGrid::HNCM_NONE: sh.assign_subset(*iter, 3); break;
			case HangingNodeRefiner_MultiGrid::HNCM_UNKNOWN: sh.assign_subset(*iter, 4); break;
			case HangingNodeRefiner_MultiGrid::HNCM_INVALID: sh.assign_subset(*iter, 5); break;
			case HangingNodeRefiner_MultiGrid::HNCM_REPLACE: sh.assign_subset(*iter, 6); break;
			case HangingNodeRefiner_MultiGrid::HNCM_NO_NBRS: sh.assign_subset(*iter, 7); break;
		}
	}

	sh.subset_info(0).name = "coarsen";
	sh.subset_info(1).name = "all";
	sh.subset_info(2).name = "partial";
	sh.subset_info(3).name = "none";
	sh.subset_info(4).name = "unknown";
	sh.subset_info(5).name = "invalid";
	sh.subset_info(6).name = "replace";
	sh.subset_info(7).name = "no_nbrs";
	sh.subset_info(8).name = "unassigned";
	AssignSubsetColors(sh);


	if(MultiGrid* mg = dynamic_cast<MultiGrid*>(&g)){
		if(mg->num<Volume>() > 0)
			SaveGridHierarchyTransformed(*mg, sh, filename, 3);
		else
			SaveGridHierarchyTransformed(*mg, sh, filename, 0.1);
	}
	else
		SaveGridToFile(g, sh, filename);
}

void HangingNodeRefiner_MultiGrid::
save_coarsen_marks_to_file(ISelector& sel, const char* filename)
{
	SaveCoarsenMarksToFile(sel, filename);
}


///	temporary method, which will be removed when debugging is done.
void HangingNodeRefiner_MultiGrid::
debug_save(ISelector& sel, const char* filePrefix)
{
	if(debugging_enabled()){
		stringstream ss;
		ss << filePrefix << "_p";
		#ifdef UG_PARALLEL
			ss << pcl::GetProcRank();
		#else
			ss << "0";
		#endif
		ss << ".ugx";
		//UG_LOG("Saving coarsen debug marks: " << ss.str() << endl);
		save_coarsen_marks_to_file(sel, ss.str().c_str());
	}
}

///	temporary method, which will be removed when debugging is done.
static void ContinuousDebugSave(ISelector& sel)
{
	static int counter = 0;
	stringstream ss;
	ss << "coarsen_debug_" << counter << "_p";
	#ifdef UG_PARALLEL
		ss << pcl::GetProcRank();
	#else
		ss << "0";
	#endif
	ss << ".ugx";
	UG_LOG("Saving coarsen debug marks: " << ss.str() << endl);
	SaveCoarsenMarksToFile(sel, ss.str().c_str());
	++counter;
}

///	temporary method, which will be removed when debugging is done.
static void ParallelLayoutDebugSave(MultiGrid& mg)
{
	static int counter = 0;
	stringstream ss;
	ss << "coarsen_parallel_layout_" << counter << "_p";
	#ifdef UG_PARALLEL
		ss << pcl::GetProcRank();
	#else
		ss << "0";
	#endif
	ss << ".ugx";
	UG_LOG("Saving coarsen parallel layout: " << ss.str() << endl);
	if(mg.num<Volume>() > 0)
		SaveParallelGridLayout(mg, ss.str().c_str(), 3);
	else
		SaveParallelGridLayout(mg, ss.str().c_str(), 0.1);
	++counter;
}


//void HangingNodeRefiner_MultiGrid::
//collect_objects_for_coarsen(bool scheduleCoarseningBeginsMessage)
//{
//	MultiGrid& mg = *m_pMG;
//	selector_t& sel = get_refmark_selector();
//
//	bool gotVols = contains_volumes();
//	bool gotFaces = contains_faces();
//	bool gotEdges = contains_edges();
//
//debug_save(sel, "coarsen_marks_01_start");
////	deselect all elements of lower dimension, since we currently can't handle them correctly...
////todo:	A more sophisticated approach could be used. E.g. only deselect elements if they
////		are connected to an element of higher dimension.
//	if(gotVols)
//		sel.clear<Face>();
//
//	if(gotFaces)
//		sel.clear<EdgeBase>();
//
//	if(gotEdges)
//		sel.clear<VertexBase>();
//
//debug_save(sel, "coarsen_marks_01.5_only_full_dim_selected");
//	restrict_selection_to_surface_coarsen_elements();
//	copy_marks_to_vmasters(true, true, true, true);
//	restrict_selection_to_coarsen_families();
//	copy_marks_to_vslaves(true, true, true, true);
//
//	SelectAssociatedEdges(sel, sel.begin<Face>(), sel.end<Face>(), HNCM_UNKNOWN);
//	SelectAssociatedEdges(sel, sel.begin<Volume>(), sel.end<Volume>(), HNCM_UNKNOWN);
//	broadcast_marks_horizontally(false, true, false);
//	copy_marks_to_vmasters(false, true, false, false);
//
//debug_save(sel, "coarsen_marks_02_restricted_to_surface_families");
//
//	typedef vector<EdgeBase*> EdgeVec;
//	EdgeVec vedges;
//	vedges.assign(sel.begin<EdgeBase>(), sel.end<EdgeBase>());
//
//	Grid::edge_traits::secure_container edges;
//	Grid::face_traits::secure_container	faces;
//	Grid::volume_traits::secure_container vols;
//
////	make sure that coarsening won't generate constraining edges with more than
////	one constrained vertex.
//	bool running = gotVols || gotFaces;
//	while(running){
//	//	classify unknown edges
//		for(EdgeVec::iterator iter = vedges.begin(); iter != vedges.end(); ++iter)
//		{
//			EdgeBase* e = *iter;
//
//			if(sel.get_selection_status(e) == HNCM_UNKNOWN){
//				size_t numSel = 0;
//				size_t numNbrs = 0;
//				if(gotVols){
//					mg.associated_elements(vols, e);
//					numNbrs = vols.size();
//					for(size_t i = 0; i < numNbrs; ++i){
//						if(sel.is_selected(vols[i]))
//							++numSel;
//					}
//				}
//				else if(gotFaces){
//					mg.associated_elements(faces, e);
//					numNbrs = faces.size();
//					for(size_t i = 0; i < numNbrs; ++i){
//						if(sel.is_selected(faces[i]))
//							++numSel;
//					}
//				}
//
//				if(numNbrs > 0){
//					if(numSel == 0)
//						sel.select(e, HNCM_NONE);
//					else if(numSel < numNbrs)
//						sel.select(e, HNCM_PARTIAL);
//					else
//						sel.select(e, HNCM_ALL);
//				}
//				else
//					sel.select(e, HNCM_NO_NBRS);
//			}
//		}
//
//		broadcast_marks_horizontally(false, true, false);
//		copy_marks_to_vmasters(false, true, false, false);
//
//	//	clear all edges which are marked with HNCM_NONE from the selection
//		for(selector_t::traits<EdgeBase>::iterator iter = sel.begin<EdgeBase>();
//			iter != sel.end<EdgeBase>();)
//		{
//			EdgeBase* e = *iter;
//			++iter;
//			if(sel.get_selection_status(e) == HNCM_NONE)
//				sel.deselect(e);
//		}
//
//	//	if a constrained edge of a constraining edge won't be fully coarsened,
//	//	then the constraining edge may not be coarsened at all.
//		bool foundInvalid = false;
//		for(selector_t::traits<ConstrainingEdge>::iterator iter = sel.begin<ConstrainingEdge>();
//			iter != sel.end<ConstrainingEdge>(); ++iter)
//		{
//			ConstrainingEdge* cge = *iter;
//			for(size_t i = 0; i < cge->num_constrained_edges(); ++i){
//				if(sel.get_selection_status(cge->constrained_edge(i)) != HNCM_ALL){
//					foundInvalid = true;
//					sel.select(cge, HNCM_INVALID);
//					break;
//				}
//			}
//		}
//
//	//	has anybody found an invalid parent? If not, exit adjustment.
//		if(!one_proc_true(foundInvalid)){
//			running = false;
//			break;
//		}
//
//	//	clear the current candidate vector
//		vedges.clear();
//
//
//	//	exchange invalid marks
//	//	we have to copy them to vmasters, since constraining ghost edges don't
//	//	have an associated constrained edge from which they could obtain their
//	//	invalid mark.
//		copy_marks_to_vmasters(false, true, false, false);
//
//		for(selector_t::traits<EdgeBase>::iterator iter = sel.begin<EdgeBase>();
//			iter != sel.end<EdgeBase>();)
//		{
//			EdgeBase* e = *iter;
//			++iter;
//			if(sel.get_selection_status(e) == HNCM_INVALID){
//				if(gotVols){
//					mg.associated_elements(vols, e);
//					for(size_t i_vol = 0; i_vol < vols.size(); ++i_vol){
//						if(!sel.is_selected(vols[i_vol]))
//							continue;
//						GeometricObject* parent = mg.get_parent(vols[i_vol]);
//						if(parent){
//							DeselectFamily(sel, mg, parent);
//							mg.associated_elements(edges, parent);
//							for(size_t i_edge = 0; i_edge < edges.size(); ++i_edge){
//								EdgeBase* e = edges[i_edge];
//								if(sel.is_selected(e))
//									sel.select(e, HNCM_UNKNOWN);
//								for(size_t i = 0; i < mg.num_children<EdgeBase>(e); ++i){
//									EdgeBase* child = mg.get_child<EdgeBase>(e, i);
//									if(sel.is_selected(child)){
//										sel.select(child, HNCM_UNKNOWN);
//									}
//								}
//							}
//
//						//	we also have to mark edge-children of side-faces as unknown
//							mg.associated_elements(faces, parent);
//							for(size_t i_face = 0; i_face < faces.size(); ++i_face){
//								Face* f = faces[i_face];
//								if(sel.is_selected(f))
//									sel.select(f, HNCM_UNKNOWN);
//								for(size_t i = 0; i < mg.num_children<EdgeBase>(f); ++i){
//									EdgeBase* child = mg.get_child<EdgeBase>(f, i);
//									if(sel.is_selected(child)){
//										sel.select(child, HNCM_UNKNOWN);
//									}
//								}
//							}
//
//						}
//					}
//				}
//				else if(gotFaces){
//					mg.associated_elements(faces, e);
//					for(size_t i_face = 0; i_face < faces.size(); ++i_face){
//						if(!sel.is_selected(faces[i_face]))
//							continue;
//						GeometricObject* parent = mg.get_parent(faces[i_face]);
//						if(parent){
//							DeselectFamily(sel, mg, parent);
//							mg.associated_elements(edges, parent);
//							for(size_t i_edge = 0; i_edge < edges.size(); ++i_edge){
//								if(sel.is_selected(edges[i_edge]))
//									sel.select(edges[i_edge], HNCM_UNKNOWN);
//								for(size_t i = 0; i < mg.num_children<EdgeBase>(edges[i_edge]); ++i){
//									EdgeBase* child = mg.get_child<EdgeBase>(edges[i_edge], i);
//									if(sel.is_selected(child)){
//										sel.select(child, HNCM_UNKNOWN);
//									}
//								}
//							}
//						}
//					}
//				}
//
//			//	mark the formerly invalid edge as unknown
//				sel.select(e, HNCM_UNKNOWN);
//			}
//		}
//
//		if(debugging_enabled()){
//			ContinuousDebugSave(sel);
//		}
//
//	//	since h-interfaces do not always exist between v-master copies,
//	//	we have to perform an additional broadcast, to make sure, that all copies
//	//	have the correct marks.
//		broadcast_marks_vertically(false, true, true, true, true);
//		if(debugging_enabled()){
//			ContinuousDebugSave(sel);
//		}
//
//
//
//	//	find edges which were marked as unknown and prepare qedges for the next run
//		for(selector_t::traits<EdgeBase>::iterator iter = sel.begin<EdgeBase>();
//			iter != sel.end<EdgeBase>(); ++iter)
//		{
//			if(sel.get_selection_status(*iter) == HNCM_UNKNOWN)
//				vedges.push_back(*iter);
//		}
//	}
//
//debug_save(sel, "coarsen_marks_03_irregularities_resolved");
//
////	finally classify faces and vertices
//	SelectAssociatedFaces(sel, sel.begin<Volume>(), sel.end<Volume>(), HNCM_ALL);
//	SelectAssociatedVertices(sel, sel.begin<EdgeBase>(), sel.end<EdgeBase>(), HNCM_ALL);
//
//	broadcast_marks_horizontally(true, false, true);
//	copy_marks_to_vmasters(true, false, true, false);
//
//	if(gotVols){
//		for(selector_t::traits<Face>::iterator iter = sel.begin<Face>();
//			iter != sel.end<Face>();)
//		{
//			Face* f = *iter;
//			++iter;
//
//			mg.associated_elements(vols, f);
//			for(size_t i = 0; i < vols.size(); ++i){
//				if(!sel.is_selected(vols[i])){
//					sel.select(f, HNCM_PARTIAL);
//					break;
//				}
//			}
//		}
//	}
//
//	for(selector_t::traits<VertexBase>::iterator iter = sel.begin<VertexBase>();
//		iter != sel.end<VertexBase>();)
//	{
//		VertexBase* v = *iter;
//		++iter;
//
//		mg.associated_elements(edges, v);
//		for(size_t i = 0; i < edges.size(); ++i){
//			if(sel.get_selection_status(edges[i]) != HNCM_ALL){
//				sel.select(v, HNCM_PARTIAL);
//				break;
//			}
//		}
//	}
//
//	broadcast_marks_horizontally(true, false, true);
//	copy_marks_to_vmasters(true, false, true, false);
//
//debug_save(sel, "coarsen_marks_04_faces_and_vertices_classified");
//
//
////	ATTENTION:	We post the message that coarsening begins at this point, since
////				we don't want REPLACE elements to be selected and thus contained
////				in the specified geometric-object-collection.
//	if(scheduleCoarseningBeginsMessage){
//		m_messageHub->post_message(GridMessage_Adaption(GMAT_HNODE_COARSENING_BEGINS,
//										m_selMarkedElements.get_geometric_objects()));
//	}
//
//
////	select unselected constraining edges and faces of selected constrained
////	edges and faces with HNCM_REPLACE to indicate, that they have to be
////	transformed to a normal edge
////	mark parents of normal and constraining edges and faces which were marked for
////	partial refinement with a replace mark
//	for(selector_t::traits<EdgeBase>::iterator
//		iter = sel.begin<EdgeBase>(); iter != sel.end<EdgeBase>(); ++iter)
//	{
//		EdgeBase* e = *iter;
//		if(GeometricObject* parent = mg.get_parent(e)){
//			bool isConstrained = e->is_constrained();
//			if((isConstrained && (sel.get_selection_status(e) == HNCM_ALL)) ||
//			    ((!isConstrained) && (sel.get_selection_status(*iter) == HNCM_PARTIAL)))
//			{
//				if(!sel.is_selected(parent))
//					sel.select(parent, HNCM_REPLACE);
//			}
//		}
//	}
//
//	for(selector_t::traits<Face>::iterator
//		iter = sel.begin<Face>(); iter != sel.end<Face>(); ++iter)
//	{
//		Face* e = *iter;
//		if(GeometricObject* parent = mg.get_parent(e)){
//			bool isConstrained = e->is_constrained();
//			if((isConstrained && (sel.get_selection_status(e) == HNCM_ALL)) ||
//			   ((!isConstrained) && (sel.get_selection_status(*iter) == HNCM_PARTIAL)))
//			{
//				if(!sel.is_selected(parent))
//					sel.select(parent, HNCM_REPLACE);
//			}
//		}
//	}
//
//	broadcast_marks_horizontally(false, true, true);
//	copy_marks_to_vmasters(false, true, true, false);
//
//debug_save(sel, "coarsen_marks_05_adjustment_done");
//}

void HangingNodeRefiner_MultiGrid::
collect_objects_for_coarsen(bool scheduleCoarseningBeginsMessage)
{
	PROFILE_FUNC();

	MultiGrid& mg = *m_pMG;
	selector_t& sel = get_refmark_selector();

	bool gotVols = contains_volumes();
	bool gotFaces = contains_faces();
	bool gotEdges = contains_edges();

debug_save(sel, "coarsen_marks_01_start");
//	deselect all elements of lower dimension, since we currently can't handle them correctly...
//todo:	A more sophisticated approach could be used. E.g. only deselect elements if they
//		are connected to an element of higher dimension.
	if(gotVols)
		sel.clear<Face>();

	if(gotFaces)
		sel.clear<EdgeBase>();

	if(gotEdges)
		sel.clear<VertexBase>();

	restrict_selection_to_surface_coarsen_elements();
	copy_marks_to_vmasters(true, true, true, true);
	restrict_selection_to_coarsen_families();
	copy_marks_to_vslaves(true, true, true, true);

	SelectAssociatedVertices(sel, sel.begin<EdgeBase>(), sel.end<EdgeBase>(), HNCM_UNKNOWN);
	SelectAssociatedVertices(sel, sel.begin<Face>(), sel.end<Face>(), HNCM_UNKNOWN);
	SelectAssociatedVertices(sel, sel.begin<Volume>(), sel.end<Volume>(), HNCM_UNKNOWN);
	broadcast_marks_horizontally(true, false, false);
	copy_marks_to_vmasters(true, false, false, false);

debug_save(sel, "coarsen_marks_02_restricted_to_surface_families");

	typedef vector<VertexBase*> VrtVec;
	VrtVec vvrts;
	vvrts.assign(sel.begin<VertexBase>(), sel.end<VertexBase>());

	Grid::vertex_traits::secure_container vrts;
	Grid::edge_traits::secure_container edges;
	Grid::face_traits::secure_container	faces;
	Grid::volume_traits::secure_container vols;

//	make sure that coarsening won't generate situations in which surface elements
//	which meet at a common vertex have a level difference of more than 1.
	bool running = gotVols || gotFaces;
	while(running){
	//	classify unknown vertices
		for(VrtVec::iterator iter = vvrts.begin(); iter != vvrts.end(); ++iter)
		{
			VertexBase* e = *iter;

			if(sel.get_selection_status(e) == HNCM_UNKNOWN){
				size_t numSel = 0;
				size_t numNbrs = 0;
				if(gotVols){
					mg.associated_elements(vols, e);
					numNbrs = vols.size();
					for(size_t i = 0; i < numNbrs; ++i){
						if(sel.is_selected(vols[i]))
							++numSel;
					}
				}
				else if(gotFaces){
					mg.associated_elements(faces, e);
					numNbrs = faces.size();
					for(size_t i = 0; i < numNbrs; ++i){
						if(sel.is_selected(faces[i]))
							++numSel;
					}
				}

				if(numNbrs > 0){
					if(numSel == 0)
						sel.select(e, HNCM_NONE);
					else if(numSel < numNbrs)
						sel.select(e, HNCM_PARTIAL);
					else
						sel.select(e, HNCM_ALL);
				}
				else
					sel.select(e, HNCM_NO_NBRS);
			}
		}

		broadcast_marks_horizontally(true, false, false);
		copy_marks_to_vmasters(true, false, false, false);

	//	clear all vertices which are marked with HNCM_NONE from the selection
		for(selector_t::traits<VertexBase>::iterator iter = sel.begin<VertexBase>();
			iter != sel.end<VertexBase>();)
		{
			VertexBase* e = *iter;
			++iter;
			if(sel.get_selection_status(e) == HNCM_NONE)
				sel.deselect(e);
		}

	//	if a vertex has a child which will not be fully coarsed then the
	//	vertex itself may not be coarsened at all (nor neighboring elements)
		bool foundInvalid = false;
		for(selector_t::traits<VertexBase>::iterator iter = sel.begin<VertexBase>();
			iter != sel.end<VertexBase>(); ++iter)
		{
			VertexBase* e = *iter;
			VertexBase* c = mg.get_child_vertex(e);
			if(c && sel.get_selection_status(c) != HNCM_ALL){
				foundInvalid = true;
				sel.select(e, HNCM_INVALID);
			}
		}

	//	has anybody marked an element as invalid? If not, exit adjustment.
		if(!one_proc_true(foundInvalid)){
			running = false;
			break;
		}

	//	clear the current candidate vector
		vvrts.clear();

	//	exchange invalid marks
	//	we have to copy them to vmasters, since in some situations there may be
	//	ghost-vertices which are shadows of vertices which lie on the rim of a level...
		copy_marks_to_vmasters(true, false, false, false);

	//	we have to deselect all families which are connected to invalid vertices
		for(selector_t::traits<VertexBase>::iterator iter = sel.begin<VertexBase>();
			iter != sel.end<VertexBase>();)
		{
			VertexBase* e = *iter;
			++iter;
			if(sel.get_selection_status(e) == HNCM_INVALID){
				if(gotVols){
					mg.associated_elements(vols, e);
					for(size_t i_vol = 0; i_vol < vols.size(); ++i_vol){
						if(!sel.is_selected(vols[i_vol]))
							continue;
						GeometricObject* parent = mg.get_parent(vols[i_vol]);
						if(parent){
							DeselectFamily(sel, mg, parent);
							mg.associated_elements(vrts, parent);
							for(size_t i_vrt = 0; i_vrt < vrts.size(); ++i_vrt){
								VertexBase* vrt = vrts[i_vrt];
								if(sel.is_selected(vrt))
									sel.select(vrt, HNCM_UNKNOWN);

								VertexBase* c = mg.get_child_vertex(vrt);
								if(c && sel.is_selected(c)){
									sel.select(c, HNCM_UNKNOWN);
								}
							}

						//	we also have to mark vertex-children of side-faces
						//	and side-edges as unknown
							mg.associated_elements(edges, parent);
							for(size_t i_edge = 0; i_edge < edges.size(); ++i_edge){
								EdgeBase* edge = edges[i_edge];
								if(sel.is_selected(edge))
									sel.select(edge, HNCM_UNKNOWN);

								VertexBase* c = mg.get_child_vertex(edge);
								if(c && sel.is_selected(c)){
									sel.select(c, HNCM_UNKNOWN);
								}
							}

							mg.associated_elements(faces, parent);
							for(size_t i_face = 0; i_face < faces.size(); ++i_face){
								Face* f = faces[i_face];
								if(sel.is_selected(f))
									sel.select(f, HNCM_UNKNOWN);

								VertexBase* c = mg.get_child_vertex(f);
								if(c && sel.is_selected(c)){
									sel.select(c, HNCM_UNKNOWN);
								}
							}

						}
					}
				}
				else if(gotFaces){
					mg.associated_elements(faces, e);
					for(size_t i_face = 0; i_face < faces.size(); ++i_face){
						if(!sel.is_selected(faces[i_face]))
							continue;
						GeometricObject* parent = mg.get_parent(faces[i_face]);
						if(parent){
							DeselectFamily(sel, mg, parent);
							mg.associated_elements(vrts, parent);
							for(size_t i_vrt = 0; i_vrt < vrts.size(); ++i_vrt){
								VertexBase* vrt = vrts[i_vrt];
								if(sel.is_selected(vrt))
									sel.select(vrt, HNCM_UNKNOWN);

								VertexBase* c = mg.get_child_vertex(vrt);
								if(c && sel.is_selected(c)){
									sel.select(c, HNCM_UNKNOWN);
								}
							}

						//	we also have to mark vertex-children of side-edges
							mg.associated_elements(edges, parent);
							for(size_t i_edge = 0; i_edge < edges.size(); ++i_edge){
								EdgeBase* edge = edges[i_edge];
								if(sel.is_selected(edge))
									sel.select(edge, HNCM_UNKNOWN);

								VertexBase* c = mg.get_child_vertex(edge);
								if(c && sel.is_selected(c)){
									sel.select(c, HNCM_UNKNOWN);
								}
							}
						}
					}
				}

			//	mark the formerly invalid vertex as unknown
				sel.select(e, HNCM_UNKNOWN);
			}
		}

		if(debugging_enabled()){
			ContinuousDebugSave(sel);
		}

	//	since h-interfaces do not always exist between v-master copies,
	//	we have to perform an additional broadcast, to make sure, that all copies
	//	have the correct marks.
		broadcast_marks_vertically(true, true, true, true, true);
		if(debugging_enabled()){
			ContinuousDebugSave(sel);
		}



	//	find vertices which were marked as unknown and prepare the candidates array for the next run
		for(selector_t::traits<VertexBase>::iterator iter = sel.begin<VertexBase>();
			iter != sel.end<VertexBase>(); ++iter)
		{
			if(sel.get_selection_status(*iter) == HNCM_UNKNOWN)
				vvrts.push_back(*iter);
		}
	}

debug_save(sel, "coarsen_marks_03_irregularities_resolved");

//	finally classify faces and edges
	SelectAssociatedFaces(sel, sel.begin<Volume>(), sel.end<Volume>(), HNCM_ALL);
	SelectAssociatedEdges(sel, sel.begin<Face>(), sel.end<Face>(), HNCM_ALL);

//	we only have to communicate face marks if volumes exist
	if(gotVols){
		broadcast_marks_horizontally(false, true, true);
		copy_marks_to_vmasters(false, true, true, false);
	}
	else{
		broadcast_marks_horizontally(false, true, false);
		copy_marks_to_vmasters(false, true, false, false);
	}

	if(gotVols){
		for(selector_t::traits<Face>::iterator iter = sel.begin<Face>();
			iter != sel.end<Face>();)
		{
			Face* f = *iter;
			++iter;

			mg.associated_elements(vols, f);
			for(size_t i = 0; i < vols.size(); ++i){
				if(!sel.is_selected(vols[i])){
					sel.select(f, HNCM_PARTIAL);
					break;
				}
			}
		}

		for(selector_t::traits<EdgeBase>::iterator iter = sel.begin<EdgeBase>();
			iter != sel.end<EdgeBase>();)
		{
			EdgeBase* e = *iter;
			++iter;

			mg.associated_elements(vols, e);
			for(size_t i = 0; i < vols.size(); ++i){
				if(!sel.is_selected(vols[i])){
					sel.select(e, HNCM_PARTIAL);
					break;
				}
			}
		}
	}
	else if(gotFaces){
		for(selector_t::traits<EdgeBase>::iterator iter = sel.begin<EdgeBase>();
			iter != sel.end<EdgeBase>();)
		{
			EdgeBase* e = *iter;
			++iter;

			mg.associated_elements(faces, e);
			for(size_t i = 0; i < faces.size(); ++i){
				if(!sel.is_selected(faces[i])){
					sel.select(e, HNCM_PARTIAL);
					break;
				}
			}
		}
	}

	if(gotVols){
		broadcast_marks_horizontally(false, true, true);
		copy_marks_to_vmasters(false, true, true, false);
	}
	else{
		broadcast_marks_horizontally(false, true, false);
		copy_marks_to_vmasters(false, true, false, false);
	}

debug_save(sel, "coarsen_marks_04_faces_and_vertices_classified");


//	ATTENTION:	We post the message that coarsening begins at this point, since
//				we don't want REPLACE elements to be selected and thus contained
//				in the specified geometric-object-collection.
	if(scheduleCoarseningBeginsMessage){
		m_messageHub->post_message(GridMessage_Adaption(GMAT_HNODE_COARSENING_BEGINS,
										m_selMarkedElements.get_geometric_objects()));
	}


//	select unselected constraining edges and faces of selected constrained
//	edges and faces with HNCM_REPLACE to indicate, that they have to be
//	transformed to a normal edge
//	mark parents of normal and constraining edges and faces which were marked for
//	partial refinement with a replace mark
	for(selector_t::traits<EdgeBase>::iterator
		iter = sel.begin<EdgeBase>(); iter != sel.end<EdgeBase>(); ++iter)
	{
		EdgeBase* e = *iter;
		if(GeometricObject* parent = mg.get_parent(e)){
			bool isConstrained = e->is_constrained();
			if((isConstrained && (sel.get_selection_status(e) == HNCM_ALL)) ||
			    ((!isConstrained) && (sel.get_selection_status(*iter) == HNCM_PARTIAL)))
			{
				if(!sel.is_selected(parent))
					sel.select(parent, HNCM_REPLACE);
			}
		}
	}

	if(gotVols){
		for(selector_t::traits<Face>::iterator
			iter = sel.begin<Face>(); iter != sel.end<Face>(); ++iter)
		{
			Face* e = *iter;
			if(GeometricObject* parent = mg.get_parent(e)){
				bool isConstrained = e->is_constrained();
				if((isConstrained && (sel.get_selection_status(e) == HNCM_ALL)) ||
				   ((!isConstrained) && (sel.get_selection_status(*iter) == HNCM_PARTIAL)))
				{
					if(!sel.is_selected(parent))
						sel.select(parent, HNCM_REPLACE);
				}
			}
		}
	}

	if(gotVols){
		broadcast_marks_horizontally(false, true, true);
		copy_marks_to_vmasters(false, true, true, false);
	}
	else{
		broadcast_marks_horizontally(false, true, false);
		copy_marks_to_vmasters(false, true, false, false);
	}

debug_save(sel, "coarsen_marks_05_adjustment_done");
}

bool HangingNodeRefiner_MultiGrid::
perform_coarsening()
{
	PROFILE_FUNC();
/*
We have to handle elements as follows:

  normal edge:
	All nbrs marked: remove
	Some nbrs marked (hnode):
	  If has (unconstrained) parent:
	    convert to constrained. If parent is normal, then make parent constraining.
	    add edge to constrained elements of parent.
	  else
	    keep edge

  constrained edge:
    All nbrs marked: remove
    Some nbrs marked (hnode): keep (important for 3d)

  constraining edge (hnode) (only marked, if associated constrained edges are marked, too):
    (All nbrs can't be marked)
	Some nbrs marked:
	  If has (unconstrained) parent:
		convert to constrained. If parent is normal, then make parent constraining.
		add edge to constrained elements of parent.
	  else
		convert to normal edge

  constraining edge (no hnode - no nbrs marked, but children are marked):
    Convert to normal edge
*/
	if(!m_pMG)
		return false;

	MultiGrid& mg = *m_pMG;
	selector_t& sel = get_refmark_selector();

//	We let collect-objects-for-coarsen schedule the message that coarsening begins
	collect_objects_for_coarsen(true);

//	if a debug file was specified, we'll now save the marks to that file
	if(!m_adjustedMarksDebugFilename.empty())
		save_coarsen_marks_to_file(sel, m_adjustedMarksDebugFilename.c_str());

//	inform derived classes that coarsening begins
	pre_coarsen();

////////////////////
//	VOLUMES
//	erase elements and introduce constraining/constrained elements as required.
	mg.erase(sel.begin<Volume>(), sel.end<Volume>());

////////////////////
//	FACES
	for(selector_t::traits<Face>::iterator iter = sel.begin<Face>();
		iter != sel.end<Face>();)
	{
		Face* elem = *iter;
		++iter;
		ISelector::status_t selState = sel.get_selection_status(elem);
		sel.deselect(elem);

		switch(selState){
			case RM_COARSEN:
			case HNCM_ALL:{
			//	this should only be set on normal or constrained faces
				assert(!elem->is_constraining());
				mg.erase(elem);
			}break;

			case HNCM_PARTIAL:{
				if(elem->is_constrained()){
				//	a constrained element where not all nbrs are selected
				//	will be kept as it is.
					continue;
				}

				UG_ASSERT(mg.parent_type(elem) == FACE,
						  "Only a face with a parent-face may be marked for "
						  "partial coarsening!");

			//	replace elem by a constrained face
				ConstrainedFace* cdf = NULL;
				switch(elem->reference_object_id()){
					case ROID_TRIANGLE:
						cdf = *mg.create_and_replace<ConstrainedTriangle>(elem);
						break;
					case ROID_QUADRILATERAL:
						cdf = *mg.create_and_replace<ConstrainedQuadrilateral>(elem);
						break;
					default:
						UG_THROW("Unknown face reference object type in "
								 "HangingNodeRefiner_MultiGrid::coarsen");
						break;
				}

				elem = cdf;


			//	the face has to be replaced by a constrained face.
			//	Volume-children don't have to be considered here, since
			//	they should all be marked with HNCM_ALL_NBRS_SELECTED
				Face* parent = dynamic_cast<Face*>(mg.get_parent(elem));

				if(parent){
				//	the parent may not be a constrained object
					UG_ASSERT(!parent->is_constrained(),
							  "An element with a constrained parent may not be marked"
							  " for partial coarsening!");

				//	check whether the parent is already constraining.
				//	if not, it will be transformed to a constraining element later on
				//	and the connection is established then.
					ConstrainingFace* cf = dynamic_cast<ConstrainingFace*>(parent);
					if(cf){
						cf->add_constrained_object(cdf);
						cdf->set_constraining_object(cf);
					}
				}
				else{
				//	There is no parent on the local process. We thus only set the
				//	type of the constraining element
					cdf->set_parent_base_object_id(mg.parent_type(cdf));
				}

			}break;

			case HNCM_REPLACE:{
				assert(!elem->is_constrained());
				if(elem->is_constraining()){
					#ifdef UG_DEBUG
					//	make sure that the associated constrained vertex will be removed
						ConstrainingFace* cf = dynamic_cast<ConstrainingFace*>(elem);
						if(cf->num_constrained_vertices() > 0){
							VertexBase* hv = cf->constrained_vertex(0);
							UG_ASSERT(sel.get_selection_status(hv) == HNCM_ALL,
									  "A constraining face may only be replaced by a normal face"
									  " during coarsening, if its constrained vertex will be removed."
									  " At: " << ElementDebugInfo(mg, elem));
						}
					#endif

				//	the constraining face has to be transformed to a normal face
					switch(elem->reference_object_id()){
						case ROID_TRIANGLE:
							mg.create_and_replace<Triangle>(elem);
							break;
						case ROID_QUADRILATERAL:
							mg.create_and_replace<Quadrilateral>(elem);
							break;
						default:
							UG_THROW("Unknown face reference object type in "
									 "HangingNodeRefiner_MultiGrid::coarsen");
							break;
					}
				}
				else{
				//	transform the face to a constraining face
					ConstrainingFace* cf = NULL;
					switch(elem->reference_object_id()){
						case ROID_TRIANGLE:
							cf = *mg.create_and_replace<ConstrainingTriangle>(elem);
							break;
						case ROID_QUADRILATERAL:
							cf = *mg.create_and_replace<ConstrainingQuadrilateral>(elem);
							break;
						default:
							UG_THROW("Unknown face reference object type in "
									 "HangingNodeRefiner_MultiGrid::coarsen");
							break;
					}

				//	make sure that a connection with constrained children is established
				//	note - associated constrained vertices and constrained edges can
				//	not exist at this point.
					for(size_t i = 0; i < mg.num_children<Face>(cf); ++i){
						if(ConstrainedFace* cdf =
						   dynamic_cast<ConstrainedFace*>(mg.get_child<Face>(cf, i)))
						{
							cf->add_constrained_object(cdf);
							cdf->set_constraining_object(cf);
						}
					}
				}

			}break;

			default:
				UG_THROW("Bad selection mark on face during "
						 "HangingNodeRefiner_MultiGrid::coarsen. Internal error.");
				break;
		}
	}

////////////////////
//	EDGES
	for(selector_t::traits<EdgeBase>::iterator iter = sel.begin<EdgeBase>();
		iter != sel.end<EdgeBase>();)
	{
		EdgeBase* elem = *iter;
		++iter;
		ISelector::status_t selState = sel.get_selection_status(elem);
		sel.deselect(elem);

		switch(selState){
			case RM_COARSEN:
			case HNCM_ALL:{
			//	this should only be set on normal or constrained edges
				assert(!elem->is_constraining());
				mg.erase(elem);
			}break;

			case HNCM_PARTIAL:{
				if(elem->is_constrained()){
				//	a constrained element where not all nbrs are selected
				//	will be kept as it is.
					continue;
				}

				UG_ASSERT(mg.parent_type(elem) != VOLUME,
						  "Children of volume elements may not be marked for"
						  " partial coarsening!");

			//	replace the edge by a constrained edge
				ConstrainedEdge* cde = *mg.create_and_replace<ConstrainedEdge>(elem);
				elem = cde;

			//	establish connection to constraining element
				GeometricObject* parent = mg.get_parent(cde);
				if(parent){
					switch(parent->base_object_id()){
						case EDGE:{
								EdgeBase* edgeParent = dynamic_cast<EdgeBase*>(parent);
								UG_ASSERT(edgeParent, "Severe type error during coarsing - bad parent type: "
										  << ElementDebugInfo(mg, elem));

							//	check whether the parent is already constraining
								ConstrainingEdge* ce = dynamic_cast<ConstrainingEdge*>(edgeParent);
								if(ce){
									ce->add_constrained_object(cde);
									cde->set_constraining_object(ce);
								}
							}break;
						case FACE:{
							//	if the parent is a face, then it the face should have been converted
							//	to a constraining face during the face-coarsening step. If it hasn't
							//	been converted, then the somethings wrong with the selection states.
								ConstrainingFace* cf = dynamic_cast<ConstrainingFace*>(parent);
								UG_ASSERT(cf, "The parent face should have been transformed"
										  " to a constraining face already: "
										  << ElementDebugInfo(mg, cf));
								cf->add_constrained_object(cde);
								cde->set_constraining_object(cf);
								continue;
							}break;
						default:
							UG_THROW("Unsupported parent type during partial coarsening: "
									 << ElementDebugInfo(mg, parent));
							break;
					}
				}
				else{
				//	since no parent is present on the local process, we have to
				//	manually set the type of the constraining element.
					cde->set_parent_base_object_id(mg.parent_type(cde));
				}
			}break;

			case HNCM_REPLACE:{
			//	this should only not be set on constrained edges
				UG_ASSERT(!elem->is_constrained(), "Edge should not be constrained: "
						  << ElementDebugInfo(mg, elem));
				if(elem->is_constraining()){
				//	the constraining edge has to be transformed to a normal edge
					#ifdef UG_DEBUG
					//	make sure that the associated constrained vertex will be removed
						ConstrainingEdge* ce = dynamic_cast<ConstrainingEdge*>(elem);
						if(ce->num_constrained_vertices() > 0){
							VertexBase* hv = ce->constrained_vertex(0);
							UG_ASSERT(sel.get_selection_status(hv) == HNCM_ALL,
									  "A constraining edge may only be replaced by a normal edge"
									  " during coarsening, if its constrained vertex will be removed."
									  " At: " << ElementDebugInfo(mg, elem));
						}
					#endif

					mg.create_and_replace<Edge>(elem);
				}
				else{
				//	transform the edge to a constraining edge
					ConstrainingEdge* cge = *mg.create_and_replace<ConstrainingEdge>(elem);
				//	make sure that a connection with constrained children is established
				//	note - associated constrained vertices and constrained edges can
				//	not exist at this point.
					for(size_t i = 0; i < mg.num_children<EdgeBase>(cge); ++i){
						if(ConstrainedEdge* cde =
						   dynamic_cast<ConstrainedEdge*>(mg.get_child<EdgeBase>(cge, i)))
						{
							cge->add_constrained_object(cde);
							cde->set_constraining_object(cge);
						}
					}
				}
			}break;

			default:
				UG_THROW("Bad selection mark on edge during "
						 "HangingNodeRefiner_MultiGrid::coarsen. Internal error: "
						 << ElementDebugInfo(mg, elem));
				break;
		}
	}

////////////////////
//	VERTICES
	for(selector_t::traits<VertexBase>::iterator iter = sel.begin<VertexBase>();
		iter != sel.end<VertexBase>();)
	{
		VertexBase* elem = *iter;
		++iter;
		ISelector::status_t selState = sel.get_selection_status(elem);
		sel.deselect(elem);

		switch(selState){
			case RM_COARSEN:
			case HNCM_ALL:{
				mg.erase(elem);
			}break;

			case HNCM_PARTIAL:{
				if(elem->is_constrained()){
				//	a constrained element where not all nbrs are selected
				//	will be kept as it is.
					continue;
				}

			//	a constrained vertex is only created if the parent is an edge or a face
				ConstrainedVertex* hv = NULL;
				char parentType = mg.parent_type(elem);
				if((parentType == EDGE) || (parentType == FACE)){
					hv = *mg.create_and_replace<ConstrainedVertex>(elem);
					elem = hv;
				}
				else
					break;

			//	associated constrained and constraining elements
				GeometricObject* parent = mg.get_parent(hv);

				if(parent){
				//	the parent may not be a constrained object
					UG_ASSERT(!parent->is_constrained(),
							  "An element with a constrained parent may not be marked"
							  " for partial coarsening: " << ElementDebugInfo(mg, parent));

				//	if the parent is an edge or a face, the vertex will be converted
				//	to a hanging-vertex
					switch(parent->base_object_id()){
						case EDGE:{
						//	the parent already has to be a constraining edge, due to
						//	the operations performed during edge-coarsening
							ConstrainingEdge* ce = dynamic_cast<ConstrainingEdge*>(parent);
							UG_ASSERT(ce, "The parent edge already has to be constraining: "
										  << ElementDebugInfo(mg, parent));
							UG_ASSERT(ce->num_constrained_vertices() == 0,
									  "The parent edge may not yet have any constrained vertices: "
									  << ElementDebugInfo(mg, parent));

							ce->add_constrained_object(hv);
							hv->set_constraining_object(ce);
							hv->set_local_coordinate_1(0.5);
						}break;

						case FACE:{
						//	the parent already has to be a constraining face, due to
						//	the operations performed during face-coarsening
							ConstrainingFace* cf = dynamic_cast<ConstrainingFace*>(parent);
							UG_ASSERT(cf, "The parent face already has to be constraining: "
										  << ElementDebugInfo(mg, parent));
							UG_ASSERT(cf->num_constrained_vertices() == 0,
									  "The parent face may not yet have any constrained vertices: "
									  << ElementDebugInfo(mg, parent));

							cf->add_constrained_object(hv);
							hv->set_constraining_object(cf);
							hv->set_local_coordinates(0.5, 0.5);
						}break;

						default:
							UG_THROW("Unsupported parent type during partial"
									 " vertex coarsening: " << ElementDebugInfo(mg, parent));
							break;
					}
				}
				else{
				//	There is no parent on the local process. We thus only set the
				//	type of the constraining element
					hv->set_parent_base_object_id(mg.parent_type(hv));
				}

			}break;

			default:
				UG_THROW("Bad selection mark on vertex during "
						 "HangingNodeRefiner_MultiGrid::coarsen. Internal error: "
						 << ElementDebugInfo(mg, elem));
				break;
		}
	}

	debug_save(sel, "coarsen_marks_06_coarsening_done");
	if(debugging_enabled()){
		ParallelLayoutDebugSave(mg);
	}

	//ContinuousDebugSave(sel);
//	inform derived classes that coarsening is done
	post_coarsen();
	clear_marks();
	m_messageHub->post_message(GridMessage_Adaption(GMAT_HNODE_COARSENING_ENDS));

	return true;
}

}// end of namespace

