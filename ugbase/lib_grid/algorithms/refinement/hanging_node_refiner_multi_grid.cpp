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
	HangingNodeRefinerBase(refCallback),
	m_pMG(NULL)
{
	add_ref_mark_adjuster(MGHNodeAdjuster::create());
}

HangingNodeRefiner_MultiGrid::
HangingNodeRefiner_MultiGrid(MultiGrid& mg,
						IRefinementCallback* refCallback) :
	HangingNodeRefinerBase(refCallback),
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
	HangingNodeRefinerBase::set_grid(mg);

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
		Selector& sel = get_refmark_selector();

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
	for(VertexBaseIterator iter = m_selMarkedElements.begin<VertexBase>();
		iter != m_selMarkedElements.end<VertexBase>(); ++iter)
	{
		UG_ASSERT(!mg.has_children(*iter), "Only vertices without children should be selected.");
		if(refinement_is_allowed(*iter)){
			VertexBase* vrt = *mg.create<Vertex>(*iter);
			if(m_refCallback)
				m_refCallback->new_vertex(vrt, *iter);
		}
	}
}

void HangingNodeRefiner_MultiGrid::
refine_constraining_edge(ConstrainingEdge* cge)
{
//	call base implementation to perform refinement
	HangingNodeRefinerBase::refine_constraining_edge(cge);

//	the constrained edge is now a normal edge
	m_pMG->create_and_replace<Edge>(cge);
}

void HangingNodeRefiner_MultiGrid::
refine_edge_with_normal_vertex(EdgeBase* e, VertexBase** newCornerVrts)
{
//	collect child corners
	std::vector<VertexBase*> childCorners;
	collect_child_corners(childCorners, e);

//	call original implementation
	HangingNodeRefinerBase::refine_edge_with_normal_vertex(e, &childCorners.front());
}

void HangingNodeRefiner_MultiGrid::
refine_edge_with_hanging_vertex(EdgeBase* e, VertexBase** newCornerVrts)
{
//	collect child corners
	std::vector<VertexBase*> childCorners;
	collect_child_corners(childCorners, e);

//	call original implementation
	HangingNodeRefinerBase::refine_edge_with_hanging_vertex(e, &childCorners.front());
}

void HangingNodeRefiner_MultiGrid::
refine_face_with_normal_vertex(Face* f, VertexBase** newCornerVrts)
{
//	collect child corners
	std::vector<VertexBase*> childCorners;
	collect_child_corners(childCorners, f);

//	call original implementation
	HangingNodeRefinerBase::refine_face_with_normal_vertex(f, &childCorners.front());
}

void HangingNodeRefiner_MultiGrid::
refine_face_with_hanging_vertex(Face* f, VertexBase** newCornerVrts)
{
//	collect child corners
	std::vector<VertexBase*> childCorners;
	collect_child_corners(childCorners, f);

//	call original implementation
	HangingNodeRefinerBase::refine_face_with_hanging_vertex(f, &childCorners.front());
}

void HangingNodeRefiner_MultiGrid::
refine_volume_with_normal_vertex(Volume* v, VertexBase** newCornerVrts)
{
//	collect child corners
	std::vector<VertexBase*> childCorners;
	collect_child_corners(childCorners, v);

//	call original implementation
	HangingNodeRefinerBase::refine_volume_with_normal_vertex(v, &childCorners.front());
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


/*
template <class TElem>
void HangingNodeRefiner_MultiGrid::
select_subsurface_elements_only()
{
//todo:	This algorithm is too conservative. IsSubSurfaceElement(..., true)
//		causes that many possible coarsening possibilities are being ignored.
//		IsSubSurfaceElement should also take a callback 'cbConsiderElem', which
//		allows to ignore all marked surface elements.
//		In a first step one then should select all surface and sub-surface elements
//		(use IsSubSurfaceElement(..., false)) and then use the version with
//		IsSubSurfaceElement(..., true) and the mentioned callback.

	MultiGrid& mg = *m_pMG;
	Selector& sel = get_refmark_selector();

	typedef typename geometry_traits<TElem>::iterator ElemIter;
	typedef typename TElem::geometric_base_object TBaseElem;

	bool removedRefinementMark = false;

//	iterate over all selected elements of this type.
//	Note that we only can advance the iterator on special occasions.
//	Also note, that there are some continue commands in this loop.
	for(ElemIter iter = sel.begin<TElem>(); iter != sel.end<TElem>();)
	{
		TBaseElem* e = *iter;

	//	remove all refinement marks
		if(get_mark(e) != RM_COARSEN){
			++iter;
			sel.deselect(e);
			removedRefinementMark = true;
			continue;
		}

	//	check whether e is in the surface level
		int numChildren = mg.num_children<TBaseElem>(e);
		if(numChildren == 0){
		//	ok it is. if it has no parent of the same type, it can't be coarsened.
		//	At least not directly
			TBaseElem* parent = dynamic_cast<TBaseElem*>(mg.get_parent(e));
			if(parent != NULL){
			//	If the parent is not already marked, we'll
			//	perform some checks to see, whether it can be coarsened.
				if(!is_marked(parent)){
				//	check whether the parent is a sub-surface element
					if(IsSubSurfaceElement(mg, parent, true)){
					//	check whether all children of the parent are marked for coarsening
						bool markParent = true;
						size_t numSiblings = mg.num_children<TBaseElem>(parent);
						for(size_t i_child = 0; i_child < numSiblings; ++i_child){
							TBaseElem* child = mg.get_child<TBaseElem>(parent, i_child);
							if(get_mark(child) != RM_COARSEN)
							{
								markParent = false;
								break;
							}
						}

						if(markParent){
						//	all are marked surface elements.
						//	The parent thus can be coarsened.
							mark(parent, RM_COARSEN);
						}
					}
				}
			}

		//	we have to deselect the surface element. First advance the iterator
		//	then perform deselection
			++iter;
			sel.deselect(e);
			continue;

		}
		else{
		//	this element does not lie on the surface. If it is not a sub-surface
		//	element, we will deselect it.
			if(!IsSubSurfaceElement(mg, e, true)){
				++iter;
				sel.deselect(e);
				continue;
			}
		}

	//	advance the iterator (*iter at this point should be a sub-surface element
		++iter;
	}

	if(removedRefinementMark){
		UG_LOG("WARNING in HangingNodeRefiner_MultiGrid::collect_objects_for_coarsen: "
				"Removed refinement marks!\n");
	}
}
*/

template <class TElem>
void HangingNodeRefiner_MultiGrid::
restrict_selection_to_surface_coarsen_elements()
{
	MultiGrid& mg = *m_pMG;
	Selector& sel = get_refmark_selector();

	for(typename Selector::traits<TElem>::iterator iter = sel.begin<TElem>();
		iter != sel.end<TElem>();)
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

template <class TElem>
void HangingNodeRefiner_MultiGrid::
restrict_selection_to_coarsen_families()
{
	typedef typename TElem::geometric_base_object	TBaseElem;
	MultiGrid& mg = *m_pMG;
	Selector& sel = get_refmark_selector();

//	coarsening is only performed for complete families. If some siblings aren't
//	selected for coarsening, the whole family may not be coarsened. All siblings
//	have to be deselected in this case.

//	mark parents, which have been processed and mark siblings, which are valid
//	candidates. If thus a surface element is encountered, which has a marked
//	parent, but is not marked itself, then it is clear that it is not a valid
//	candidate.
	mg.begin_marking();
	for(typename Selector::traits<TElem>::iterator iter = sel.begin<TElem>();
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


/*
template <class TElem>
void HangingNodeRefiner_MultiGrid::
adjust_coarsen_marks_on_side_elements()
{
	if(!TElem::CAN_BE_SIDE)
		return;

	typedef typename TElem::geometric_base_object	TBaseElem;
	MultiGrid& mg = *m_pMG;
	Selector& sel = get_refmark_selector();

	std::vector<typename TElem::sideof*> nbrs;

//	if the grid doesn't contain any elements of higher dimension, then we're
//	already done
	if(mg.num<typename TElem::sideof>() > 0){
		for(typename Selector::traits<TElem>::iterator iter = sel.begin<TElem>();
			iter != sel.end<TElem>();)
		{
			TElem* e = *iter;
			++iter;

			CollectAssociated(nbrs, mg, e);

		//	check whether all nbrs are selected
			bool allSelected = true;
			for(size_t i = 0; i < nbrs.size(); ++i){
				if(get_mark(nbrs[i]) != RM_COARSEN){
					allSelected = false;
					break;
				}
			}

			if(!allSelected){
			//	if the element is constraining or constrained, then we'll deselect it.
			//	if it doesn't have a parent, then we'll deselect it too.
				if((!mg.has_parent(e)) || e->is_constraining() || e->is_constrained())
					sel.deselect(e);
				else{
				//	apply a mark that the element shall be converted to a constrained edge
					mark_for_hnode_coarsening(e);
				}
			}
		}
	}
}
*/
/*
template<class TElem>
static void DeselectFamily(ISelector& sel, MultiGrid& mg, TElem* elem)
{
	typedef typename TElem::geometric_base_object	TBaseElem;
	sel.deselect(elem);
	size_t numChildren = mg.num_children<TBaseElem>(elem);

	for(size_t i = 0; i < numChildren; ++i)
		DeselectFamily(sel, mg, mg.get_child<TBaseElem>(elem, i));

//	here we also indirectly take care of low-dimensional children of elem.
	for(size_t i = 0; i < elem->num_sides(); ++i)
		DeselectFamily(sel, mg, mg.get_side(elem, i));
}
*/
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

template <class TElem>
bool HangingNodeRefiner_MultiGrid::
deselect_invalid_coarsen_families()
{
	typedef typename TElem::geometric_base_object	TBaseElem;
	MultiGrid& mg = *m_pMG;
	Selector& sel = get_refmark_selector();

//	finally deselect families, which contain an element connected to an unselected
//	constraining side. This is only important for volumes and faces.
//	we'll use marks to collect subsurface elements
	mg.begin_marking();

	vector<TBaseElem*> doomedFamilies;
	vector<typename TBaseElem::side*> sides;
	for(typename Selector::template traits<TElem>::iterator
		iter = sel.begin<TElem>(); iter != sel.end<TElem>(); ++iter)
	{
	//	ignore elements with mark HNCM_NO_NBRS_SELECTED, since they won't be
	//	coarsened anyways
		if(sel.get_selection_status(*iter) == HNCM_NO_NBRS_SELECTED)
			continue;

	//	check if the parent has the same base type
		TBaseElem* parent = dynamic_cast<TBaseElem*>(mg.get_parent(*iter));
		if(parent){
			if(mg.is_marked(parent))
				continue;

			CollectAssociated(sides, mg, *iter);
		//	check if a side is a constraining object
			for(size_t i = 0; i < sides.size(); ++i){
				if(!sel.is_selected(sides[i])){
				//	the family is doomed! We have to deselect it.
					doomedFamilies.push_back(parent);
					mg.mark(parent);
					break;
				}
			}
		}
	}

//	now deselect all doomed families
	for(size_t i = 0; i < doomedFamilies.size(); ++i){
		DeselectFamily(sel, mg, doomedFamilies[i]);
	}

	mg.end_marking();
	return !doomedFamilies.empty();
}

template <class TElem>
void HangingNodeRefiner_MultiGrid::
classify_selection()
{
	if(!TElem::CAN_BE_SIDE)
		return;

	typedef typename TElem::sideof	TSideOf;

	MultiGrid& mg = *m_pMG;
	Selector& sel = get_refmark_selector();

//	if there are no elements of higher dimension, then return immediatly.
	if(mg.num<TSideOf>() == 0){
		return;
	}

	std::vector<TSideOf*> nbrs;

//	iterate over all selected elements and assign the selection status based
//	on how many elements are selected in the direct neighborhood.
	for(typename Selector::traits<TElem>::iterator iter = sel.begin<TElem>();
		iter != sel.end<TElem>(); ++iter)
	{
		TElem* e = *iter;
		CollectAssociated(nbrs, mg, e);

		bool nbrsPartiallySelected = false;
		size_t numSel = 0;
		for(size_t i = 0; i < nbrs.size(); ++i){
			Selector::status_t s = sel.get_selection_status(nbrs[i]);
			//if((s == RM_COARSEN) || (s == HNCM_ALL_NBRS_SELECTED))
			if(sel.is_selected(nbrs[i]) && (s != HNCM_NO_NBRS_SELECTED)){
				++numSel;
				if(s == HNCM_SOME_NBRS_SELECTED)
					nbrsPartiallySelected = true;
			}
		}

		if((numSel == nbrs.size()) && (!nbrsPartiallySelected)) // if nbrs.size() == 0, thats fine
			sel.select(e, HNCM_ALL_NBRS_SELECTED);
		else if(numSel == 0)
			sel.select(e, HNCM_NO_NBRS_SELECTED);
		else
			sel.select(e, HNCM_SOME_NBRS_SELECTED);
	}
}

template <class TElem>
void HangingNodeRefiner_MultiGrid::
deselect_isolated_sides()
{
	if(!TElem::CAN_BE_SIDE)
		return;

	typedef typename TElem::sideof	TSideOf;

	MultiGrid& mg = *m_pMG;
	Selector& sel = get_refmark_selector();

//	if there are no elements of higher dimension, then the element can not be
//	considered as a side
	if(mg.num<TSideOf>() == 0){
		return;
	}

//	iterate over all selected elements of the given type and deselect it,
	for(typename Selector::traits<TElem>::iterator iter = sel.begin<TElem>();
		iter != sel.end<TElem>();)
	{
		TElem* e = *iter;
		++iter;

		if(sel.get_selection_status(e) == HNCM_NO_NBRS_SELECTED)
			sel.deselect(e);
	}
}

template <class TElem>
void HangingNodeRefiner_MultiGrid::
deselect_uncoarsenable_parents()
{
	typedef typename TElem::geometric_base_object	TBaseElem;
	MultiGrid& mg = *m_pMG;
	Selector& sel = get_refmark_selector();

	for(typename Selector::traits<TElem>::iterator iter = sel.begin<TElem>();
		iter != sel.end<TElem>();)
	{
		TElem* e = *iter;
		++iter;

	//	if the element has children and if one of the children may not be
	//	coarsened (if HNCM_ALL_NBRS_SELECTED is not specified), then the
	//	element itself may not be coarsened, too.
		if(mg.has_children<TBaseElem>(e)){
			for(size_t i = 0; i < mg.num_children<TBaseElem>(e); ++i){
				TBaseElem* child = mg.get_child<TBaseElem>(e, i);
				if(sel.get_selection_status(child) != HNCM_ALL_NBRS_SELECTED){
					sel.deselect(e);
					break;
				}
			}
		}
	}
}


void HangingNodeRefiner_MultiGrid::
save_coarsen_marks_to_file(Selector& sel, const char* filename)
{
	assert(sel.grid());
	Grid& g = *sel.grid();

	SubsetHandler sh(g);

	AssignGridToSubset(g, sh, 4);
	for(VertexBaseIterator iter = g.vertices_begin(); iter != g.vertices_end(); ++iter){
		Selector::status_t status = sel.get_selection_status(*iter);
		switch(status){
			case RM_COARSEN: sh.assign_subset(*iter, 0); break;
			case HNCM_ALL_NBRS_SELECTED: sh.assign_subset(*iter, 1); break;
			case HNCM_SOME_NBRS_SELECTED: sh.assign_subset(*iter, 2); break;
			case HNCM_NO_NBRS_SELECTED: sh.assign_subset(*iter, 3); break;
		}
	}

	for(EdgeBaseIterator iter = g.edges_begin(); iter != g.edges_end(); ++iter){
		Selector::status_t status = sel.get_selection_status(*iter);
		switch(status){
			case RM_COARSEN: sh.assign_subset(*iter, 0); break;
			case HNCM_ALL_NBRS_SELECTED: sh.assign_subset(*iter, 1); break;
			case HNCM_SOME_NBRS_SELECTED: sh.assign_subset(*iter, 2); break;
			case HNCM_NO_NBRS_SELECTED: sh.assign_subset(*iter, 3); break;
		}
	}

	for(FaceIterator iter = g.faces_begin(); iter != g.faces_end(); ++iter){
		Selector::status_t status = sel.get_selection_status(*iter);
		switch(status){
			case RM_COARSEN: sh.assign_subset(*iter, 0); break;
			case HNCM_ALL_NBRS_SELECTED: sh.assign_subset(*iter, 1); break;
			case HNCM_SOME_NBRS_SELECTED: sh.assign_subset(*iter, 2); break;
			case HNCM_NO_NBRS_SELECTED: sh.assign_subset(*iter, 3); break;
		}
	}

	for(VolumeIterator iter = g.volumes_begin(); iter != g.volumes_end(); ++iter){
		Selector::status_t status = sel.get_selection_status(*iter);
		switch(status){
			case RM_COARSEN: sh.assign_subset(*iter, 0); break;
			case HNCM_ALL_NBRS_SELECTED: sh.assign_subset(*iter, 1); break;
			case HNCM_SOME_NBRS_SELECTED: sh.assign_subset(*iter, 2); break;
			case HNCM_NO_NBRS_SELECTED: sh.assign_subset(*iter, 3); break;
		}
	}

	AssignSubsetColors(sh);
	sh.subset_info(0).name = "coarsen";
	sh.subset_info(1).name = "all_nbrs_selected";
	sh.subset_info(2).name = "some_nbrs_selected";
	sh.subset_info(3).name = "no_nbrs_selected";
	sh.subset_info(4).name = "nothing selected";


	if(MultiGrid* mg = dynamic_cast<MultiGrid*>(&g))
		SaveGridHierarchyTransformed(*mg, sh, filename, 0.1);
	else
		SaveGridToFile(g, sh, filename);
}


///	temporary method, which will be removed when debugging is done.
void HangingNodeRefiner_MultiGrid::
debug_save(Selector& sel, const char* filename)
{
//	save_coarsen_marks_to_file(sel, filename);
}

void HangingNodeRefiner_MultiGrid::
collect_objects_for_coarsen()
{
/* This is where we want to go:
 * - normal marks on vertices, if all associated elements have normal mark
 *
 * - normal mark on normal edges if all associated faces are marked normal
 * - hnode mark on normal edges if not all (but at least one) associated face is marked normal
 * - normal mark on constrained edge if all associated faces are marked normal
 * - normal mark on constraining edge, if associated constrained are marked,
 *   but no associated faces are marked.
 * - hnode mark on constraining edge, if constrained edges are marked and at least
 *   one associated face is marked.
 *
 * - normal mark on normal faces, if all associated volumes have a mark,
 * - hnode mark on normal faces, if not all but at least one associated volume has a mark.
 * - normal mark on constrained face, if all associated volumes are marked.
 * - normal mark on constraining face, if no ass. volume is marked, but associated
 *   constrained faces are marked.
 * - hnode mark on constraining face, if at least one ass. vol is marked and constrained
 *   faces are marked.
 *
 * - normal marks on all volumes which shall be coarsened.
 *
 * One has to be carful, since only whole families may be coarsened at once
 * (either all children of a parent or none may be removed).
 *
 *
 * This is how we get there:
 * - Select surface coarsen families together with associated elements.
 *   Deselect all others at the same time.
 * - Classify selection by applying selection-marks enumerated in HNodeCoarsenMarks.
 *   This is a virtual method. Here parallel communication may be performed in
 *   derived classes.
 * - Deselect constrained objects as required (first faces, then edges, then vertices)
 * - Deselect constraining objects as required (first faces, then edges)
 * - Deselect coarsen families which have missing sides (first
 */

//	this selector holds all refmarks
	Selector& sel = get_refmark_selector();

//	store whether the grid contains edges, faces or volumes.
//	we explicitly do this, since a parallel version might have to communicate...
	bool containsEdges = contains_edges();
	bool containsFaces = contains_faces();
	bool containsVolumes = contains_volumes();

debug_save(sel, "debug_save_0_initial_selection.ugx");

//	first we'll shrink the selection so that only surface elements are selected
	restrict_selection_to_surface_coarsen_elements<Volume>();
	restrict_selection_to_surface_coarsen_elements<Face>();
	restrict_selection_to_surface_coarsen_elements<EdgeBase>();
	restrict_selection_to_surface_coarsen_elements<VertexBase>();

debug_save(sel, "debug_save_1_restricted_to_surface_elems.ugx");

//	restrict to coarsen family
	restrict_selection_to_coarsen_families<Volume>();
	restrict_selection_to_coarsen_families<Face>();
	restrict_selection_to_coarsen_families<EdgeBase>();
	restrict_selection_to_coarsen_families<VertexBase>();

debug_save(sel, "debug_save_2_restricted_to_coarsen_families.ugx");

//	now select all associated geometric objects. This is required so that
//	elements in between hihgher dimensional coarsen elements will be deleted,
//	too.
	SelectAssociatedGeometricObjects(sel, RM_COARSEN);

//debug_save(sel, "debug_save_3_associated_objects_selected.ugx");

//todo:	This iteration is due to a lazy implementation: deselection of
//		invalid coarsen families can lead to a recursion. One could handle
//		such a case using queues with only involved elements. However - this
//		implementation is a lot less programming work and hopefully shouldn't
//		lead to significant performance losses. It also should lend itself
//		very well to parallelization.

	bool bRunning = true;
	while(bRunning){
		//UG_LOG("classifying selection...\n");
	//	we now have to classify the selection to mark, whether a selected element
	//	is surrounded by only marked, by some marked or by no marked elements of
	//	next higher dimension.
	//	we then deselect elements which do not have selected associated elements of
	//	the next higher dimension.
	//	note that we do this for faces first, then for edges and finally for vertices.
		if(containsVolumes){
			classify_selection<Face>();
		//	broadcast face marks
			broadcast_face_coarsen_marks();

			deselect_isolated_sides<Face>();
			deselect_uncoarsenable_parents<Face>();

		// If a constrained element is selected with HNCM_ALL_NBRS_SELECTED and its associated
		// constraining element is not selected, then we have to convert the associated
		// constraining element to a normal element. we'll thus select it with
		// HNCM_NO_NBRS_SELECTED.
		// Note that it is sufficient if we know that one constrained element is marked with
		// HNCM_ALL_NBRS_SELECTED, since the other one should have the same mark in this case
		// (since we always select whole families...).
			for(Selector::traits<ConstrainedTriangle>::iterator iter = sel.begin<ConstrainedTriangle>();
				iter != sel.end<ConstrainedTriangle>(); ++iter)
			{
				ConstrainedFace* cf = *iter;
				if(sel.get_selection_status(cf) == HNCM_ALL_NBRS_SELECTED){
					if(!sel.is_selected(cf->get_constraining_object()))
						sel.select(cf->get_constraining_object(), HNCM_NO_NBRS_SELECTED);
				}
			}

			for(Selector::traits<ConstrainedQuadrilateral>::iterator iter = sel.begin<ConstrainedQuadrilateral>();
				iter != sel.end<ConstrainedQuadrilateral>(); ++iter)
			{
				ConstrainedFace* cf = *iter;
				if(sel.get_selection_status(cf) == HNCM_ALL_NBRS_SELECTED){
					if(!sel.is_selected(cf->get_constraining_object()))
						sel.select(cf->get_constraining_object(), HNCM_NO_NBRS_SELECTED);
				}
			}
		}

	//	handle edges
		if(containsFaces){
			classify_selection<EdgeBase>();
			broadcast_edge_coarsen_marks();

			debug_save(sel, "debug_save_3a_classified_selection.ugx");
			deselect_isolated_sides<EdgeBase>();
			debug_save(sel, "debug_save_3b_deleted_isolated_sides.ugx");
			deselect_uncoarsenable_parents<EdgeBase>();
			debug_save(sel, "debug_save_3c_deleted_uncoarsenable_parents.ugx");

			// see commentary for constrained traingles / quadrilaterals above
			for(Selector::traits<ConstrainedEdge>::iterator iter = sel.begin<ConstrainedEdge>();
				iter != sel.end<ConstrainedEdge>(); ++iter)
			{
				ConstrainedEdge* ce = *iter;
				if(sel.get_selection_status(ce) == HNCM_ALL_NBRS_SELECTED){
					if(!sel.is_selected(ce->get_constraining_object()))
						sel.select(ce->get_constraining_object(), HNCM_NO_NBRS_SELECTED);
				}
			}
		}

	//	handle vertices
		if(containsEdges){
			classify_selection<VertexBase>();
			broadcast_vertex_coarsen_marks();

			deselect_isolated_sides<VertexBase>();
			deselect_uncoarsenable_parents<VertexBase>();
		}

		{
			static int counter = 0;
			stringstream ss;
			ss << "debug_save_4_adjusted_side_selections" << counter << ".ugx";
			++counter;
			debug_save(sel, ss.str().c_str());
		}

	//	now deselect coarsen families, which are adjacent to unselected constraining sides
	//	Those may not be coarsened and are regarded as invalid.
		bool deselectedFamilies = deselect_invalid_coarsen_families<Face>();
		{
			static int counter = 0;
			stringstream ss;
			ss << "debug_save_5_deselected_invalid_coarsen_families_face" << counter << ".ugx";
			++counter;
			debug_save(sel, ss.str().c_str());
		}
		deselectedFamilies |= deselect_invalid_coarsen_families<Volume>();

		{
			static int counter = 0;
			stringstream ss;
			ss << "debug_save_6_deselected_invalid_coarsen_families_vol" << counter << ".ugx";
			++counter;
			debug_save(sel, ss.str().c_str());
		}

	//	this allows us to check whether another iteration is enforced by
	//	this or by any other process (important in a parallel environment).
		bRunning = continue_collect_objects_for_coarsen(deselectedFamilies);
	}
	debug_save(sel, "debug_save_7_final_selection.ugx");
}

//void HangingNodeRefiner_MultiGrid::
//assign_hnode_coarsen_marks()
//{
//	MultiGrid& mg = *m_pMG;
//	Selector& sel = get_refmark_selector();
//////////////////////
////	FACES
//	if(mg.num<Volume>() > 0){
//		for(Selector::traits<Face>::iterator iter = sel.begin<Face>();
//			iter != sel.end<Face>();)
//		{
//			Face* elem = *iter;
//			if((sel.get_selection_status(elem) == HNCM_SOME_NBRS_SELECTED)
//				&& !elem->is_constrained())
//			{
//				Face* parent = dynamic_cast<Face*>(mg.get_parent(elem));
//				if(parent){
//					UG_ASSERT(!parent->is_constrained(), "Parent may not be constrained");
//					mark_for_hnode_refinement(elem);
//					if(!parent->is_constraining())
//						mark_for_hnode_refinement(parent);
//				}
//				else{
//				//todo: what if parent is on another process (elem is a v-slave)?
//					if(elem->is_constraining()){
//						mark_for_hnode_refinement(elem);
//					}
//					continue;
//				}
//			}
//		}
//	}
//
//////////////////////
////	FACES
//	if(mg.num<Face>() > 0){
//		for(Selector::traits<EdgeBase>::iterator iter = sel.begin<EdgeBase>();
//			iter != sel.end<EdgeBase>();)
//		{
//			EdgeBase* elem = *iter;
//			if((sel.get_selection_status(elem) == HNCM_SOME_NBRS_SELECTED)
//				&& !elem->is_constrained())
//			{
//				GeometricObject* parent = mg.get_parent(elem);
//				if(parent){
//					UG_ASSERT(parent->base_object_id() != VOLUME,
//							"child-edges of volumes may not be marked for partial coarsening");
//					UG_ASSERT(!parent->is_constrained(), "Parent may not be constrained");
//
//					mark_for_hnode_refinement(elem);
//					if(!parent->is_constraining())
//						mark_for_hnode_refinement(parent);
//				}
//				else{
//				//todo: what if parent is on another process (elem is a v-slave)?
//					if(elem->is_constraining()){
//						mark_for_hnode_refinement(elem);
//					}
//					continue;
//				}
//			}
//		}
//	}
//}

bool HangingNodeRefiner_MultiGrid::
perform_coarsening()
{
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
	Selector& sel = get_refmark_selector();

	collect_objects_for_coarsen();
//	assign_hnode_coarsen_marks();

//	if a debug file was specified, we'll now save the marks to that file
	if(!m_adjustedMarksDebugFilename.empty())
		save_coarsen_marks_to_file(sel, m_adjustedMarksDebugFilename.c_str());

//	if the selector is empty, we'll return false
	if(sel.empty())
		return false;

//	inform derived classes that coarsening begins
	pre_coarsen();

////////////////////
//	VOLUMES
//	erase elements and introduce constraining/constrained elements as required.
	mg.erase(sel.begin<Volume>(), sel.end<Volume>());

////////////////////
//	FACES
	for(Selector::traits<Face>::iterator iter = sel.begin<Face>();
		iter != sel.end<Face>();)
	{
		Face* elem = *iter;
		++iter;
		ISelector::status_t selState = sel.get_selection_status(elem);
		sel.deselect(elem);

		switch(selState){
			case RM_COARSEN:
			case HNCM_ALL_NBRS_SELECTED:{
			//	this should only be set on normal or constrained faces
				assert(!elem->is_constraining());
				mg.erase(elem);
			}break;

			case HNCM_SOME_NBRS_SELECTED:{
				if(elem->is_constrained()){
				//	a constrained element where not all nbrs are selected
				//	will be kept as it is.
					continue;
				}

			//	the face has to be replaced by a constrained face.
			//	Volume-children don't have to be considered here, since
			//	they should all be marked with HNCM_ALL_NBRS_SELECTED
				Face* parent = dynamic_cast<Face*>(mg.get_parent(elem));

				if(!parent){
				//	the face doesn't have a parent and thus can't be constrained.
				//	convert it to a normal face instead.
					if(elem->is_constraining()){
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

				//	we're done for this element
					continue;
				}

			//	the parent may not be a constrained object
				assert(!parent->is_constrained());

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

			//	elem is no longer valid!!!
				elem = NULL;

			//	check whether the parent is already constraining
				ConstrainingFace* cf = dynamic_cast<ConstrainingFace*>(parent);
				if(!cf){
				//	no it is not. We first have to transform it to a constraining face
					switch(parent->reference_object_id()){
						case ROID_TRIANGLE:
							cf = *mg.create_and_replace<ConstrainingTriangle>(parent);
							break;
						case ROID_QUADRILATERAL:
							cf = *mg.create_and_replace<ConstrainingQuadrilateral>(parent);
							break;
						default:
							UG_THROW("Unknown face reference object type in "
									 "HangingNodeRefiner_MultiGrid::coarsen");
							break;
					}
				}

				cf->add_constrained_object(cdf);
				cdf->set_constraining_object(cf);
			}break;

			case HNCM_NO_NBRS_SELECTED:{
			//	this should only be set on constraining faces
				assert(elem->is_constraining());

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
			}break;

			default:
				UG_THROW("Bad selection mark on face during "
						 "HangingNodeRefiner_MultiGrid::coarsen. Internal error.");
				break;
		}
	}

////////////////////
//	EDGES
	for(Selector::traits<EdgeBase>::iterator iter = sel.begin<EdgeBase>();
		iter != sel.end<EdgeBase>();)
	{
		EdgeBase* elem = *iter;
		++iter;
		ISelector::status_t selState = sel.get_selection_status(elem);
		sel.deselect(elem);

		switch(selState){
			case RM_COARSEN:
			case HNCM_ALL_NBRS_SELECTED:{
			//	this should only be set on normal or constrained edges
				assert(!elem->is_constraining());
				mg.erase(elem);
			}break;

			case HNCM_SOME_NBRS_SELECTED:{
				if(elem->is_constrained()){
				//	a constrained element where not all nbrs are selected
				//	will be kept as it is.
					continue;
				}

			//	the edge has to be replaced by a constrained edge.
			//	Volume-children don't have to be considered here, since
			//	they should all be marked with HNCM_ALL_NBRS_SELECTED
				GeometricObject* parent = mg.get_parent(elem);
				if(!parent){
				//	the edge doesn't have a parent and thus can't be constrained.
				//	convert it to a normal edge instead.
					if(elem->is_constraining())
						mg.create_and_replace<Edge>(elem);

				//	we're done for this element
					continue;
				}

			//	the parent may not be a constrained object
				assert(!parent->is_constrained());

			//	replace elem by a constrained edge
				ConstrainedEdge* cde = *mg.create_and_replace<ConstrainedEdge>(elem);

			//	elem is no longer valid!!!
				elem = NULL;

			//	if the parent is a face, then it the face should have been converted
			//	to a constraining face during the face-coarsening step. If it hasn't
			//	been converted, then the somethings wrong with the selection states.
				if(parent->base_object_id() == FACE){
					ConstrainingFace* cf = dynamic_cast<ConstrainingFace*>(parent);
					assert(cf);

					cf->add_constrained_object(cde);
					cde->set_constraining_object(cf);
					continue;
				}
				else{
					EdgeBase* edgeParent = dynamic_cast<EdgeBase*>(parent);
					assert(edgeParent);

				//	check whether the parent is already constraining
					ConstrainingEdge* ce = dynamic_cast<ConstrainingEdge*>(edgeParent);
					if(!ce){
					//	no it is not. We first have to transform it to a constraining edge
						ce = *mg.create_and_replace<ConstrainingEdge>(edgeParent);
					}

					ce->add_constrained_object(cde);
					cde->set_constraining_object(ce);
				}
			}break;

			case HNCM_NO_NBRS_SELECTED:{
			//	this should only be set on constraining edges
				assert(elem->is_constraining());

			//	the constraining edge has to be transformed to a normal edge
				mg.create_and_replace<Edge>(elem);
			}break;

			default:
				UG_THROW("Bad selection mark on edge during "
						 "HangingNodeRefiner_MultiGrid::coarsen. Internal error.");
				break;
		}
	}

////////////////////
//	VERTICES
	for(Selector::traits<VertexBase>::iterator iter = sel.begin<VertexBase>();
		iter != sel.end<VertexBase>();)
	{
		VertexBase* elem = *iter;
		++iter;
		ISelector::status_t selState = sel.get_selection_status(elem);
		sel.deselect(elem);

		switch(selState){
			case RM_COARSEN:
			case HNCM_ALL_NBRS_SELECTED:{
				mg.erase(elem);
			}break;

			case HNCM_SOME_NBRS_SELECTED:{
				if(elem->is_constrained()){
				//	a constrained element where not all nbrs are selected
				//	will be kept as it is.
					continue;
				}

			//	depending on its parent, the vertex may have to be replaced by
			//	a hanging vertex
				GeometricObject* parent = mg.get_parent(elem);

				if(!parent){
				//	the vertex doesn't have a parent and thus can't be constrained.
				//	we're done for this element
					continue;
				}

			//	the parent may not be a constrained object
				assert(!parent->is_constrained());

			//	if the parent is an edge or a face, the vertex will be converted
			//	to a hanging-vertex
				switch(parent->base_object_id()){
					case EDGE:{
					//	the parent already has to be a constraining edge, due to
					//	the operations performed during edge-coarsening
						ConstrainingEdge* ce = dynamic_cast<ConstrainingEdge*>(parent);
						assert(ce);
						assert(ce->num_constrained_vertices() == 0);

						ConstrainedVertex* hv = *mg.create_and_replace<ConstrainedVertex>(elem);
						ce->add_constrained_object(hv);
						hv->set_constraining_object(ce);
						hv->set_local_coordinate_1(0.5);
					}break;

					case FACE:{
					//	the parent already has to be a constraining face, due to
					//	the operations performed during face-coarsening
						ConstrainingFace* cf = dynamic_cast<ConstrainingFace*>(parent);
						assert(cf);
						assert(cf->num_constrained_vertices() == 0);

						ConstrainedVertex* hv = *mg.create_and_replace<ConstrainedVertex>(elem);
						cf->add_constrained_object(hv);
						hv->set_constraining_object(cf);
						hv->set_local_coordinates(0.5, 0.5);
					}break;

					default:
					//	child vertices of volumes lie inside of coarsen-families.
					//	HNCM_ALL_NBRS_SELECTED thus has to be applied.
						assert((parent->base_object_id() != VOLUME));
						break;
				}

			}break;

			default:
				UG_THROW("Bad selection mark on vertex during "
						 "HangingNodeRefiner_MultiGrid::coarsen. Internal error.");
				break;
		}
	}

//	inform derived classes that coarsening is done
	post_coarsen();

	return true;
}

/* Old adjust_selection code...
DebugSave(sel, "debug_save_0_initial_selection.ugx");

//	first we'll shrink the selection so that only surface elements are selected
	restrict_selection_to_surface_coarsen_elements<Volume>();
	restrict_selection_to_surface_coarsen_elements<Face>();
	restrict_selection_to_surface_coarsen_elements<EdgeBase>();
	restrict_selection_to_surface_coarsen_elements<VertexBase>();

DebugSave(sel, "debug_save_1_restricted_to_surface_elems.ugx");

//	restrict to coarsen family
	restrict_selection_to_coarsen_families<Volume>();
	restrict_selection_to_coarsen_families<Face>();
	restrict_selection_to_coarsen_families<EdgeBase>();
	restrict_selection_to_coarsen_families<VertexBase>();

DebugSave(sel, "debug_save_2_restricted_to_coarsen_families.ugx");

//	now select all associated geometric objects. This is required so that
//	elements in between hihgher dimensional coarsen elements will be deleted,
//	too.
	SelectAssociatedGeometricObjects(sel, RM_COARSEN);

DebugSave(sel, "debug_save_3_associated_objects_selected.ugx");


//	side elements which connect to unselected elements either have to be deselected
//	(if they are constrained / constraining objects) or have to be marked for
//	hnode coarsening (if they are normal objects).
	adjust_coarsen_marks_on_side_elements<Face>();
	adjust_coarsen_marks_on_side_elements<EdgeBase>();
	adjust_coarsen_marks_on_side_elements<VertexBase>();

DebugSave(sel, "debug_save_4_adjusted_coarsen_marks_on_side_elements.ugx");

//	now select constraining edges and faces of selected constrained ones.
//	we can select the constraining object without further checks, since the
//	preparations above assure, that either all its constrained objects are
//	selected or that none is.
//	All constraining objects selected here will have to be converted to constrained
//	ones later on. we thus apply a hnode mark
	for(Selector::traits<ConstrainedEdge>::iterator iter = sel.begin<ConstrainedEdge>();
		iter != sel.end<ConstrainedEdge>(); ++iter)
	{
		ConstrainedEdge* ce = *iter;
		sel.select(ce->get_constraining_object(), RM_COARSEN);
		mark_for_hnode_coarsening(ce->get_constraining_object());
	}

	for(Selector::traits<ConstrainedTriangle>::iterator iter = sel.begin<ConstrainedTriangle>();
		iter != sel.end<ConstrainedTriangle>(); ++iter)
	{
		ConstrainedFace* cf = *iter;
		sel.select(cf->get_constraining_object(), RM_COARSEN);
		mark_for_hnode_coarsening(cf->get_constraining_object());
	}

	for(Selector::traits<ConstrainedQuadrilateral>::iterator iter = sel.begin<ConstrainedQuadrilateral>();
		iter != sel.end<ConstrainedQuadrilateral>(); ++iter)
	{
		ConstrainedFace* cf = *iter;
		sel.select(cf->get_constraining_object(), RM_COARSEN);
		mark_for_hnode_coarsening(cf->get_constraining_object());
	}

DebugSave(sel, "debug_save_5_selected_constraining_objects.ugx");

//	now deselect coarsen families, which are adjacent to unselected constraining
//	elements. Those may not be coarsened and are regarded as invalid.
	deselect_invalid_coarsen_families<Volume>();
	deselect_invalid_coarsen_families<Face>();
DebugSave(sel, "debug_save_6_deselected_invalid_coarsen_families.ugx");
*/
}// end of namespace
