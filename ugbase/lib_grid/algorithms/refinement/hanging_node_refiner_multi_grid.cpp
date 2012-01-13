// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 11.01.2011 (m,d,y)
 
#include "common/assert.h"
#include "lib_grid/algorithms/multi_grid_util.h"
#include "lib_grid/algorithms/selection_util.h"
#include "hanging_node_refiner_multi_grid.h"

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

}

HangingNodeRefiner_MultiGrid::
HangingNodeRefiner_MultiGrid(MultiGrid& mg,
						IRefinementCallback* refCallback) :
	HangingNodeRefinerBase(refCallback),
	m_pMG(NULL)
{
	set_grid(&mg);
}

HangingNodeRefiner_MultiGrid::
~HangingNodeRefiner_MultiGrid()
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
			|| ConstrainingEdge::type_match(elem);
}

bool
HangingNodeRefiner_MultiGrid::
refinement_is_allowed(Face* elem)
{
	return (!m_pMG->has_children(elem))
			|| ConstrainingFace::type_match(elem);
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
collect_objects_for_refine()
{
	if(!m_pMG)
		throw(UGError("HangingNodeRefiner_MultiGrid::collect_objects_for_refine: No grid assigned."));

	MultiGrid& mg = *m_pMG;

//	first call base implementation
	HangingNodeRefinerBase::collect_objects_for_refine();

//	now select all associated vertices of marked objects,
//	since we have to create new vertices in the next levels of the hierarchies.
//	only vertices which do not already have child vertices are selected.
	for(EdgeBaseIterator iter = m_selMarkedElements.begin<EdgeBase>();
		iter != m_selMarkedElements.end<EdgeBase>(); ++iter)
	{
		for(size_t i = 0; i < (*iter)->num_vertices(); ++i){
			if(!mg.has_children((*iter)->vertex(i)))
				mark((*iter)->vertex(i));
		}
	}

	for(FaceIterator iter = m_selMarkedElements.begin<Face>();
		iter != m_selMarkedElements.end<Face>(); ++iter)
	{
		for(size_t i = 0; i < (*iter)->num_vertices(); ++i){
			if(!mg.has_children((*iter)->vertex(i)))
				mark((*iter)->vertex(i));
		}
	}

	for(VolumeIterator iter = m_selMarkedElements.begin<Volume>();
		iter != m_selMarkedElements.end<Volume>(); ++iter)
	{
		for(size_t i = 0; i < (*iter)->num_vertices(); ++i){
			if(!mg.has_children((*iter)->vertex(i)))
				mark((*iter)->vertex(i));
		}
	}
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

void HangingNodeRefiner_MultiGrid::
refine_constraining_edge(ConstrainingEdge* cge)
{
//	call base implementation to perform refinement
	HangingNodeRefinerBase::refine_constraining_edge(cge);

//	the constrained edge is now a normal edge
	m_pMG->create_and_replace<Edge>(cge);
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
		UG_ASSERT(mg.num_children<TBaseElem>(e) > 0,
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
		}
		else{
		//	the parent of this element has a different type or it doesn't exist
		//	at all. The element is thus not regarded as a valid candidate.
			sel.deselect(e);
		}
	}
	mg.end_marking();
}

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
				if(get_mark(nbrs[i]) == RM_COARSEN)
					allSelected = false;
			}

			if(!allSelected)
				sel.deselect(e);
		}
	}
}

void HangingNodeRefiner_MultiGrid::
collect_objects_for_coarsen()
{
//	first we'll shrink the selection so that only surface elements are selected
	restrict_selection_to_surface_coarsen_elements<Volume>();
	restrict_selection_to_surface_coarsen_elements<Face>();
	restrict_selection_to_surface_coarsen_elements<EdgeBase>();
	restrict_selection_to_surface_coarsen_elements<VertexBase>();

//	restrict to coarsen family
	restrict_selection_to_coarsen_families<Volume>();
	restrict_selection_to_coarsen_families<Face>();
	restrict_selection_to_coarsen_families<EdgeBase>();
	restrict_selection_to_coarsen_families<VertexBase>();

//	now select all associated geometric objects. This is required so that
//	elements in between hihgher dimensional coarsen elements will be deleted,
//	too.
	SelectAssociatedGeometricObjects(get_refmark_selector(), RM_COARSEN);

//todo:	This would be a good time to communicate between processes, so that
//		all processes know, whether constrained / constraining edges are
//		to be coarsened and whether normal edges will become constrained
//		edges after the coarsening.

//	deselect side elements which connect to unselected elements, since those
//	may not be deleted.
	adjust_coarsen_marks_on_side_elements<Face>();
	adjust_coarsen_marks_on_side_elements<EdgeBase>();
	adjust_coarsen_marks_on_side_elements<VertexBase>();

//todo: check whether sides are selected or whether their children are selected.
//		if not, whether its children are selected. If not, the whole family
//		may not be refined.
//		Note that this can be quite difficult. Probably one has to move through
//		the levels from top-1 to bottom.
}

bool HangingNodeRefiner_MultiGrid::
coarsen()
{
	if(!m_pMG)
		return false;

	collect_objects_for_coarsen();

//todo: erase elements and introduce constraining/constrained elements as required.

//	...

//todo:	return true
	return false;
}

}// end of namespace
