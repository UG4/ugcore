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
/*
template <class TElem>
void restrict_selection_to_surface_coarsen_elements()
{
	MultiGrid& mg = *m_pMG;
	Selector& sel = get_refmark_selector();

	for(typename Selector::traits<TElem>::iterator iter = sel.begin<TElem>();
		iter != sel.end<TElem>();)
	{
		TElem* e = *iter;
		++iter;

	//	make sure that only coarsen-marks are applied
		if(sel.get_mark(e) != RM_COARSEN){
			sel.deselect(e);
			continue;
		}

	//	make sure that the element is a surface element
		if(mg.has_children(e))
			sel.deselect(e);
	}
}

template <class TElem>
void deselect_undeletable_coarsen_elements()
{
	typedef typename TElem::GeometricBaseObject	TBaseElem;
	MultiGrid& mg = *m_pMG;
	Selector& sel = get_refmark_selector();

	vector<VertexBase*> vrts;
	vector<EdgeBase*> edges;
	vector<Face*> faces;
	vector<Volume*> vols;

	for(typename Selector::traits<TElem>::iterator iter = sel.begin<TElem>();
		iter != sel.end<TElem>();)
	{
		TElem* e = *iter;
		++iter;

		size_t numChildren = mg.num_children<TBaseElem>(e);
		if(numChildren > 0){
		//	we have to deselect the element
			sel.deselect(e);

		//	check if one of the children of the same base type is not selected
			bool gotOne = false;
			for(size_t i = 0; i < numChildren; ++i){
				if(!sel.is_selected(mg.get_child<TBaseElem>(e, i))){
					gotOne = true;
					break;
				}
			}

			if(gotOne){
			//	the element and all associated elements may not be erased.
			//	Note that also elements of higher dimension may not be erased,
			//	since otherwise neighbored surface elements could lie in
			//	levels further apart than 1.
				CollectAssociated(vrts, mg, e);
				CollectAssociated(edges, mg, e);
				CollectAssociated(faces, mg, e);
				CollectAssociated(vols, mg, e);

			//	Note that we'll deselect e a second time here. This is however
			//	not a problem
				for(size_t i = 0; i < vrts.size(); ++i)
					sel.deselect(vrts[i]);
				for(size_t i = 0; i < edges.size(); ++i)
					sel.deselect(edges[i]);
				for(size_t i = 0; i < faces.size(); ++i)
					sel.deselect(faces[i]);
				for(size_t i = 0; i < vols.size(); ++i)
					sel.deselect(vols[i]);
			}
		}
	}
}

template <class TElem>
void restrict_selection_to_coarsen_families()
{
	typedef typename TElem::GeometricBaseObject	TBaseElem;
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

	//	get the parent
		TBaseElem* parent = dynamic_cast<TBaseElem*>(mg.get_parent(e));
		if(parent){
		//	check whether all children of e of type TBaseElem are marked
		}
		else{
		//	the parent of this element has a different type or it doesn't exist
		//	at all. The element is thus not regarded as a valid candidate.
			sel.deselect(e);
		}
	}
	mg.end_marking();
}
*/
void HangingNodeRefiner_MultiGrid::
collect_objects_for_coarsen()
{
//	first we'll shrink the selection so that only surface elements are selected
/*	restrict_selection_to_surface_coarsen_elements<Volume>();
	restrict_selection_to_surface_coarsen_elements<Face>();
	restrict_selection_to_surface_coarsen_elements<EdgeBase>();
	restrict_selection_to_surface_coarsen_elements<VertexBase>();
*/
//	now select all associated geometric objects. Those will be needed, since
//	we have to check whether objects in an adaptive multigrid may be removed at all.
//	With this extended selection, the new surface level can be easily checked.
/*	SelectAssociatedGeometricObjects(get_refmark_selector(), RM_COARSEN);
*/
//	now we'll deselect all elements which have a children. If a child is unselected,
//	associated elements of higher dimension may not be removed, since this would lead
//	to a jump in the hierarchy higher than 1.


//	since an element may only be removed during coarsening, if all its siblings
//	are also removed, we'll now deselect all those, which may not be removed
/*	restrict_selection_to_valid_candidates<Volume>();
	restrict_selection_to_valid_candidates<Face>();
	restrict_selection_to_valid_candidates<EdgeBase>();
	// no need to check vertices
*/
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
