//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m09 d10

#include "multi_grid.h"

#ifndef __H__LIB_GRID__MULTI_GRID_IMPL__
#define __H__LIB_GRID__MULTI_GRID_IMPL__

namespace ug
{
/*
template <class TChild, class TElem>
int MultiGrid::num_children(TElem* elem)
{

	int sharedPipeSec = elem->shared_pipe_section();
	assert(sharedPipeSec != -1 && "bad shared pipe section!");

	return get_elem_info(elem)->m_children.num_elements(sharedPipeSec);
}

//	child access
template <class TChild, class TElem>
int MultiGrid::get_children(std::vector<TChild*>& vChildrenOut, TElem* elem)
{

	int sharedPipeSec = geometry_traits<TChild>::SHARED_PIPE_SECTION;
	assert(sharedPipeSec != -1 && "bad shared pipe section!");

	MGElementInfo* elemInfo = get_elem_info(elem);

	typename SectionContainer::iterator iter =
				elemInfo->m_children.section_begin(sharedPipeSec);

	typename SectionContainer::iterator iterEnd =
			elemInfo->m_children.section_end(sharedPipeSec);

	int numChildren = elemInfo->m_children.num_elements(sharedPipeSec);

	if(vChildrenOut.capacity() < numChildren)
		vChildrenOut.reserve(numChildren);

	vChildrenOut.clear();

	for(; iter != iterEnd; ++iter)
	{
		vChildrenOut.push_back((TChild*) *iter);
	}

	return vChildrenOut.size();

}
*/

template <class TElem, class TParent>
void MultiGrid::element_created(TElem* elem, TParent* pParent)
{
//	if hierarchical_insertion is enabled, the element will be put
//	into the next higher level of pParents level.

	int level = 0;
	if(pParent)
	{
	//	the element is inserted into a new layer.
		level = get_level(pParent) + 1;
	}

//	register parent and child
	typename mginfo_traits<TElem>::info_type& info = get_info(elem);
	info.m_pParent = pParent;
	if(pParent)
	{
	//	add the element to the parents children list
		typename mginfo_traits<TParent>::info_type& parentInfo = get_info(pParent);
		parentInfo.add_child(elem);
		
	//	set the new status
		switch(parentInfo.m_state)
		{
			case MGES_FIXED:		set_state(elem, MGES_FIXED); break;
			case MGES_CONSTRAINED:	set_state(elem, MGES_CONSTRAINED); break;
			case MGES_CONSTRAINING:	set_state(elem, MGES_CONSTRAINED); break;
		}
	}
	else {
		set_state(elem, MGES_NORMAL);
	}

//	put the element into the hierarchy
	m_hierarchy.assign_subset(elem, level);
}

template <class TElem>
void MultiGrid::element_to_be_erased(TElem* elem)
{
//	we have to remove the elements children as well.
	typename mginfo_traits<TElem>::info_type& info = get_info(elem);
	info.erase_all_children(*this);
}

template <class TElem, class TParent>
void MultiGrid::element_to_be_erased(TElem* elem, TParent* pParent)
{
//	unregister the element from its parent.
	typename mginfo_traits<TParent>::info_type& parentInfo = get_info(pParent);
	parentInfo.remove_child(elem);
	element_to_be_erased(elem);
}

/*
template <class TElem>
void MultiGrid::element_to_be_replaced(TElem* elemOld, TElem* elemNew)
{
}
*/
}//	end of namespace

#endif
