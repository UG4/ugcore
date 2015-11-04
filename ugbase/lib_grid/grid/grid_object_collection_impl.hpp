#ifndef __H__LIB_GRID__GEOMETRIC_OBJECT_COLLECTION_IMPL__
#define __H__LIB_GRID__GEOMETRIC_OBJECT_COLLECTION_IMPL__

#include <cassert>
#include "grid_object_collection.h"
#include "element_storage.h"

namespace ug
{
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	GridObjectCollection

////////////////////////////////////////////////////////////////////////
//	get_container

template <class TGeomObj> inline
const typename ElementStorage<typename geometry_traits<TGeomObj>::grid_base_object>::
SectionContainer*
GridObjectCollection::get_container(size_t level) const
{
	return SectionContainerSelector<typename geometry_traits<TGeomObj>::grid_base_object>::
			section_container(m_levels[level].vrtContainer, m_levels[level].edgeContainer,
							  m_levels[level].faceContainer, m_levels[level].volContainer);
}

template <class TGeomObj> inline
typename ElementStorage<typename geometry_traits<TGeomObj>::grid_base_object>::
SectionContainer*
GridObjectCollection::get_container(size_t level)
{
	return SectionContainerSelector<typename geometry_traits<TGeomObj>::grid_base_object>::
			section_container(m_levels[level].vrtContainer, m_levels[level].edgeContainer,
							  m_levels[level].faceContainer, m_levels[level].volContainer);
}

////////////////////////////////////////////////////////////////////////
//	begin
template <class TGeomObj>
typename geometry_traits<TGeomObj>::const_iterator
GridObjectCollection::begin(size_t level) const
{
	return iterator_cast<typename geometry_traits<TGeomObj>::const_iterator>
		(get_container<TGeomObj>(level)->section_begin(
			geometry_traits<TGeomObj>::CONTAINER_SECTION));
}

////////////////////////////////////////////////////////////////////////
//	end
template <class TGeomObj>
typename geometry_traits<TGeomObj>::const_iterator
GridObjectCollection::end(size_t level) const
{
	return iterator_cast<typename geometry_traits<TGeomObj>::const_iterator>
		(get_container<TGeomObj>(level)->section_end(
			geometry_traits<TGeomObj>::CONTAINER_SECTION));
}

////////////////////////////////////////////////////////////////////////
//	begin
template <class TGeomObj>
typename geometry_traits<TGeomObj>::iterator
GridObjectCollection::begin(size_t level)
{
	return iterator_cast<typename geometry_traits<TGeomObj>::iterator>
		(get_container<TGeomObj>(level)->section_begin(
			geometry_traits<TGeomObj>::CONTAINER_SECTION));
}

////////////////////////////////////////////////////////////////////////
//	end
template <class TGeomObj>
typename geometry_traits<TGeomObj>::iterator
GridObjectCollection::end(size_t level)
{
	return iterator_cast<typename geometry_traits<TGeomObj>::iterator>
		(get_container<TGeomObj>(level)->section_end(
			geometry_traits<TGeomObj>::CONTAINER_SECTION));
}

////////////////////////////////////////////////////////////////////////
//	element numbers
template <class TGeomObj>
size_t
GridObjectCollection::num(size_t level) const
{
	int secIndex = geometry_traits<TGeomObj>::CONTAINER_SECTION;

	if(secIndex == -1)
		return get_container<TGeomObj>(level)->num_elements();

	return get_container<TGeomObj>(level)->num_elements(secIndex);
}

//	GridObjectCollection
template <class TGeomObj>
size_t
GridObjectCollection::num() const
{
	size_t counter = 0;
	for(size_t i = 0; i < m_levels.size(); ++i)
		counter += num<TGeomObj>(i);
	return counter;
}

}

#endif
