// created by Sebastian Reiter
// y09 m07 d31
// s.b.reiter@googlemail.com

#ifndef __H__LIB_GRID__GEOMETRIC_OBJECT_COLLECTION_IMPL__
#define __H__LIB_GRID__GEOMETRIC_OBJECT_COLLECTION_IMPL__

#include <cassert>
#include "geometric_object_collection.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	GeometricObjectCollection

////////////////////////////////////////////////////////////////////////
//	Iterators
/*
////////////////////////////////////////////////////////////////////////
//	begin
template <class TGeomObj>
typename geometry_traits<TGeomObj>::iterator
GeometricObjectCollection::begin()
{
	assert(geometry_traits<TGeomObj>::BASE_OBJECT_TYPE_ID != -1
		&& "ERROR. Invalid base_object_type of GeomObjType!");

	return iterator_cast<typename geometry_traits<TGeomObj>::iterator>
		(m_pSectionContainers[geometry_traits<TGeomObj>::BASE_OBJECT_TYPE_ID]->section_begin(geometry_traits<TGeomObj>::SHARED_PIPE_SECTION));
}

////////////////////////////////////////////////////////////////////////
//	end
template <class TGeomObj>
typename geometry_traits<TGeomObj>::iterator
GeometricObjectCollection::end()
{
	assert(geometry_traits<TGeomObj>::BASE_OBJECT_TYPE_ID != -1
		&& "ERROR. Invalid base_object_type of GeomObjType!");
		
	return iterator_cast<typename geometry_traits<TGeomObj>::iterator>
		(m_pSectionContainers[geometry_traits<TGeomObj>::BASE_OBJECT_TYPE_ID]->section_end(geometry_traits<TGeomObj>::SHARED_PIPE_SECTION));
}


////////////////////////////////////////////////////////////////////////
//	element numbers
template <class TGeomObj>
uint
GeometricObjectCollection::num()
{
	assert(geometry_traits<TGeomObj>::BASE_OBJECT_TYPE_ID != -1
		&& "ERROR in Grid::num(). Invalid base_object_type of GeomObjType!");

	int objType = geometry_traits<TGeomObj>::BASE_OBJECT_TYPE_ID;
	int secIndex = geometry_traits<TGeomObj>::SHARED_PIPE_SECTION;

	if(secIndex == -1)
		return m_pSectionContainers[objType]->num_elements();

	return m_pSectionContainers[objType]->num_elements(secIndex);
}
*/

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	get_container
template <class TGeomObj> inline
GeometricObjectSectionContainer*
GeometricObjectCollection::get_container(size_t level)
{
	int objType = geometry_traits<TGeomObj>::BASE_OBJECT_TYPE_ID;
	assert(objType >= 0 && objType < 5 && "ERROR. Invalid base_object_type of GeomObjType!");
		
	return m_levels[level].pSectionContainers[objType];
}

////////////////////////////////////////////////////////////////////////
//	begin
template <class TGeomObj>
typename geometry_traits<TGeomObj>::iterator
GeometricObjectCollection::begin(size_t level)
{
	return iterator_cast<typename geometry_traits<TGeomObj>::iterator>
		(get_container<TGeomObj>(level)->section_begin(
			geometry_traits<TGeomObj>::SHARED_PIPE_SECTION));
}

////////////////////////////////////////////////////////////////////////
//	end
template <class TGeomObj>
typename geometry_traits<TGeomObj>::iterator
GeometricObjectCollection::end(size_t level)
{
	return iterator_cast<typename geometry_traits<TGeomObj>::iterator>
		(get_container<TGeomObj>(level)->section_end(
			geometry_traits<TGeomObj>::SHARED_PIPE_SECTION));
}


////////////////////////////////////////////////////////////////////////
//	element numbers
template <class TGeomObj>
size_t
GeometricObjectCollection::num(size_t level)
{
	int secIndex = geometry_traits<TGeomObj>::SHARED_PIPE_SECTION;

	if(secIndex == -1)
		return get_container<TGeomObj>(level)->num_elements();

	return get_container<TGeomObj>(level)->num_elements(secIndex);
}

//	GeometricObjectCollection
template <class TGeomObj>
size_t
GeometricObjectCollection::num()
{
	size_t counter = 0;
	for(size_t i = 0; i < m_levels.size(); ++i)
		counter += num<TGeomObj>(i);
	return counter;
}

}

#endif
