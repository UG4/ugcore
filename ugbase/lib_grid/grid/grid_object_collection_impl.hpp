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
