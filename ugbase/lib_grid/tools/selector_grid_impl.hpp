/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
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

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	...
////////////////////////////////////////////////////////////////////////

#ifndef __H__LIBGRID__SELECTOR_GRID_IMPL__
#define __H__LIBGRID__SELECTOR_GRID_IMPL__

#include <cassert>

namespace ug
{

template <typename TElem>
inline int
Selector::get_section_index() const
{
	return geometry_traits<TElem>::CONTAINER_SECTION;
}

template <typename TElem>
inline void
Selector::clear()
{
	if(m_pGrid){
	//	mark all elements as deselected
		typename geometry_traits<TElem>::iterator iter;
		for(iter = begin<TElem>(); iter != end<TElem>(); ++iter)
			mark_deselected(*iter);

	//	clear the section
		const int sInd = get_section_index<TElem>();
		if(sInd < 0)
			section_container<TElem>().clear();
		else
			section_container<TElem>().clear_section(sInd);
	}
}

template <typename TElem>
inline size_t
Selector::num() const
{
	const int sInd = get_section_index<TElem>();
	if(sInd < 0)
		return section_container<TElem>().num_elements();
	else
		return section_container<TElem>().num_elements(sInd);
}

inline size_t
Selector::num() const
{
	return num<Vertex>() + num<Edge>() + num<Face>() + num<Volume>();
}

//	empty
inline bool 
Selector::empty() const
{
	return num() == 0;
}

template <typename TElem>
inline bool 
Selector::empty() const
{
	return num<TElem>() == 0;
}

//	begin
template <typename TElem>
inline typename geometry_traits<TElem>::iterator
Selector::begin()
{
	const int sInd = get_section_index<TElem>();
	if(sInd < 0)
		return iterator_cast<typename geometry_traits<TElem>::iterator>(
								section_container<TElem>().begin());
	else
		return iterator_cast<typename geometry_traits<TElem>::iterator>(
					section_container<TElem>().section_begin(sInd));
}

//	const begin
template <typename TElem>
inline typename geometry_traits<TElem>::const_iterator
Selector::begin() const
{
	const int sInd = get_section_index<TElem>();
	if(sInd < 0)
		return iterator_cast<typename geometry_traits<TElem>::const_iterator>(
								section_container<TElem>().begin());
	else
		return iterator_cast<typename geometry_traits<TElem>::const_iterator>(
					section_container<TElem>().section_begin(sInd));
}

//	end
template <typename TElem>
inline typename geometry_traits<TElem>::iterator
Selector::end()
{
	const int sInd = get_section_index<TElem>();
	if(sInd < 0)
		return iterator_cast<typename geometry_traits<TElem>::iterator>(
									section_container<TElem>().end());
	else
		return iterator_cast<typename geometry_traits<TElem>::iterator>(
								section_container<TElem>().section_end(sInd));
}

//	const end
template <typename TElem>
inline typename geometry_traits<TElem>::const_iterator
Selector::end() const
{
	const int sInd = get_section_index<TElem>();
	if(sInd < 0)
		return iterator_cast<typename geometry_traits<TElem>::const_iterator>(
									section_container<TElem>().end());
	else
		return iterator_cast<typename geometry_traits<TElem>::const_iterator>(
								section_container<TElem>().section_end(sInd));
}

template <typename TElem>
TElem*
Selector::front()
{
	const int sInd = get_section_index<TElem>();
	return static_cast<TElem*>(section_container<TElem>().front(sInd));
}

template <typename TElem>
TElem*
Selector::back()
{
	const int sInd = get_section_index<TElem>();
	return static_cast<TElem*>(section_container<TElem>().back(sInd));
}

////////////////////////////////////////
//	for compatibility with MGSelector


inline size_t Selector::
num_levels() const
{
	return 1;
}

inline uint Selector::
num(size_t) const
{
	return (uint)num();
}

template <typename TElem>
inline size_t Selector::
num(size_t) const
{
	return num<TElem>();
}

inline bool Selector::
empty(size_t) const
{
	return empty();
}

template <typename TElem>
inline bool Selector::
empty(size_t) const
{
	return empty<TElem>();
}

template <typename TElem>
inline typename geometry_traits<TElem>::iterator
Selector::begin(size_t)
{
	return begin<TElem>();
}

//	end
///	calls end<TElem>();
template <typename TElem>
inline typename geometry_traits<TElem>::iterator
Selector::end(size_t)
{
	return end<TElem>();
}

template <typename TElem>
typename Grid::traits<TElem>::SectionContainer&
Selector::
section_container()
{
	return SectionContainerSelector<typename geometry_traits<TElem>::grid_base_object>::
			section_container(m_vertices, m_edges, m_faces, m_volumes);
}


template <typename TElem>
const typename Grid::traits<TElem>::SectionContainer&
Selector::
section_container() const
{
	return SectionContainerSelector<typename geometry_traits<TElem>::grid_base_object>::
			section_container(m_vertices, m_edges, m_faces, m_volumes);
}

}//	end of namespace

#endif
