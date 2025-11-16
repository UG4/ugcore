/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__ntree_iterator__
#define __H__UG__ntree_iterator__

#include <cassert>

namespace ug{

///	this iterator is used by the ntree class to provide access to the elements of a given node
template <class elem_t, class entry_t>
class const_ntree_element_iterator
{
	public:
		using this_type = const_ntree_element_iterator;
		using iterator_category = std::forward_iterator_tag;
		using difference_type = size_t;
		using pointer = elem_t*;
		using value_type = elem_t;
		using reference = value_type&;

		const_ntree_element_iterator() : m_entries(nullptr), m_entryInd(s_invalidIndex)	{}
		const_ntree_element_iterator(const entry_t* entries, size_t entryInd) :
			m_entries(entries), m_entryInd(entryInd)	{}

		this_type operator ++()							{increment(); return *this;}
		this_type operator ++(int unused)				{this_type i = *this; increment(); return i;}

		bool operator ==(const this_type& iter) const	{return equal(iter);}
		bool operator !=(const this_type& iter) const	{return !equal(iter);}

		value_type operator *()	const
		{
			assert(m_entryInd != s_invalidIndex);
			return m_entries[m_entryInd].elem;
		}

	private:
		inline bool equal(const this_type& other) const
		{
			return m_entryInd == other.m_entryInd;
		}

		void increment()
		{
			assert(m_entries);
			assert(m_entryInd != s_invalidIndex);
			m_entryInd = m_entries[m_entryInd].nextEntryInd;
		}

	///	marks an index as invalid
		static constexpr size_t s_invalidIndex = -1;

		const entry_t*	m_entries;
		size_t			m_entryInd;
};

}// end of namespace

#endif
