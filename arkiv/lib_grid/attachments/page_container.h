/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__page_container__
#define __H__UG__page_container__
/*
#include <vector>
#include <memory>
#include "common/types.h"

namespace ug
{

template <class T, int MAX_PAGE_SIZE = 4096,
		  class Allocator = std::allocator<T> >
class PageContainer
{
	public:
		using allocator_type = Allocator;
		using size_type = typename Allocator::size_type;
		using difference_type = typename Allocator::difference_type;
		using reference = typename Allocator::reference;
		using const_reference = typename Allocator::const_reference;

	public:
		PageContainer();

		PageContainer(const PageContainer& pc);

		~PageContainer();

		PageContainer& operator = (const PageContainer& pc);

		inline size_t size() const;
		inline size_t capacity() const;

		void resize(size_t size, const T& val = T());

		void reserve(size_t size);

		void clear();

		inline T& operator [] (size_t ind);
		inline const T& operator [] (size_t ind) const;

		void swap(PageContainer& pc) noexcept;

	private:
		void assign_container(const PageContainer& pc);

	///	returns the page in which the data for the given index lies
		inline T* get_page(size_t ind) const;

	///	returns the index of the page in which the data for the given index lies
		inline size_t get_page_index(size_t ind) const;

	///	returns the offset that a index has in its page
		inline size_t get_page_offset(size_t ind) const;

	private:
		std::vector<T*>	m_pages;
		const size_t	m_numPageEntries;
		size_t			m_size;
		Allocator		m_alloc;
};

}//	end of namespace

////////////////////////////////
//	include implementation
#include "page_container_impl.h"
*/
#endif
