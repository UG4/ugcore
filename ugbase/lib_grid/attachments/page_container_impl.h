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

#ifndef __H__UG__page_container_impl__
#define __H__UG__page_container_impl__

#include <algorithm>
#include "page_container.h"

//	a temporary include
#include "common/log.h"

namespace ug
{

template <class T, int MAX_PAGE_SIZE, class Allocator>
PageContainer<T, MAX_PAGE_SIZE, Allocator>::
PageContainer() :
	m_numPageEntries((size_t)(MAX_PAGE_SIZE / sizeof(T))),
	m_size(0)
{
}

template <class T, int MAX_PAGE_SIZE, class Allocator>
PageContainer<T, MAX_PAGE_SIZE, Allocator>::
PageContainer(const PageContainer& pc) :
	m_numPageEntries((size_t)(MAX_PAGE_SIZE / sizeof(T)))
{
	assign_container(pc);
}

template <class T, int MAX_PAGE_SIZE, class Allocator>
PageContainer<T, MAX_PAGE_SIZE, Allocator>::
~PageContainer()
{
	clear();
	for(size_t i = 0; i < m_pages.size(); ++i){
		m_alloc.deallocate(m_pages[i], m_numPageEntries);
	}
}

template <class T, int MAX_PAGE_SIZE, class Allocator>
PageContainer<T, MAX_PAGE_SIZE, Allocator>&
PageContainer<T, MAX_PAGE_SIZE, Allocator>::
operator=(const PageContainer& pc)
{
	assign_container(pc);
	return *this;
}

template <class T, int MAX_PAGE_SIZE, class Allocator>
size_t PageContainer<T, MAX_PAGE_SIZE, Allocator>::
size() const
{
	return m_size;
}

template <class T, int MAX_PAGE_SIZE, class Allocator>
size_t PageContainer<T, MAX_PAGE_SIZE, Allocator>::
capacity() const
{
	return m_pages.size() * m_numPageEntries;
}

template <class T, int MAX_PAGE_SIZE, class Allocator>
void PageContainer<T, MAX_PAGE_SIZE, Allocator>::
resize(size_t size, const T& val)
{
	using namespace std;
//UG_LOG("    reserving...\n");
//	allocate memory
	reserve(size);
//UG_LOG("    constructing...\n");
//	call constructor on new objects
//	to do this with optimal performance,
//	we'll iterate over the pages directly
	while(m_size < size){
		T* page = get_page(m_size);
		size_t offset = get_page_offset(m_size);
		const size_t maxI = min(m_numPageEntries, offset + size - m_size);

		for(size_t i = offset; i < maxI; ++i)
			m_alloc.construct(page + i, val);

		m_size += maxI - offset;
	}
//UG_LOG("    destroying...\n");
//	if resize shrinks the data array, we have to call the destructors of
//	deleted objects. At this point size <= m_size.
	while(m_size > size){
		T* page = get_page(m_size - 1);
		size_t maxI = get_page_offset(m_size - 1) + 1;
		size_t minI = 0;
		const size_t diff = m_size - size;
		if(maxI > diff)
			minI = (maxI) - diff;

		for(size_t i = minI; i < maxI; ++i)
			m_alloc.destroy(page + i);

		m_size -= (maxI - minI);
	}
//UG_LOG("    done...\n");
}

template <class T, int MAX_PAGE_SIZE, class Allocator>
void PageContainer<T, MAX_PAGE_SIZE, Allocator>::
reserve(size_t size)
{
	using namespace std;

	//UG_LOG("*num page entries: " << m_numPageEntries << ", required size: " << size << endl);
	while(m_pages.size() * m_numPageEntries < size){
		//UG_LOG("**allocating " << m_numPageEntries << " elems of size " << sizeof(T) << endl);
		T* buf = m_alloc.allocate(m_numPageEntries);
		//UG_LOG("**adding a new page (current num pages: " << m_pages.size() << ")\n");
		m_pages.push_back(buf);
	}
	//UG_LOG("*done\n");
}

template <class T, int MAX_PAGE_SIZE, class Allocator>
void PageContainer<T, MAX_PAGE_SIZE, Allocator>::
clear()
{
	resize(0);
}

template <class T, int MAX_PAGE_SIZE, class Allocator>
T& PageContainer<T, MAX_PAGE_SIZE, Allocator>::
operator[](size_t ind)
{
	assert(ind < m_size);
	return get_page(ind)[get_page_offset(ind)];
}

template <class T, int MAX_PAGE_SIZE, class Allocator>
const T& PageContainer<T, MAX_PAGE_SIZE, Allocator>::
operator[](size_t ind) const
{
	assert(ind < m_size);
	return get_page(ind)[get_page_offset(ind)];
}

template <class T, int MAX_PAGE_SIZE, class Allocator>
void PageContainer<T, MAX_PAGE_SIZE, Allocator>::
swap(PageContainer& pc)
{
	m_pages.swap(pc.m_pages);
	Allocator talloc = m_alloc;
	m_alloc = pc.m_alloc;
	pc.m_alloc = talloc;

	size_t tmp = m_size;
	m_size = pc.m_size;
	pc.m_size = tmp;
}

template <class T, int MAX_PAGE_SIZE, class Allocator>
void PageContainer<T, MAX_PAGE_SIZE, Allocator>::
assign_container(const PageContainer& pc)
{
	using namespace std;
	clear();

	reserve(pc.m_size);

//	copy all entries with optimal performance
	while(m_size < pc.m_size){
		T* page = get_page(m_size);
		T* srcPage = pc.get_page(m_size);
		size_t offset = get_page_offset(m_size);
		const size_t maxI = min(m_numPageEntries, offset + pc.m_size - m_size);

		for(size_t i = offset; i < maxI; ++i)
			m_alloc.construct(page + i, srcPage[i]);

		m_size += maxI;
	}
}

template <class T, int MAX_PAGE_SIZE, class Allocator>
inline T* PageContainer<T, MAX_PAGE_SIZE, Allocator>::
get_page(size_t ind) const
{
	assert(get_page_index(ind) < m_pages.size());
	return m_pages[get_page_index(ind)];
}

template <class T, int MAX_PAGE_SIZE, class Allocator>
inline size_t PageContainer<T, MAX_PAGE_SIZE, Allocator>::
get_page_index(size_t ind) const
{
	return ind / m_numPageEntries;
}

template <class T, int MAX_PAGE_SIZE, class Allocator>
inline size_t PageContainer<T, MAX_PAGE_SIZE, Allocator>::
get_page_offset(size_t ind) const
{
	return ind % m_numPageEntries;
}

}//	end of namespace

#endif
