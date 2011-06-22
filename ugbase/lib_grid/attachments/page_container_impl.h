// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 22.06.2011 (m,d,y)

#ifndef __H__UG__page_container_impl__
#define __H__UG__page_container_impl__

#include <algorithm>
#include "page_container.h"

namespace ug
{

template <class T, int MAX_PAGE_SIZE, class Allocator>
PageContainer::
PageContainer() :
	m_numPageEntries((size_t)(MAX_PAGE_SIZE / sizeof(T)))
{
}

template <class T, int MAX_PAGE_SIZE, class Allocator>
PageContainer::
PageContainer(const PageContainer& pc)
{
	assign_container(pc);
}

template <class T, int MAX_PAGE_SIZE, class Allocator>
PageContainer::
~PageContainer()
{
	clear();
	for(size_t i = 0; i < m_pages.size(); ++i){
		delete[] m_pages[i];
	}
}

template <class T, int MAX_PAGE_SIZE, class Allocator>
PageContainer& PageContainer::
operator=(const PageContainer& pc)
{
	assign_container(pc);
	return *this;
}

template <class T, int MAX_PAGE_SIZE, class Allocator>
size_t PageContainer::
size() const
{
	return m_size;
}

template <class T, int MAX_PAGE_SIZE, class Allocator>
size_t PageContainer::
capacity() const
{
	return m_capacity;
}

template <class T, int MAX_PAGE_SIZE, class Allocator>
void PageContainer::
resize(size_t size)
{
	using namespace std;

//	allocate memory
	reserve(size);

//	call constructor on new objects
//	to do this with optimal performance,
//	we'll iterate over the pages directly
	while(m_size < size){
		T* page = get_page(m_size);
		size_t offset = get_page_offset(m_size);
		const size_t maxI = min(m_numPageEntries, offset + size - m_size);

		for(size_t i = offset; i < maxI; ++i)
			m_alloc.construct(page[i], val);

		m_size += maxI;
	}

//	if resize shrinks the data array, we have to call the destructors of
//	deleted objects. At this point size <= m_size.
	while(m_size > size){
		T* page = get_page(m_size);
		size_t maxI = get_page_offset(m_size) + 1;
		const size_t minI = 0;
		const size_t diff = m_size - size;
		if(offset > diff)
			minI = offset - diff + 1;

		for(size_t i = minI; i < maxI; ++i)
			m_alloc.destroy(page[i]);

		m_size -= (maxI - minI);
	}
}

template <class T, int MAX_PAGE_SIZE, class Allocator>
void PageContainer::
reserve(size_t size)
{
	while(m_pages.size() * m_numPageEntries < size){
		T* buf = m_alloc.allocate(m_numPageEntries);
		m_pages.push_back(buf);
	//todo:....
		size_t offset = get_page_offset(m_size);
		const size_t maxI = min(m_numPageEntries, offset + size - m_size);

		for(size_t i = offset; i < maxI; ++i)
			m_alloc.construct(page[i], val);

		m_size += maxI;
	}
}

template <class T, int MAX_PAGE_SIZE, class Allocator>
void PageContainer::
clear()
{

}

template <class T, int MAX_PAGE_SIZE>
T& PageContainer::
operator[](size_t ind)
{

}

template <class T, int MAX_PAGE_SIZE>
const T& PageContainer::
operator[](size_t ind) const
{

}

template <class T, int MAX_PAGE_SIZE>
void PageContainer::
swap(PageContainer& pc)
{

}

template <class T, int MAX_PAGE_SIZE>
void PageContainer::
assign_container(const PageContainer& pc)
{

}

template <class T, int MAX_PAGE_SIZE>
inline T* PageContainer::
get_page(size_t ind)
{
	assert(get_page_index(ind) < m_pages.size());
	return m_pages[get_page_index(ind)];
}

template <class T, int MAX_PAGE_SIZE>
inline size_t PageContainer::
get_page_index(size_t ind)
{
	return ind / m_numPageEntries;
}

template <class T, int MAX_PAGE_SIZE>
inline size_t PageContainer::
get_page_offset(size_t ind)
{
	return ind % numPageEntries;
}

}//	end of namespace

#endif
