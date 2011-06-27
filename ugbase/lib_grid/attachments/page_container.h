// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 22.06.2011 (m,d,y)

#ifndef __H__UG__page_container__
#define __H__UG__page_container__

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
		typedef Allocator allocator_type;
		typedef typename Allocator::size_type size_type;
		typedef typename Allocator::difference_type difference_type;
		typedef typename Allocator::reference reference;
		typedef typename Allocator::const_reference const_reference;

	public:
		PageContainer();

		PageContainer(const PageContainer& pc);

		~PageContainer();

		PageContainer& operator=(const PageContainer& pc);

		inline size_t size() const;
		inline size_t capacity() const;

		void resize(size_t size, const T& val = T());

		void reserve(size_t size);

		void clear();

		inline T& operator[](size_t ind);
		inline const T& operator[](size_t ind) const;

		void swap(PageContainer& pc);

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

#endif
