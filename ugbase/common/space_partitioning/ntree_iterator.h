// created by Sebastian Reiter
// s.b.reiter@gmail.com
// Sep 5, 2013

#ifndef __H__UG__ntree_iterator__
#define __H__UG__ntree_iterator__

#include <cassert>

namespace ug{

///	this iterator is used by the ntree class to provide access to the elements of a given node
template <class elem_t, class entry_t>
class const_ntree_element_iterator
{
	public:
		typedef const_ntree_element_iterator	this_type;
		typedef std::forward_iterator_tag		iterator_category;
		typedef size_t							difference_type;
		typedef elem_t*							pointer;
		typedef elem_t							value_type;
		typedef value_type&						reference;

		const_ntree_element_iterator() : m_entries(NULL), m_entryInd(s_invalidIndex)	{}
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
		static const size_t s_invalidIndex = -1;

		const entry_t*	m_entries;
		size_t			m_entryInd;
};

}// end of namespace

#endif
