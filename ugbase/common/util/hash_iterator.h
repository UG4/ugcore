#ifndef __H__UG__hash_iterator__
#define __H__UG__hash_iterator__

#include <cassert>
namespace ug{

///	this iterator is used by the hash class to provide access to the elements of a given key
template <class TKey, class TValue, class TEntry>
class hash_iterator
{
	public:
		typedef hash_iterator					this_type;
		typedef std::forward_iterator_tag		iterator_category;
		typedef size_t							difference_type;
		typedef TValue*							pointer;
		typedef TValue							value_type;
		typedef TValue&							reference;

		typedef TKey 	key_t;
		typedef TValue 	value_t;
		typedef TEntry	entry_t;


		hash_iterator() : m_entries(NULL), m_entryInd(s_invalidIndex)	{}
		hash_iterator(const key_t& key, const entry_t* entries, size_t entryInd) :
			m_key(key), m_entries(entries), m_entryInd(entryInd)	{}

		this_type operator ++()							{increment(); return *this;}
		this_type operator ++(int unused)				{this_type i = *this; increment(); return i;}

		bool operator ==(const this_type& iter) const	{return equal(iter);}
		bool operator !=(const this_type& iter) const	{return !equal(iter);}

		value_type& operator *()
		{
			assert(m_entryInd != s_invalidIndex);
			return m_entries[m_entryInd].value;
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
			m_entryInd = m_entries[m_entryInd].next;
			while(m_entryInd != s_invalidIndex){
				if(m_entries[m_entryInd].key == m_key)
					break;
				m_entryInd = m_entries[m_entryInd].next;
			}
		}

	///	marks an index as invalid
		static const size_t s_invalidIndex = -1;

		key_t			m_key;
		const entry_t*	m_entries;
		size_t			m_entryInd;
};

}// end of namespace

#endif
