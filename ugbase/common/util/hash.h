// created by Sebastian Reiter
// s.b.reiter@gmail.com
// oct 2013

#ifndef __H__UG__hash__
#define __H__UG__hash__

#include <vector>
#include <utility>
#include "hash_iterator.h"
#include "hash_function.h"

namespace ug{

///	An associative container for key-value pairs, which provides fast access using hash-keys
/**	\addtogroup ugbase_common_util
 */
template <class TKey, class TValue>
class Hash
{
	private:
		struct Entry;

	public:
		typedef TKey	key_t;
		typedef TValue	value_t;
		typedef hash_iterator<key_t, value_t, Entry>		iterator;
//		typedef const_hash_iterator<key_t, value_t, Entry>	const_iterator;

		Hash();
		Hash(size_t hashSize);

		void resize_hash(size_t size);
		size_t hash_size() const;

//		void enable_auto_hash_resize(float maxElemToHashSizeRatio,
//									 float bestElemToHashSizeRatio);
//
//		void disable_auto_hash_resize();


	///	Reserves memory for key-value-pair storage.
	/** \sa capacity*/
		void reserve(size_t size);

	///	returns the capacity of the container in which key-value-pairs are stored.
	/**	Use reserve to adjust this capacity. If the capacity would be too small
	 * to insert a new element, it will be automatically increased on element insertion.
	 * \sa reserve*/
		size_t capacity() const;

	///	returns the number of key-value-pairs currently stored in the hash
		size_t size() const;

		void clear();

		bool empty() const;
		bool has_entry(const key_t& key) const;

		value_t& get_entry(const key_t& key);
		const value_t& get_entry(const key_t& key) const;
		bool get_entry(value_t& valOut, const key_t& key) const;

		void insert(const key_t& key, const value_t& val);
		void erase(const key_t& key);

		iterator begin(const key_t& key);
		iterator end(const key_t& key);

//		const_iterator begin(const key_t& key) const;
//		const_iterator end(const key_t& key) const;

	private:
		size_t hash_index(const key_t& key) const;
		size_t find_entry(const key_t& key) const;
		inline size_t invalid_index() const	{return -1;}

		struct Entry{
			key_t	key;
			value_t	value;
			size_t	next;

			Entry()	{}
			Entry(const key_t& k, const value_t& v) :
				key(k), value(v), next(-1)
			{}
		};

		std::vector<Entry>	m_entries;
		size_t				m_numEntries;
		size_t				m_firstUnusedEntry;

	/**	each entry holds a pair of indices pointing to the first and the
	 * last entry for the given hash-key.*/
		std::vector<std::pair<size_t, size_t> >	m_hashList;

};

}// end of namespace


#include "hash_impl.hpp"

#endif
