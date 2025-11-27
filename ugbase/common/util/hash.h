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
template <typename TKey, typename TValue>
class Hash
{
	private:
		struct Entry;

	public:
		using key_t = TKey;
		using value_t = TValue;
		using iterator = hash_iterator<key_t, value_t, Entry>;
		// using const_iterator = const_hash_iterator<key_t, value_t, Entry>;

		Hash();
		explicit Hash(size_t hashSize);

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

		static inline constexpr size_t invalid_index() {return -1;}

		struct Entry{
			key_t	key;
			value_t	value;
			size_t	next;

			Entry()	= default;
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
