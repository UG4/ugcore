#ifndef __H__UG__hash_impl__
#define __H__UG__hash_impl__

#include <cassert>
#include <algorithm>
#include <limits>
#include "common/error.h"
#include "hash.h"

namespace ug{

template <class TKey, class TValue>
Hash<TKey, TValue>::
Hash() :
	m_numEntries(0),
	m_firstUnusedEntry(invalid_index())
{
	resize_hash(1);
}


template <class TKey, class TValue>
Hash<TKey, TValue>::
Hash(size_t hashSize) :
	m_numEntries(0),
	m_firstUnusedEntry(invalid_index())
{
	resize_hash(hashSize);
}


template <class TKey, class TValue>
void Hash<TKey, TValue>::
resize_hash(size_t size)
{
//todo: This method could be more efficient by reusing memory in m_hashList...

//	create a new hash and insert all existing items into it
//	when creating the hash-list, we'll reserve some additional memory to
//	avoid frequent reallocations.
//	note: during the first resize no additional memory will be allocated.
	using namespace std;
	vector<pair<size_t, size_t> > newHashList;
	newHashList.resize(max<size_t>(1, size),
					   pair<size_t, size_t>(invalid_index(), invalid_index()));

	if(m_numEntries > 0){
		for(size_t i = 0; i < m_hashList.size(); ++i){
			size_t curInd = m_hashList[i].first;
			while(curInd != invalid_index()){
				Entry& e = m_entries[curInd];
				size_t hi = hash_key(e.key) % size;

				if(newHashList[hi].first == invalid_index())
					newHashList[hi].first = newHashList[hi].second = curInd;
				else{
					m_entries[newHashList[hi].second].next = curInd;
					newHashList[hi].second = curInd;
				}

				curInd = e.next;
				e.next = invalid_index();
			}
		}
	}

	m_hashList.swap(newHashList);
}


template <class TKey, class TValue>
size_t Hash<TKey, TValue>::
hash_size() const
{
	return m_hashList.size();
}


template <class TKey, class TValue>
void Hash<TKey, TValue>::
reserve(size_t size)
{
	m_entries.reserve(size);
}


template <class TKey, class TValue>
size_t Hash<TKey, TValue>::
capacity() const
{
	return m_entries.capacity();
}


template <class TKey, class TValue>
void Hash<TKey, TValue>::
clear()
{
	using namespace std;
	m_entries.clear();
	m_numEntries = 0;
	m_firstUnusedEntry = invalid_index();
	m_hashList.assign(m_hashList.size(),
					  make_pair<size_t, size_t>(invalid_index(), invalid_index()));
}


template <class TKey, class TValue>
bool Hash<TKey, TValue>::
empty() const
{
	return m_numEntries == 0;
}


template <class TKey, class TValue>
bool Hash<TKey, TValue>::
has_entry(const key_t& key) const
{
	return find_entry(key) != invalid_index();
}


template <class TKey, class TValue>
TValue& Hash<TKey, TValue>::
get_entry(const key_t& key)
{
	size_t eind = find_entry(key);
	assert((eind != invalid_index()) && "No such entry. Check existance with has_entry first!");

	if(eind == invalid_index()){
		UG_THROW("No entry exists for the specified key. Please call 'has_entry' first"
				" or use the alternate version of 'get_entry', which returns a bool.");
	}
	return m_entries[eind].value;
}


template <class TKey, class TValue>
const TValue& Hash<TKey, TValue>::
get_entry(const key_t& key) const
{
	size_t eind = find_entry(key);
	assert((eind != invalid_index()) && "No such entry. Check existance with has_entry first!");

	if(eind == invalid_index()){
		UG_THROW("No entry exists for the specified key. Please call 'has_entry' first"
				" or use the alternate version of 'get_entry', which returns a bool.");
	}
	return m_entries[eind].value;
}

template <class TKey, class TValue>
bool Hash<TKey, TValue>::
get_entry(TValue& valOut, const key_t& key) const
{
	size_t eind = find_entry(key);

	if(eind == invalid_index())
		return false;

	valOut = m_entries[eind].value;
	return true;
}


template <class TKey, class TValue>
void Hash<TKey, TValue>::
insert(const key_t& key, const value_t& val)
{
	size_t eind = m_firstUnusedEntry;
	if(eind == invalid_index()){
		eind = m_entries.size();
		m_entries.push_back(Entry(key, val));
	}
	else{
		m_firstUnusedEntry = m_entries[eind].next;
		m_entries[eind] = Entry(key, val);
	}

	size_t hi = hash_index(key);

	if(m_hashList[hi].first == invalid_index()){
		m_hashList[hi].first = m_hashList[hi].second = eind;
	}
	else{
		m_entries[m_hashList[hi].second].next = eind;
		m_hashList[hi].second = eind;
	}

	++m_numEntries;
}


template <class TKey, class TValue>
void Hash<TKey, TValue>::
erase(const key_t& key)
{
	size_t hi = hash_index(key);
	size_t prevInd = invalid_index();
	size_t curInd = m_hashList[hi].first;
	while(curInd != invalid_index()){
		if(m_entries[curInd].key == key)
			break;
		prevInd = curInd;
		curInd = m_entries[curInd].next;
	}

	if(curInd != invalid_index()){
	//	if curInd was the first entry for this hash-index, we have to adjust
	//	the beginning of the sequence. If not, we have to adjust the next pointer
	//	of the previous entry.
		if(prevInd == invalid_index())
			m_hashList[hi].first = m_entries[curInd].next;
		else
			m_entries[prevInd].next = m_entries[curInd].next;

	//	if curInd was the last entry for this hash-index, we have to adjust
	//	the end of the sequence.
		if(m_entries[curInd].next == invalid_index())
			m_hashList[hi].second = prevInd;
	}

//	curInd has to be added to the list of unused entries:
	m_entries[curInd].next = m_firstUnusedEntry;
	m_firstUnusedEntry = curInd;
	--m_numEntries;
}


template <class TKey, class TValue>
typename Hash<TKey, TValue>::iterator Hash<TKey, TValue>::
begin(const key_t& key)
{
	if(empty())
		return end(key);
	return iterator(key, &m_entries.front(), m_hashList[hash_index(key)].first);
}


template <class TKey, class TValue>
typename Hash<TKey, TValue>::iterator Hash<TKey, TValue>::
end(const key_t& key)
{
	return iterator(key, NULL, invalid_index());
}

template <class TKey, class TValue>
size_t Hash<TKey, TValue>::
hash_index(const key_t& key) const
{
	return hash_key(key) % m_hashList.size();
}

template <class TKey, class TValue>
size_t Hash<TKey, TValue>::
find_entry(const key_t& key) const
{
	size_t hi = hash_index(key);
	size_t curInd = m_hashList[hi].first;
	while(curInd != invalid_index()){
		if(m_entries[curInd].key == key)
			return curInd;
		curInd = m_entries[curInd].next;
	}
	return invalid_index();
}

}// end of namespace

#endif
