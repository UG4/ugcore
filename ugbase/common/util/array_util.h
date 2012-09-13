//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m09 d21

#ifndef __H__UTIL__ARRAY_UTIL__
#define __H__UTIL__ARRAY_UTIL__

#include <algorithm>

namespace ug
{

///	removes the first occurance of the specified entry.
/**
 * Runs in O(size).
 * copies all entries after the specified one to their predecessor.
 * Be sure that TType supports operator= and operator==
 * \return new size
 */
template <class TType>
int ArrayEraseEntry(TType* array, const TType& entry, size_t size)
{
//	find the entry
	size_t i;
	for(i = 0; i < size; ++i)
	{
		if(array[i] == entry)
			break;
	}

//	proceed if the entry has been found
	if(i >= size)
		return size;

//	copy elements
	for(; i < size - 1; ++i)
		array[i] = array[i+1];

//	done
	return size - 1;
}

///	Swaps the first entry with the given value with the last entry in the list
/**
 * Runs in O(size).
 * Iterates through all entries until the given one is found.
 * Be sure that TType supports operator= and operator==
 */
template <class TType>
void ArraySwapWithLast(TType* array, const TType& entry, size_t size)
{
	using namespace std;
//	find the entry
	size_t i;
	for(i = 0; i < size; ++i){
		if(array[i] == entry)
			break;
	}

//	proceed if the entry has been found and if it is not already the last
	if(i + 1>= size)
		return;

//	swap elements
	swap(array[i], array[size - 1]);
}

///	replaces the first occurance of oldEntry with newEntry
/**
 * Runs in O(size).
 * Be sure that TType supports operator= and operator==
 * \return true if oldEntry was found, false if not.
 */
template <class TType>
int ArrayReplaceEntry(TType* array, const TType& newEntry,
					  const TType& oldEntry, size_t size)
{
//	find the entry
	size_t i;
	for(i = 0; i < size; ++i)
	{
		if(array[i] == oldEntry){
			array[i] = newEntry;
			return true;
		}
	}

	return false;
}

}//	end of namespace

#endif
