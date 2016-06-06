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

#ifndef __H__UTIL__ARRAY_UTIL__
#define __H__UTIL__ARRAY_UTIL__

#include <algorithm>

namespace ug
{

/// \addtogroup ugbase_common_types
/// \{

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
	return (int)size - 1;
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
bool ArrayReplaceEntry(TType* array, const TType& newEntry,
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

// end group ugbase_common_types
/// \}

}//	end of namespace

#endif
