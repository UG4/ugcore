//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m09 d21

#ifndef __H__UTIL__ARRAY_UTIL__
#define __H__UTIL__ARRAY_UTIL__

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
int ArrayEraseEntry(TType* array, const TType& entry, uint size)
{
//	find the entry
	uint i;
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

}//	end of namespace

#endif
