//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m07 d21

#ifndef __H__LIB_GRID__SUBSET_UTIL_IMPL__
#define __H__LIB_GRID__SUBSET_UTIL_IMPL__

#include "subset_util.h"

namespace ug
{
////////////////////////////////////////////////////////////////////////
//	FindFirstFreeSubset
template <class TElem>
int GetMaxSubsetIndex(SubsetHandler& sh)
{
//	go from back to front
	for(int i = (int)sh.num_subsets() - 1; i >= 0; --i)
	{
		if(sh.num_elements<TElem>(i) > 0)
		{
		//	this is the highest subset that contains elements of type TElem
			return i;
		}
	}

//	no subset contains elements of type TElem.
	return -1;
}

////////////////////////////////////////////////////////////////////////
//	MakeSubsetsConsecutive
template <class TElem>
void MakeSubsetsConsecutive(SubsetHandler& sh)
{
//	TODO: this algo could be slightly improved regarding runtime.

//	iterate through all subsets.
	for(uint i = 0; i < sh.num_subsets(); ++i)
	{
	//	check whether the subset is empty
		if(sh.num_elements<TElem>(i) == 0)
		{
		//	it is. find the next filled one.
			for(uint j = i + 1; j < sh.num_subsets(); ++j)
			{
				if(sh.num_elements<TElem>(j) > 0)
				{
				//	this subset is filled. move it to position i.
					sh.move_subset(j, i);
					break;
				}
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////
//	AssignAssociatedVerticesToSubset
template <class TIterator>
void AssignAssociatedVerticesToSubset(ISubsetHandler& sh, TIterator elemsBegin,
										TIterator elemsEnd, int subsetIndex)
{
//	iterate through the elements
	for(;elemsBegin != elemsEnd; elemsBegin++)
	{
		typename TIterator::value_type elem = *elemsBegin;
		uint numVrts = elem->num_vertices();
	//	iterate through the vertices of elem and assign them
		for(uint i = 0; i < numVrts; ++i)
			sh.assign_subset(elem->vertex(i), subsetIndex);
	}
}

}//	end of namespace

#endif
