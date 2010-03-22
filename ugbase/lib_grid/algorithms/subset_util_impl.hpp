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

////////////////////////////////////////////////////////////////////////
template <class TElem, class TSubsetHandler>
void AssignAssociatedVerticesToSubsets(TSubsetHandler& sh,
									const ISubsetHandler& srcIndHandler)
{
	typedef typename geometry_traits<TElem>::iterator iterator;
	for(size_t l  = 0; l < sh.num_levels(); ++l){
		for(size_t si = 0; si < sh.num_subsets(); ++si){
			for(iterator iter = sh.template begin<TElem>(si, l);
				iter != sh.template end<TElem>(si, l); ++iter)
			{
				TElem* e = *iter;
				for(size_t i = 0; i < e->num_vertices(); ++i)
				{
					VertexBase* vrt = e->vertex(i);
					sh.assign_subset(vrt, srcIndHandler.get_subset_index(vrt));
				}
			}
		}
	}
}
									
////////////////////////////////////////////////////////////////////////
template <class TElem, class TSubsetHandler>
void AssignAssociatedEdgesToSubsets(TSubsetHandler& sh,
									const ISubsetHandler& srcIndHandler)
{
	typedef typename geometry_traits<TElem>::iterator iterator;
	std::vector<EdgeBase*> vEdges;
	
	for(size_t l  = 0; l < sh.num_levels(); ++l){
		for(size_t si = 0; si < sh.num_subsets(); ++si){
			for(iterator iter = sh.template begin<TElem>(si, l);
				iter != sh.template end<TElem>(si, l); ++iter)
			{
				TElem* e = *iter;
				CollectEdges(vEdges, *sh.get_assigned_grid(), e);
				
				for(size_t i = 0; i < vEdges.size(); ++i)
				{
					EdgeBase* edge = vEdges[i];
					sh.assign_subset(edge, srcIndHandler.get_subset_index(edge));
				}
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////
template <class TElem, class TSubsetHandler>
void AssignAssociatedFacesToSubsets(TSubsetHandler& sh,
									const ISubsetHandler& srcIndHandler)
{
	typedef typename geometry_traits<TElem>::iterator iterator;
	std::vector<Face*> vFaces;
	
	for(size_t l  = 0; l < sh.num_levels(); ++l){
		for(size_t si = 0; si < sh.num_subsets(); ++si){
			for(iterator iter = sh.template begin<TElem>(si, l);
				iter != sh.template end<TElem>(si, l); ++iter)
			{
				TElem* e = *iter;
				CollectFaces(vFaces, *sh.get_assigned_grid(), e);
				
				for(size_t i = 0; i < vFaces.size(); ++i)
				{
					Face* f = vFaces[i];
					sh.assign_subset(f, srcIndHandler.get_subset_index(f));
				}
			}
		}
	}
}

///	helper with with dummy-param for compile-time function selection.
template <class TElem, class TSubsetHandlerDest, class TSubsetHandlerSrc>
void AssignAssociatedLowerDimElemsToSubsets(TSubsetHandlerDest& sh,
									const TSubsetHandlerSrc& srcIndHandler,
									const Volume&)
{
//	we have to find all associated elements of lower dimension.
	if(srcIndHandler.template num<Face>() > 0)
		AssignAssociatedFacesToSubsets<TElem>(sh, srcIndHandler);
	if(srcIndHandler.template num<EdgeBase>() > 0)
		AssignAssociatedEdgesToSubsets<TElem>(sh, srcIndHandler);
	if(srcIndHandler.template num<VertexBase>() > 0)
		AssignAssociatedVerticesToSubsets<TElem>(sh, srcIndHandler);
}

///	helper with with dummy-param for compile-time function selection.
template <class TElem, class TSubsetHandlerDest, class TSubsetHandlerSrc>
void AssignAssociatedLowerDimElemsToSubsets(TSubsetHandlerDest& sh,
									const TSubsetHandlerSrc& srcIndHandler,
									const Face&)
{
//	we have to find all associated elements of lower dimension.
	if(srcIndHandler.template num<EdgeBase>() > 0)
		AssignAssociatedEdgesToSubsets<TElem>(sh, srcIndHandler);
	if(srcIndHandler.template num<VertexBase>() > 0)
		AssignAssociatedVerticesToSubsets<TElem>(sh, srcIndHandler);
}

///	helper with with dummy-param for compile-time function selection.
template <class TElem, class TSubsetHandlerDest, class TSubsetHandlerSrc>
void AssignAssociatedLowerDimElemsToSubsets(TSubsetHandlerDest& sh,
									const TSubsetHandlerSrc& srcIndHandler,
									const EdgeBase&)
{
//	we have to find all associated elements of lower dimension.
	if(srcIndHandler.template num<VertexBase>() > 0)
		AssignAssociatedVerticesToSubsets<TElem>(sh, srcIndHandler);
}

////////////////////////////////////////////////////////////////////////
template <class TElem, class TSubsetHandlerDest, class TSubsetHandlerSrc>
void AssignAssociatedLowerDimElemsToSubsets(TSubsetHandlerDest& sh,
									const TSubsetHandlerSrc& srcIndHandler)
{
	AssignAssociatedLowerDimElemsToSubsets<TElem>(sh,
											srcIndHandler, TElem());
}

}//	end of namespace

#endif
