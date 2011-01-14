//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m12 d03

#ifndef __H__LIB_GRID__MISC_UTIL__IMPL__
#define __H__LIB_GRID__MISC_UTIL__IMPL__

#include "misc_util.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
//	CalculateCenter
template <class TIterator, class TAAPosVRT>
typename TAAPosVRT::ValueType
CalculateCenter(TIterator begin, TIterator end, TAAPosVRT& aaPos)
{
	int counter = 0;
	typename TAAPosVRT::ValueType center;
	VecSet(center, 0);
	for(TIterator iter = begin; iter != end; ++iter, ++counter)
		VecAdd(center, center, CalculateCenter(*iter, aaPos));
		
	if(counter > 0)
		VecScale(center, center, 1./(number)counter);
		
	return center;
}

////////////////////////////////////////////////////////////////////////
//	FindByCoordinate
template<class TElem, class TVertexPositionAttachmentAccessor>
TElem* FindByCoordinate(const typename TVertexPositionAttachmentAccessor::ValueType& coord,
						typename geometry_traits<TElem>::iterator iterBegin,
						typename geometry_traits<TElem>::iterator iterEnd,
						TVertexPositionAttachmentAccessor& aaPosVRT)
{
	if(iterBegin == iterEnd)
		return NULL;

	typename geometry_traits<TElem>::iterator iter = iterBegin;
	TElem* bestElem = *iter;
	number bestDistSq = VecDistanceSq(coord, CalculateCenter(bestElem, aaPosVRT));
	iter++;

	while(iter != iterEnd)
	{
		number distSq = VecDistanceSq(coord, CalculateCenter(*iter, aaPosVRT));
		if(distSq < bestDistSq)
		{
			bestDistSq = distSq;
			bestElem = *iter;
		}

		++iter;
	}

	return bestElem;
}

////////////////////////////////////////////////////////////////////////
//	NumSharedVertices
template <class TElemPtr1, class TElemPtr2>
size_t NumSharedVertices(Grid& grid, TElemPtr1 elem1, TElemPtr2 elem2)
{
	grid.begin_marking();
//	first mark all vertices of elem1
	for(size_t i = 0; i < elem1->num_vertices(); ++i)
		grid.mark(elem1->vertex(i));

//	now count how many of vertex 2 are marked.
	size_t counter = 0;
	for(size_t i = 0; i < elem2->num_vertices(); ++i){
		if(grid.is_marked(elem2->vertex(i)))
			++counter;
	}
	
	grid.end_marking();
	
	return counter;
}

////////////////////////////////////////////////////////////////////////
//	EraseElements
template <class TElem>
void EraseElements(Grid& grid, typename geometry_traits<TElem>::iterator iterBegin,
						typename geometry_traits<TElem>::iterator iterEnd)
{
//	be careful to not invalidate the iterators.
	while(iterBegin != iterEnd)
	{
		TElem* e = *iterBegin;
		iterBegin++;
		grid.erase(e);
	}
}

////////////////////////////////////////////////////////////////////////
//	AssignIndices
template <class TElem>
void AssignIndices(typename geometry_traits<TElem>::iterator iterBegin,
					typename geometry_traits<TElem>::iterator iterEnd,
					Grid::AttachmentAccessor<TElem, AInt>& aaInt)
{
	int index = 0;
	while(iterBegin != iterEnd)
	{
		aaInt[*iterBegin] = index++;
		iterBegin++;
	}
}

}//	end of namespace

#endif
