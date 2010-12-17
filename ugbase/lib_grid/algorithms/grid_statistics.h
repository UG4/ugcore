//	created by Martin Stepniewski, Sebastian Reiter
//	mastep@gmx.de
//	y09 m11 d11

#ifndef __H__UG__GRID_STATISTICS__
#define __H__UG__GRID_STATISTICS__

#include <vector>
#include "lib_grid/lg_base.h"

namespace ug
{
	
//**********************************************************************
//								declarations
//**********************************************************************

////////////////////////////////////////////////////////////////////////
//	AssignTetrahedronAttributesByAspectRatio - mstepnie
/// assigns tetrahedral elements of a grid to subsets respecting their aspect ratio
bool AssignTetrahedronAttributesByAspectRatio(Grid& grid,
											  SubsetHandler& shVolume,
											  AInt& aTetrahedronAspectRatioClass,
											  std::vector<double>& offsets);

////////////////////////////////////////////////////////////////////////
///	assigns a subset based on the quality of the given element.
/**
 * Currently only faces are supported.
 *
 * \param intervals contains the intervals which define into which subset
 *					an element goes. Numbers have to be sorted, starting at
 *					0 and ending at 1 (0 and 1 should be contained in intervals).
 */
/*
template <class TIterator>
bool AssignSubsetByQuality(Grid& grid, SubsetHandler& sh,
						   TIterator elemsBegin, TIterator elemsEnd,
						   std::vector<number> intervals)
{
	sh.clear();
	if(intervals.empty()){
		sh.assign_subset(elemsBegin, elemsEnd, 0);
		return true;
	}

//	access position
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPosition);
	
//	iterate over all elements
	for(TIterator iter = elemsBegin; iter != elemsEnd; ++iter)
	{
		typename TIterator::value_type elem = *iter;
		number quality = FaceQuality(elem, aaPos);
		//...
	}
}
*/
}//	end of namespace

#endif
