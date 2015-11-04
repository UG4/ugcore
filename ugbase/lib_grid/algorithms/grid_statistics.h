//	created by Martin Stepniewski, Sebastian Reiter
//	mastep@gmx.de
//	y09 m11 d11

#ifndef __H__UG__GRID_STATISTICS__
#define __H__UG__GRID_STATISTICS__

#include <algorithm>
#include <limits>
#include <vector>
#include "lib_grid/lg_base.h"

#ifdef UG_PARALLEL
	#include "lib_grid/parallelization/distributed_grid.h"
	#include "pcl/pcl_process_communicator.h"
#endif

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
											  std::vector<double>& offsets,
											  Grid::VertexAttachmentAccessor<APosition>& aaPos);

////////////////////////////////////////////////////////////////////////
///	assigns a subset based on the quality of the given element.
/**
 * Currently only faces are supported.
 *
 * \param intervals contains the intervals which define into which subset
 *					an element goes. Numbers have to be sorted, starting at
 *					0 and ending at 1 (0 and 1 should be contained in intervals).
 */
template <class TIterator>
bool AssignSubsetsByQuality(Grid& grid, SubsetHandler& sh,
						   TIterator elemsBegin, TIterator elemsEnd,
						   std::vector<number> intervals)
{
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

		size_t newInd = -1;
		for(size_t i = 0; i < intervals.size(); ++i){
			if(intervals[i] < quality)
				newInd = i;
			else
				break;
		}

		sh.assign_subset(elem, newInd);
	}
	return true;
}


template <class TIterator, class TAAPos>
void PrintElementEdgeRatios(Grid& grid, TIterator elemsBegin, TIterator elemsEnd,
							TAAPos aaPos)
{
	using namespace std;
	typedef typename PtrToValueType<typename TIterator::value_type>::base_type	elem_t;
	UG_COND_THROW(elem_t::BASE_OBJECT_ID == VERTEX || elem_t::BASE_OBJECT_ID == EDGE,
				  "Can't evaluate anisotropy statistics for vertices or edges.");

	Grid::edge_traits::secure_container	edges;

	number minRatio = 1;
	number maxRatio = 0;
	number avRatio = 0;
	vector<number>	ratios;
	for(TIterator i_elem = elemsBegin; i_elem != elemsEnd; ++i_elem){
		elem_t* elem = *i_elem;
		
		#ifdef UG_PARALLEL
			if(grid.distributed_grid_manager()->is_ghost(elem))
				continue;
		#endif

		grid.associated_elements(edges, elem);
		number shortest = numeric_limits<double>::max();
		number longest = 0;
		for_each_in_vec(Edge* e, edges){
			number l = EdgeLength(e, aaPos);
			shortest = min(shortest, l);
			longest = max(longest, l);
		}end_for;

		number ratio = 0;
		if(longest > 0)
			ratio = shortest / longest;

		minRatio = min(minRatio, ratio);
		maxRatio = max(maxRatio, ratio);
		avRatio += ratio;
		ratios.push_back(ratio);
	}

	int num = (int)ratios.size();
	#ifdef UG_PARALLEL
		pcl::ProcessCommunicator com;
		minRatio = com.allreduce(minRatio, PCL_RO_MIN);
		maxRatio = com.allreduce(maxRatio, PCL_RO_MAX);
		num = com.allreduce(num, PCL_RO_SUM);
		avRatio = com.allreduce(avRatio, PCL_RO_SUM);
	#endif
	
	if(num == 0){
		UG_LOG("---\n");
	}
	else{
		avRatio /= (number)num;
		UG_LOG("min: " << minRatio << ",  max: " << maxRatio << ",  av: " << avRatio);
		
		if(num > 1){
			number sdSum = 0;
			for(size_t i = 0; i < ratios.size(); ++i)
				sdSum += sq(avRatio - ratios[i]);

			#ifdef UG_PARALLEL
				sdSum = com.allreduce(sdSum, PCL_RO_SUM);
			#endif

			number sd = sqrt(sdSum / ((number)num - 1));
			UG_LOG(",  sd: " << sd);
		}
		UG_LOG(endl);
	}
}

}//	end of namespace

#endif
