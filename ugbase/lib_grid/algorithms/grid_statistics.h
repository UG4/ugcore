/*
 * Copyright (c) 2009-2015:  G-CSC, Goethe University Frankfurt
 * Authors: Martin Stepniewski, Sebastian Reiter
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
template <typename TIterator>
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


template <typename TIterator, typename TAAPos>
void PrintElementEdgeRatios(Grid& grid, TIterator elemsBegin, TIterator elemsEnd,
							TAAPos aaPos)
{
	using namespace std;
	using elem_t = typename PtrToValueType<typename TIterator::value_type>::base_type;
	UG_COND_THROW(elem_t::BASE_OBJECT_ID == GridBaseObjectId::VERTEX || elem_t::BASE_OBJECT_ID == GridBaseObjectId::EDGE,
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
		for(size_t _vfeI = 0; _vfeI < edges.size(); ++_vfeI){ Edge* e = edges[_vfeI];{
			number l = EdgeLength(e, aaPos);
			shortest = min(shortest, l);
			longest = max(longest, l);
		}};

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
