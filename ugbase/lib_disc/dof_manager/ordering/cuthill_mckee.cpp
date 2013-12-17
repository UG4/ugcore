/*
 * cuthill_mckee.cpp
 *
 *  Created on: 21.03.2011
 *      Author: andreasvogel
 */

#include "common/common.h"
#include "cuthill_mckee.h"
#include <algorithm>
#include <vector>
#include <queue>
#include "common/profiler/profiler.h"
#include "lib_disc/domain.h"

namespace ug{

/// help class to provide compare operator for indices based on their degree
/**
 * This class is used to provide an ordering for indices. The ordering relation
 * is based on the connectivity-degree, i.e. on the number of connections the
 * index has. The indices with less connections are ordered first.
 */
struct CompareDegree {
///	constructor, passing field with connections for each index
	CompareDegree(const std::vector<std::vector<size_t> >& vInfo) : m_vCon(vInfo) {}

///	comparison operator
	bool operator() (size_t i,size_t j)
	{
		UG_ASSERT(i < m_vCon.size(), "Invalid index.");
		UG_ASSERT(j < m_vCon.size(), "Invalid index.");
		return (m_vCon[i].size() < m_vCon[j].size());
	}

private:
///	storage field of connections of each index
	const std::vector<std::vector<size_t> >& m_vCon;
};

// computes ordering using Cuthill-McKee algorithm
void ComputeCuthillMcKeeOrder(std::vector<size_t>& vNewIndex,
                              std::vector<std::vector<size_t> >& vvConnection,
                              bool bReverse)
{
	PROFILE_FUNC();
//	list of sorted (will be filled later)
	std::vector<size_t> vNewOrder;

//	create flag list to remember already handled indices
	std::vector<bool> vHandled(vvConnection.size(), false);

//	Sort neighbours by degree (i.e. by number of neighbours those have)
	CompareDegree myCompDegree(vvConnection);
	for(size_t i = 0; i < vvConnection.size(); ++i)
	{
	//	indices with no adjacent indices are marked as handled (and skipped)
		if(vvConnection[i].size() == 0){
			vHandled[i] = true;
		}
	//	sort adjacent index by degree
		else {
			std::sort(	vvConnection[i].begin(),
			          	vvConnection[i].end(), myCompDegree);
		}
	}

//	start with first index
	size_t firstNonHandled = 0;

//	handle all indices
	while(true)
	{
	//	find first unhandled index
		size_t i_notHandled = firstNonHandled;
		for(; i_notHandled < vHandled.size(); ++i_notHandled){
			if(!vHandled[i_notHandled]) {firstNonHandled = i_notHandled; break;}
		}
	//	check if one unhandled vertex left
		if(i_notHandled == vHandled.size()) break;

	//	Find node with smallest degree for all remaining indices
		size_t start = firstNonHandled;
		for(size_t i = start+1; i < vHandled.size(); ++i)
		{
			if(!vHandled[i] && vvConnection[i].size() < vvConnection[start].size())
				start = i;
		}

	//	Add start vertex to mapping
		vNewOrder.push_back(start);
		vHandled[start] = true;

	//	Create queue of adjacent vertices
		std::queue<size_t> qAdjacent;
		for(size_t i = 0; i < vvConnection[start].size(); ++i)
		{
			const size_t ind = vvConnection[start][i];

			if(!vHandled[ind] && ind != start)
				qAdjacent.push(ind);
		}

	//	add adjacent vertices to mapping
		while(!qAdjacent.empty())
		{
		//	get next index
			const size_t front = qAdjacent.front();

		//	if not handled
			if(!vHandled[front])
			{
			//	Add to mapping
				vNewOrder.push_back(front);
				vHandled[front] = true;

			//	add adjacent to queue
				for(size_t i = 0; i < vvConnection[front].size(); ++i)
				{
					const size_t ind = vvConnection[front][i];

					if(!vHandled[ind] && ind != front)
						qAdjacent.push(ind);
				}
			}

		//	pop index
			qAdjacent.pop();
		}
	}

// 	Create list of mapping
	vNewIndex.clear(); vNewIndex.resize(vvConnection.size(), (size_t)-1);

//	write new indices into out array
	size_t cnt = 0;
	if(bReverse)
	{
		for(size_t oldInd = 0; oldInd < vvConnection.size(); ++oldInd)
		{
		//	skip non-sorted indices
			if(vvConnection[oldInd].size() == 0) continue;

		//	get old index
			UG_ASSERT(cnt < vNewOrder.size(), "cnt: "<<cnt<<", ordered: "<<vNewOrder.size())
			const size_t newInd = vNewOrder[vNewOrder.size() - 1 - cnt]; ++cnt;
			UG_ASSERT(newInd < vNewIndex.size(), "newInd: "<<newInd<<", size: "<<vNewIndex.size())

		//	set new index to order
			vNewIndex[newInd] = oldInd;
		}
	}
	else
	{
		for(size_t oldInd = 0; oldInd < vvConnection.size(); ++oldInd)
		{
		//	skip non-sorted indices
			if(vvConnection[oldInd].size() == 0) continue;

		//	get old index
			UG_ASSERT(cnt < vNewOrder.size(), "cnt: "<<cnt<<", ordered: "<<vNewOrder.size())
			const size_t newInd = vNewOrder[cnt++];
			UG_ASSERT(newInd < vNewIndex.size(), "newInd: "<<newInd<<", size: "<<vNewIndex.size())

		//	set new index to order
			vNewIndex[newInd] = oldInd;
		}
	}

//	check if all ordered indices have been written
	if(cnt != vNewOrder.size())
		UG_THROW("OrderCuthillMcKee: No all indices sorted, that must be sorted: "
				<<cnt<<" written, but should write: "<<vNewOrder.size());

//	fill non-sorted indices
	for(size_t i = 1; i < vNewIndex.size(); ++i)
	{
		if(vNewIndex[i] == (size_t)-1) vNewIndex[i] = vNewIndex[i-1] + 1;
	}
}

void OrderCuthillMcKee(DoFDistribution& dofDistr, bool bReverse)
{
	PROFILE_FUNC();
//	get adjacency graph
	std::vector<std::vector<size_t> > vvConnection;
	try{
		dofDistr.get_connections(vvConnection);
	}
	UG_CATCH_THROW("OrderCuthillMcKee: No adjacency graph available.");

//	get mapping for cuthill-mckee order
	std::vector<size_t> vNewIndex;
	ComputeCuthillMcKeeOrder(vNewIndex, vvConnection, bReverse);

//	reorder indices
	dofDistr.permute_indices(vNewIndex);
}

template <typename TDomain>
void OrderCuthillMcKee(ApproximationSpace<TDomain>& approxSpace, bool bReverse)
{
	std::vector<SmartPtr<DoFDistribution> > vDD = approxSpace.dof_distributions();

	for(size_t i = 0; i < vDD.size(); ++i)
		OrderCuthillMcKee(*vDD[i], bReverse);
}

#ifdef UG_DIM_1
template void OrderCuthillMcKee<Domain1d>(ApproximationSpace<Domain1d>& approxSpace, bool bReverse);
#endif
#ifdef UG_DIM_2
template void OrderCuthillMcKee<Domain2d>(ApproximationSpace<Domain2d>& approxSpace, bool bReverse);
#endif
#ifdef UG_DIM_3
template void OrderCuthillMcKee<Domain3d>(ApproximationSpace<Domain3d>& approxSpace, bool bReverse);
#endif

}
