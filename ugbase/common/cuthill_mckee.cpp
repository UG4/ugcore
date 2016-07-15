#include "common/common.h"
#include "cuthill_mckee.h"
#include <algorithm>
#include <vector>
#include <queue>
#include "common/profiler/profiler.h"

namespace ug{

bool CheckPermutationBijective(const std::vector<size_t> &perm)
{
	std::vector<size_t> invPerm;
	invPerm.resize(perm.size());
	bool bId = true;
	for(size_t i=0; i<perm.size(); i++) invPerm[i] = (size_t) (-1);

	for(size_t i=0; i<perm.size(); i++)
	{
		UG_COND_THROW(invPerm[perm[i]] != (size_t) (-1), "not a bijective permutation "
			"(double mapping to index " << perm[i] << " by indices " << invPerm[perm[i]] << " and " << i << ")!");
		bId = bId && perm[i] == i;
		invPerm[perm[i]] = i;
	}
	return bId;
}

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
                              bool bReverse,
							  bool bPreserveConsec)
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
	//	edit (mbreit, 10-11-2015): using stable_sort here because implementations
	//	of std::sort seem to vary among different platforms with regard to the
	//	sorting outcome of entries where myCompDegree is "=", causing differences
	//	in convergence behavior in ILU-T
		else {
			std::stable_sort(	vvConnection[i].begin(),
			          	vvConnection[i].end(), myCompDegree);
		}
	}

	// also sort vvConnection itself, this is extremely useful if there are many
	// identity rows, as finding the "start" index again and again will take
	// REALLY much time in that case
	std::vector<size_t> sorting(vvConnection.size());
	size_t szSort = sorting.size();
	for (size_t i = 0; i < szSort; ++i)
		sorting[i] = i;
	std::stable_sort(sorting.begin(), sorting.end(), myCompDegree);

//	start with first index
	size_t firstNonHandled = 0;

//	handle all indices
	while(true)
	{
	//	find first unhandled index
		size_t i_notHandled = firstNonHandled;
		for(; i_notHandled < szSort; ++i_notHandled){
			if(!vHandled[sorting[i_notHandled]]) {firstNonHandled = i_notHandled; break;}
		}

	//	check if one unhandled vertex left
		if(i_notHandled == szSort) break;

/* // This is no longer necessary as firstUnhandled automatically
      has smallest degree due to sorting.
	//	Find node with smallest degree for all remaining indices
		size_t start = firstNonHandled;
		for(size_t i = start+1; i < vHandled.size(); ++i)
		{
			if(!vHandled[i] && vvConnection[i].size() < vvConnection[start].size())
				start = i;
		}
*/
		size_t start = sorting[firstNonHandled];

	//	Create queue of adjacent vertices
		std::queue<size_t> qAdjacent;
		qAdjacent.push(start);

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

					if(!vHandled[ind])
						qAdjacent.push(ind);
				}
			}

		//	pop index
			qAdjacent.pop();
		}
	}

// 	Create list of mapping
	vNewIndex.clear(); vNewIndex.resize(vvConnection.size(), (size_t)-1);

	if (bPreserveConsec)
	{
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
			UG_THROW("OrderCuthillMcKee: Not all indices sorted that must be sorted: "
					<<cnt<<" written, but should write: "<<vNewOrder.size());

	//	fill non-sorted indices (preserving consecutive indexing)
		for(size_t i = 1; i < vNewIndex.size(); ++i)
		{
			if(vNewIndex[i] == (size_t)-1) vNewIndex[i] = vNewIndex[i-1] + 1;
		}
	}
	else
	{
		size_t newOrdSz = vNewOrder.size();
		if (bReverse)
		{
			for (size_t i = 0; i < newOrdSz; ++i)
			{
				size_t oldInd = vNewOrder[newOrdSz - 1 - i];
				vNewIndex[oldInd] = i;
			}
		}
		else
		{
			for (size_t i = 0; i < newOrdSz; ++i)
			{
				size_t oldInd = vNewOrder[i];
				vNewIndex[oldInd] = i;
			}
		}

		// move unconnected indices to the bottom of the ordering
		for (size_t i = 0; i < vNewIndex.size(); ++i)
		{
			if (vNewIndex[i] == (size_t)-1)
				vNewIndex[i] = newOrdSz++;
		}
	}

	//CheckPermutationBijective(vNewIndex);
}

}
