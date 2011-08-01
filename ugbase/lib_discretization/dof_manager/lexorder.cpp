/*
 * cuthill_mckee.cpp
 *
 *  Created on: 21.03.2011
 *      Author: andreasvogel
 */

#include "common/common.h"
#include "lexorder.h"
#include <algorithm>
#include <vector>
#include <queue>
#include <utility>

namespace ug{

/// utility class to provide compare operator for indices based on their position
/**
 * This class is used to provide an ordering for indices. The ordering relation
 * is based on the geometric position.
 *  */




/*

template<int dim>
struct ComparePosition {
///	constructor, passing field with connections for each index
	ComparePosition(const std::vector<std::pair<MathVector<dim>, size_t> >& vPos) : m_vPos(vPos) {}

///	comparison operator
	bool operator() (size_t i,size_t j)
	{
		UG_ASSERT(i < m_vPos.size(), "Invalid index.");
		UG_ASSERT(j < m_vPos.size(), "Invalid index.");
		return ComparePosDim(m_vPos[i].first, m_vPos[j].first);
	}

private:
	///	storage field for teh position of each index
	const std::vector<MathVector<dim> >& m_vPos;

	// Order for 1D
	bool ComparePosDim(MathVector<1> p1, MathVector<1> p2)
	{return (p1[0]<p2[0]);};

	// Order for 2D
	bool ComparePosDim(MathVector<2> p1, MathVector<2> p2)
	{return (p1[0]==p2[0] && p1[1]<p2[1]);};

	// Order for 3D
	bool ComparePosDim(MathVector<3> p1, MathVector<3> p2)
	{return (p1[0]==p2[0] && p1[1]==p2[1] && p1[2]<p2[2]);};
};

// computes ordering using Cuthill-McKee algorithm
template<int dim>
bool ComputeLexicographicOrder(std::vector<size_t>& vNewIndex,
					const std::vector<std::pair<MathVector<dim>, size_t> >& vPos,
                    int order, int sign)
{
//	list of sorted indices
	const size_t numVec = vPos.size();
	std::vector<size_t> vNewOrder(numVec);
	if(vNewOrder.size()< numVec)
		{
			UG_LOG("ERROR in OrderLex: "
					"Out of memory.\n");
			return false;
		}

// sort indices based on their position
	ComparePosition<dim> myCompPos(vPos);
	std::sort(vPos.begin(), vPos.end(), myCompPos);


	for (size_t i=0; i<numVec; ++i)
		vNewOrder[vPos.second] = i;


//	we're done
	return true;
}

*/

}
