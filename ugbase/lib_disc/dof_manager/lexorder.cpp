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


// Order for 1D
template<int dim>
bool ComparePosDim(const std::pair<MathVector<dim>, size_t> &p1,
                   const std::pair<MathVector<dim>, size_t> &p2)
{return false;}

template<>
bool ComparePosDim<1>(const std::pair<MathVector<1>, size_t> &p1,
                      const std::pair<MathVector<1>, size_t> &p2)
{
	return (p1.first[0]<p2.first[0]);
};

// Order for 2D
template<>
bool ComparePosDim<2>(const std::pair<MathVector<2>, size_t> &p1,
                      const std::pair<MathVector<2>, size_t> &p2)
{
	if (p1.first[1]==p2.first[1]) {
		return (p1.first[0] < p2.first[0]);
	}
	else {
		return (p1.first[1] < p2.first[1]);
	}
};

// Order for 3D
template<>
bool ComparePosDim<3>(const std::pair<MathVector<3>, size_t> &p1,
                      const std::pair<MathVector<3>, size_t> &p2)
{
	if (p1.first[2]==p2.first[2]){
		if (p1.first[1]==p2.first[1]) {
			return p1.first[0] < p2.first[0];
		}
		else {
			return (p1.first[1] < p2.first[1]);
		}
	}
	else{
		return (p1.first[2] < p2.first[2]);
	}
};


template<int dim>
void ComputeLexicographicOrder(std::vector<size_t>& vNewIndex,
                               std::vector<std::pair<MathVector<dim>, size_t> >& vPos)
{
//	a) order all indices
	if(vNewIndex.size() == vPos.size()){
	//  sort indices based on their position
		std::sort(vPos.begin(), vPos.end(), ComparePosDim<dim>);

	//	write mapping
		for (size_t i=0; i < vPos.size(); ++i)
			vNewIndex[vPos[i].second] = i;
	}
//	b) only some indices to order
	else{
	//	create copy of pos
		std::vector<std::pair<MathVector<dim>, size_t> > vPosOrig(vPos);

	//  sort indices based on their position
		std::sort(vPos.begin(), vPos.end(), ComparePosDim<dim>);

	//	write mapping
		for (size_t i=0; i < vNewIndex.size(); ++i)
			vNewIndex[i] = i;
		for (size_t i=0; i < vPos.size(); ++i)
			vNewIndex[vPos[i].second] = vPosOrig[i].second;
	}
}


template void ComputeLexicographicOrder<1>(std::vector<size_t>& vNewIndex, std::vector<std::pair<MathVector<1>, size_t> >& vPos);
template void ComputeLexicographicOrder<2>(std::vector<size_t>& vNewIndex, std::vector<std::pair<MathVector<2>, size_t> >& vPos);
template void ComputeLexicographicOrder<3>(std::vector<size_t>& vNewIndex, std::vector<std::pair<MathVector<3>, size_t> >& vPos);

}
