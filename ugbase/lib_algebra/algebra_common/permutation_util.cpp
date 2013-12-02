/*
 * permutation_util.cpp
 *
 *  Created on: 02.12.2013
 *      Author: mrupp
 */

#include "permutation_util.h"


namespace ug{


bool GetInversePermutation(const std::vector<size_t> &perm, std::vector<size_t> &invPerm)
{
	invPerm.resize(perm.size());
	bool bId = true;
	for(size_t i=0; i<perm.size(); i++) invPerm[i] = (size_t) (-1);

	for(size_t i=0; i<perm.size(); i++)
	{
		UG_COND_THROW(invPerm[perm[i]] != (size_t) (-1), "not a bijective permutation!");
		bId = bId && perm[i] == i;
		invPerm[perm[i]] = i;
	}
	return bId;
}


}
