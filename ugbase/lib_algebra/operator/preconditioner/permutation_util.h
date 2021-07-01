/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
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

#ifndef __H__UG__CPU_ALGEBRA__PERMUTATION_UTIL__DUPLICATE
#define __H__UG__CPU_ALGEBRA__PERMUTATION_UTIL__DUPLICATE

#include "common/common.h"
#include "common/profiler/profiler.h"
#include "common/error.h"
#include <vector>

namespace ug{
/**
 * Function to return a permutation of a matrix
 * @param[out] PA the permuted matrix PA(perm[r], perm[c]) = A(r, c)
 * @param[in] A the input matrix
 * @param[in] perm array mapping i -> perm[i]
 */
template<typename TMatrix>
static void permute_matrix(TMatrix &PA, const TMatrix &A, const std::vector<size_t> &perm)
{
	PROFILE_FUNC_GROUP("algebra");
	PA.resize_and_clear(A.num_rows(), A.num_cols());

	for(size_t r=0; r<A.num_rows(); r++)
	{
		size_t Pr=perm[r];
		for(typename TMatrix::const_row_iterator it = A.begin_row(r); it != A.end_row(r); ++it)
		{
			size_t Pc = perm[it.index()];
			PA(Pr, Pc) = it.value();
		}
	}
}

/**
 * Function to compute a permutation of a vector
 * @param[out] Pv the permuted vector: Pv[perm[i]] = v[i]
 * @param[in] v the input vector
 * @param[in] perm array mapping oldindices i -> new index perm[i]
 */
template<typename TVector>
static void permute_vector(TVector &Pv, const TVector &v, const std::vector<size_t> &perm)
{
	if(Pv.size() != v.size()) Pv.resize(v.size());
	for(size_t i=0; i<v.size(); i++)
		Pv[ perm[i] ] = v[i];
}


/**
 * Function to calculate the inverse of a permutation
 * @param[in] perm array mapping  i -> perm[i]
 * @param[out] invPerm array mapping  i -> invPerm[i] so that perm[ invPerm[i] ] = i
 * @return true if perm[i] == i forall i.
 */

bool inverse_permutation(const std::vector<size_t> &perm, std::vector<size_t> &invPerm)
{
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


/// @}
} // end namespace ug


#endif
