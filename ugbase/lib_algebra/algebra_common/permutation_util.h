/*
 * permutation_util.h
 *
 *  Created on: 02.12.2013
 *      Author: mrupp
 */

#ifndef __H__UG__CPU_ALGEBRA__PERMUTATION_UTIL__
#define __H__UG__CPU_ALGEBRA__PERMUTATION_UTIL__


#include "common/profiler/profiler.h"
#include "common/error.h"
#include "lib_disc/dof_manager/ordering/cuthill_mckee.h"
#include <vector>

namespace ug{
/**
 * Function to return a permutation of a matrix
 * @param[out] PA the permuted matrix PA(perm[r], perm[c]) = A(r, c)
 * @param[in] A the input matrix
 * @param[in] perm array mapping i -> perm[i]
 */
template<typename TMatrix>
static void SetMatrixAsPermutation(TMatrix &PA, const TMatrix &A, std::vector<size_t> &perm)
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
static void SetVectorAsPermutation(TVector &Pv, const TVector &v, std::vector<size_t> &perm)
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
bool GetInversePermutation(const std::vector<size_t> &perm, std::vector<size_t> &invPerm);

/**
 * @param mat 			A sparse matrix
 * @param newIndex		the cuthill-mckee ordered new indices
 */
template<typename TSparseMatrix>
void GetCuthillMcKeeOrder(const TSparseMatrix &mat, std::vector<size_t> &newIndex)
{
	std::vector<std::vector<size_t> > neighbors;
	neighbors.resize(mat.num_rows());

	for(size_t i=0; i<mat.num_rows(); i++)
	{
		for(typename TSparseMatrix::const_row_iterator i_it = mat.begin_row(i); i_it != mat.end_row(i); ++i_it)
			neighbors[i].push_back(i_it.index());
	}

	ComputeCuthillMcKeeOrder(newIndex, neighbors, false);
}
/// @}
} // end namespace ug


#endif
