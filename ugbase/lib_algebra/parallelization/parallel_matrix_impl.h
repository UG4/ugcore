/*
 * parallel_matrix_impl.h
 *
 *  Created on: 19.10.2010
 *      Author: A. Vogel
 */

#ifndef __H__LIB_ALGEBRA__PARALLELIZATION__PARALLEL_MATRIX_IMPL__
#define __H__LIB_ALGEBRA__PARALLELIZATION__PARALLEL_MATRIX_IMPL__

#include "parallel_matrix.h"

namespace ug
{

template <typename TMatrix>
bool
ParallelMatrix<TMatrix>::
change_storage_type(ParallelStorageType type)
{
	// check that communicator exists
	if(m_pCommunicator == NULL)
	{
		UG_LOG("No communicator set. Cannot change storage type.\n");
		return false;
	}

	// can only change if current state is defined
	if(has_storage_type(PST_UNDEFINED))
	{
		UG_LOG("Current Storage Type is undefined. Cannot change storage type.\n");
		return false;
	}

	// if already in that type
	if(has_storage_type(type))
		return true;

//	todo: Implement more.
	UG_LOG("Currently no storage conversion supported.");
	return false;
}

// calculate res = A x
template <typename TMatrix>
template<typename TPVector>
bool
ParallelMatrix<TMatrix>::
apply(TPVector &res, const TPVector &x) const
{
//	check types combinations
	int type = -1;
	if(this->has_storage_type(PST_ADDITIVE)   && x.has_storage_type(PST_CONSISTENT)) type = 0;
	if(this->has_storage_type(PST_CONSISTENT) && x.has_storage_type(PST_ADDITIVE)) type = 1;
	if(this->has_storage_type(PST_CONSISTENT) && x.has_storage_type(PST_CONSISTENT)) type = 2;

//	if no admissible type is found, return error
	if(type == -1)
	{
		UG_LOG("ERROR in 'ParallelMatrix::apply' (b = A*x): Wrong storage type of Matrix/Vector: Possibilities are:\n"
				"    - A is PST_ADDITIVE and x is PST_CONSISTENT\n"
				"    - A is PST_CONSISTENT and x is PST_ADDITIVE\n");
		return false;
	}

//	apply on single process vector
	if(!TMatrix::apply(res, x))
		return false;

//	set outgoing vector to additive storage
	switch(type)
	{
		case 0: res.set_storage_type(PST_ADDITIVE); break;
		case 1: res.set_storage_type(PST_ADDITIVE); break;
		case 2: res.set_storage_type(PST_CONSISTENT); break;
	}

//	we're done.
	return true;
}

// calculate res = A.T x
template <typename TMatrix>
template<typename TPVector>
bool
ParallelMatrix<TMatrix>::
apply_transposed(TPVector &res, const TPVector &x) const
{
//	check types combinations
	int type = -1;
	if(this->has_storage_type(PST_ADDITIVE)   && x.has_storage_type(PST_CONSISTENT)) type = 0;
	if(this->has_storage_type(PST_CONSISTENT) && x.has_storage_type(PST_ADDITIVE)) type = 1;
	if(this->has_storage_type(PST_CONSISTENT) && x.has_storage_type(PST_CONSISTENT)) type = 2;

//	if no admissible type is found, return error
	if(type == -1)
	{
		UG_LOG("ERROR in 'ParallelMatrix::apply_transposed' (b = A^T*x): Wrong storage type of Matrix/Vector: Possibilities are:\n"
				"    - A is PST_ADDITIVE and x is PST_CONSISTENT\n"
				"    - A is PST_CONSISTENT and x is PST_ADDITIVE\n");
		return false;
	}

//	apply on single process vector
	if(!TMatrix::apply_transposed(res, x))
		return false;

//	set outgoing vector to additive storage
	switch(type)
	{
		case 0: res.set_storage_type(PST_ADDITIVE); break;
		case 1: res.set_storage_type(PST_ADDITIVE); break;
		case 2: res.set_storage_type(PST_CONSISTENT); break;
	}

//	we're done.
	return true;
}

// calculate res -= A x
template <typename TMatrix>
template<typename TPVector>
bool
ParallelMatrix<TMatrix>::
matmul_minus(TPVector &res, const TPVector &x) const
{
//	check types combinations
	int type = -1;
	if(this->has_storage_type(PST_ADDITIVE) && x.has_storage_type(PST_CONSISTENT) && res.has_storage_type(PST_ADDITIVE)) type = 0;

//	if no admissible type is found, return error
	if(type == -1)
	{
		UG_LOG("ERROR in 'ParallelMatrix::matmul_minus' (b -= A*x): Wrong storage type of Matrix/Vector: Possibilities are:\n"
				"    - A is PST_ADDITIVE and x is PST_CONSISTENT and b is PST_ADDITIVE\n");
		return false;
	}

//	apply on single process vector
	if(!TMatrix::matmul_minus(res, x))
		return false;

//	set outgoing vector to additive storage (it could have been PST_UNIQUE before)
	switch(type)
	{
		case 0: res.set_storage_type(PST_ADDITIVE); break;
	}

//	we're done.
	return true;
}




} // end namespace ug

#endif /* __H__LIB_ALGEBRA__PARALLELIZATION__PARALLEL_MATRIX_IMPL__ */
