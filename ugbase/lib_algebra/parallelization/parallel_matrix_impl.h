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
		{UG_LOG("No communicator set. Cannot change storage type.\n"); return false;}

	// can only change if current state is defined
	if(has_storage_type(PST_UNDEFINED)) return false;

	// if already in that type
	if(has_storage_type(type)) return true;

//	todo: currently only one storage type allowed. Implement more.
	return false;
}

} // end namespace ug

#endif /* __H__LIB_ALGEBRA__PARALLELIZATION__PARALLEL_MATRIX_IMPL__ */
