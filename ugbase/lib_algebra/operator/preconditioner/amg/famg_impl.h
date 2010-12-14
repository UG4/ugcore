/**
 * \file amg_impl.h
 *
 * \author Martin Rupp
 *
 * \date 16.11.2010
 *
 * implementation file for amg
 *
 * Goethe-Center for Scientific Computing 2009-2010.
 */

#ifndef __H__LIB_ALGEBRA__AMG__FAMG_IMPL_H__
#define __H__LIB_ALGEBRA__AMG__FAMG_IMPL_H__

//#include "sparsematrix_util.h"

//#include "famg_nodeinfo.h"
#include "stopwatch.h"
#include "common/assert.h"
//#include "maxheap.h"

namespace ug{

void c_create_AMG_level(SparseMatrix<double> &AH, SparseMatrix<double> &R, const SparseMatrix<double> &A,
		SparseMatrix<double> &P, cAMG_helper &amghelper, int level);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// createFAMGLevel:
//-------------------------
/**
 * create FAMG matrix R, P, and AH = R A P
 * \param AH
 * \param R
 * \param A
 * \param P
 * \param level
 */
template<typename TAlgebra>
void famg<TAlgebra>::create_AMG_level(matrix_type &AH, SparseMatrix<double> &R, const matrix_type &A,
		SparseMatrix<double> &P, int level)
{
	c_create_AMG_level(AH, R, A, P, amghelper, level);

}

template<typename TAlgebra>
famg<TAlgebra>::famg() : amg_base<TAlgebra>()
{
}

template<typename TAlgebra>
famg<TAlgebra>::~famg()
{

}

template<typename TAlgebra>
void famg<TAlgebra>::tostring() const
{
	amg_base<TAlgebra>::tostring();
	UG_LOG("FAMG Preconditioner:\n");
}


} // namespace ug

#endif //  __H__LIB_ALGEBRA__AMG__FAMG_IMPL_H__
