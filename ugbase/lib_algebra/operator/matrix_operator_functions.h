/**	This file contains some methods that forward operations on
 * IMatrixOperator to the concrete matrix methods. It can be removed
 * as soon as IMatrixOperator is derived from its Matrix.
 */

#ifndef __H__LIB_ALGEBRA__MATRIX_OPERATOR_FUNCTIONS__
#define __H__LIB_ALGEBRA__MATRIX_OPERATOR_FUNCTIONS__

#include "operator_interface.h"

namespace ug
{
/*
template <typename X, typename Y, typename M>
void MatAdd(	IMatrixOperator<X, Y, M>& opOut,
				IMatrixOperator<X, Y, M>& op1,
				IMatrixOperator<X, Y, M>& op2)
{
	typedef IMatrixOperator<X, Y, M> MatrixOperator;

	typedef typename MatrixOperator::matrix_type Matrix;

	Matrix& matOut = opOut.get_matrix();
	Matrix& mat1 = op1.get_matrix();
	Matrix& mat2 = op2.get_matrix();

//	matOut = mat1 + mat2;
	MatAdd<X, M>(matOut, mat1, mat2);
}
*/

template <typename X, typename Y, typename M>
void MatIdentity(	IMatrixOperator<X, Y, M>& opOut)
{
	typedef IMatrixOperator<X, Y, M> MatrixOperator;

	typedef typename MatrixOperator::matrix_type Matrix;

	Matrix& matOut = opOut.get_matrix();
	size_t numRows = matOut.num_rows();
	size_t numCols = matOut.num_cols();

	matOut.resize(0, 0);
	matOut.resize(numRows, numCols);

	for(size_t i = 0; i < numRows; ++i)
		matOut(i, i) = 1.0;
}


template <typename X, typename Y, typename M>
void MatAdd( IMatrixOperator<X, Y, M>& res, number alpha1, const IMatrixOperator<X, Y, M>& A1, number alpha2, const IMatrixOperator<X, Y, M>& A2)
{
	typedef IMatrixOperator<X, Y, M> MatrixOperator;

	typedef typename MatrixOperator::matrix_type Matrix;

	Matrix& matRes = res.get_matrix();
	Matrix& matA1 = A1.get_matrix();
	Matrix& matA2 = A2.get_matrix();

	MatAdd(matRes, alpha1, matA1, alpha2, matA2);
}


}//	end of namespace

#endif
