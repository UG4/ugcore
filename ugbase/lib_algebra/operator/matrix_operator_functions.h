/**
 * \file lib_algebra/operator/matrix_operator_functions.h
 * This file contains some methods that forward operations on
 * IMatrixOperator to the concrete matrix methods. It can be removed
 * as soon as IMatrixOperator is derived from its Matrix.
 */

#ifndef __H__LIB_ALGEBRA__MATRIX_OPERATOR_FUNCTIONS__
#define __H__LIB_ALGEBRA__MATRIX_OPERATOR_FUNCTIONS__

#include "lib_algebra/operator/interface/operator.h"

namespace ug
{
/*
template <typename X, typename Y, typename M>
void MatAdd(	IMatrixOperator<M, X, Y>& opOut,
				IMatrixOperator<M, X, Y>& op1,
				IMatrixOperator<M, X, Y>& op2)
{
	typedef IMatrixOperator<M, X, Y> MatrixOperator;

	typedef typename MatrixOperator::matrix_type Matrix;

	Matrix& matOut = opOut.get_matrix();
	Matrix& mat1 = op1.get_matrix();
	Matrix& mat2 = op2.get_matrix();

//	matOut = mat1 + mat2;
	MatAdd<X, M>(matOut, mat1, mat2);
}
*/

template <typename X, typename Y, typename M>
void MatIdentity(	MatrixOperator<M, X, Y>& opOut)
{
	PROFILE_FUNC_GROUP("algebra");
	typedef MatrixOperator<M, X, Y> MatrixOperator;

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
void MatAdd( MatrixOperator<M, X, Y>& res, number alpha1, MatrixOperator<M, X, Y>& A1, number alpha2, MatrixOperator<M, X, Y>& A2)
{
	PROFILE_FUNC_GROUP("algebra");
	typedef MatrixOperator<M, X, Y> MatrixOperator;

	typedef typename MatrixOperator::matrix_type Matrix;

	Matrix& matRes = res.get_matrix();
	Matrix& matA1 = A1.get_matrix();
	Matrix& matA2 = A2.get_matrix();
	MatAdd(matRes, alpha1, matA1, alpha2, matA2);
}

template <typename X, typename Y, typename M>
void MatScale( MatrixOperator<M, X, Y>& A, number alpha)
{
	PROFILE_FUNC_GROUP("algebra");
	typedef MatrixOperator<M, X, Y> MatrixOperator;
	typedef typename MatrixOperator::matrix_type Matrix;
	Matrix& matA = A.get_matrix();

	matA.scale(alpha);
}

}//	end of namespace

#endif
