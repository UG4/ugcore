/*
 * Copyright (c) 2011-2014:  G-CSC, Goethe University Frankfurt
 * Author: Konstantinos Xylouris
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

/**
 * \file lib_algebra/operator/matrix_operator_functions.h
 * This file contains some methods that forward operations on
 * IMatrixOperator to the concrete matrix methods. It can be removed
 * as soon as IMatrixOperator is derived from its Matrix.
 */

#ifndef __H__LIB_ALGEBRA__MATRIX_OPERATOR_FUNCTIONS__
#define __H__LIB_ALGEBRA__MATRIX_OPERATOR_FUNCTIONS__

#include "lib_algebra/operator/interface/matrix_operator.h"

namespace ug {


template <typename X, typename Y, typename M>
void MatIdentity(	MatrixOperator<M, X, Y>& opOut)
{
	PROFILE_FUNC_GROUP("algebra");
	using MatrixOperator = MatrixOperator<M, X, Y>;

	using Matrix = typename MatrixOperator::matrix_type;

	Matrix& matOut = opOut.get_matrix();
	size_t numRows = matOut.num_rows();
	size_t numCols = matOut.num_cols();

	matOut.resize_and_clear(numRows, numCols);

	for(size_t i = 0; i < numRows; ++i)
		matOut(i, i) = 1.0;
}


template <typename X, typename Y, typename M>
void MatAdd( MatrixOperator<M, X, Y>& res, number alpha1, MatrixOperator<M, X, Y>& A1, number alpha2, MatrixOperator<M, X, Y>& A2)
{
	PROFILE_FUNC_GROUP("algebra");
	using MatrixOperator = MatrixOperator<M, X, Y>;

	using Matrix = typename MatrixOperator::matrix_type;

	Matrix& matRes = res.get_matrix();
	Matrix& matA1 = A1.get_matrix();
	Matrix& matA2 = A2.get_matrix();
	MatAdd(matRes, alpha1, matA1, alpha2, matA2);
}

template <typename X, typename Y, typename M>
void MatScale( MatrixOperator<M, X, Y>& A, number alpha)
{
	PROFILE_FUNC_GROUP("algebra");
	using MatrixOperator = MatrixOperator<M, X, Y>;
	using Matrix = typename MatrixOperator::matrix_type;
	Matrix& matA = A.get_matrix();

	matA.scale(alpha);
}

template <typename X, typename Y, typename M>
void MatTranspose( MatrixOperator<M, X, Y>& AT,  MatrixOperator<M, X, Y>& A)
{
	PROFILE_FUNC_GROUP("algebra");
	using MatrixOperator = MatrixOperator<M, X, Y>;
	using Matrix = typename MatrixOperator::matrix_type;

	Matrix& matA = A.get_matrix();
	Matrix& matAT = AT.get_matrix();

	matAT.set_as_transpose_of(matA);
}

}//	end of namespace

#endif
