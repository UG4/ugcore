/*
 * Copyright (c) 2014-2015:  G-CSC, Goethe University Frankfurt
 * Author: Raphael Prohl
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

#ifndef __H__UG__COMMON__MATH_TENSOR_FUNCTIONS_COMMON_IMPL__
#define __H__UG__COMMON__MATH_TENSOR_FUNCTIONS_COMMON_IMPL__

#include <cmath>
#include <iostream>
#include <iomanip>
#include <cassert>
#include "math_tensor.h"
#include "common/assert.h"
#include "common/static_assert.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
// transformation of a tensors
////////////////////////////////////////////////////////////////////////////////

// transformation of a tensor of rank 2 (R^dim x R^dim) into a vector (R^(dim^2))
template <std::size_t TDim, std::size_t TDimSQ>
void
Tens2ToVec(DenseVector<FixedArray1<number, TDimSQ> > &vec,
		const MathMatrix<TDim, TDim>& tens2)
{
	static constexpr size_t dim = TDim;
	static constexpr size_t dimSQ = TDimSQ;
	if (dimSQ != dim * dim)
		UG_THROW("Tens2ToVec::Invalid dimensions: " << dimSQ << " is not the square of " << dim << "! \n");

	for (size_t i = 0; i < dim; ++i)
		for (size_t j = 0; j < dim; ++j){
			vec[i*dim + j] = tens2[i][j];
		}
}

//	transformation of a vector (R^(dim^2)) into a tensor of rank 2 (R^dim x R^dim)
template <std::size_t TDim, std::size_t TDimSQ>
void
VecToTens2(MathMatrix<TDim, TDim>& tens2,
		const DenseVector<FixedArray1<number, TDimSQ> > &vec)
{
	static constexpr size_t dim = TDim;
	static constexpr size_t dimSQ = TDimSQ;
	if (dimSQ != dim * dim)
		UG_THROW("VecToTens2::Invalid dimensions: " << dimSQ << " is not the square of " << dim << "! \n");

	for (size_t i = 0; i < dim; ++i)
		for (size_t j = 0; j < dim; ++j){
			tens2[i][j] = vec[i*dim + j];
		}
}

//	transformation of a tensor of rank 4 (R^dim x R^dim x R^dim x R^dim)
//	into a matrix (R^(dim^2) x R^(dim^2))
template <std::size_t TDim, std::size_t TDimSQ>
void
Tens4ToMat(DenseMatrix<FixedArray2<number, TDimSQ, TDimSQ> > &mat,
		const MathTensor4<TDim, TDim, TDim, TDim>& tens4)
{
	static constexpr size_t dim = TDim;
	static constexpr size_t dimSQ = TDimSQ;
	if (dimSQ != dim * dim)
		UG_THROW("Tens4ToMat::Invalid dimensions: " << dimSQ << " is not the square of " << dim << "! \n");

	for (size_t i = 0; i < dim; ++i)
		for (size_t j = 0; j < dim; ++j)
			for (size_t k = 0; k < dim; ++k)
				for (size_t l = 0; l < dim; ++l){
					mat(i*dim + j, k*dim + l) = tens4[i][j][k][l];
				}
}


// 	transformation of a matrix (R^(dim^2) x R^(dim^2))
//	into a tensor of rank 4 (R^dim x R^dim x R^dim x R^dim)
template <std::size_t TDim, std::size_t TDimSQ>
void
MatToTens4(MathTensor4<TDim, TDim, TDim, TDim>& tens4,
		const DenseMatrixInverse<FixedArray2<number, TDimSQ, TDimSQ> > &mat)
{
	static constexpr size_t dim = TDim;
	static constexpr size_t dimSQ = TDimSQ;
	if (dimSQ != dim * dim)
		UG_THROW("MatToTens4::Invalid dimensions: " << dimSQ << " is not the square of " << dim << "! \n");

	for (size_t i = 0; i < dim; ++i)
		for (size_t j = 0; j < dim; ++j)
			for (size_t k = 0; k < dim; ++k)
				for (size_t l = 0; l < dim; ++l){
					tens4[i][j][k][l] = mat(i*dim + j, k*dim + l);
				}
}


////////////////////////////////////////////////////////////////////////////////
// Addition of tensors
////////////////////////////////////////////////////////////////////////////////

//	adds a fourth order tensor:
//	tens4_out = tens4a + tens4b
template <std::size_t TDim>
void
Tens4Add(MathTensor4<TDim, TDim, TDim, TDim>& tens4_out,
		const MathTensor4<TDim, TDim, TDim, TDim>& tens4a,
		const MathTensor4<TDim, TDim, TDim, TDim>& tens4b)
{
	static constexpr size_t dim = TDim;

	for (size_t i = 0; i < dim; ++i)
		for (size_t j = 0; j < dim; ++j)
			for (size_t k = 0; k < dim; ++k)
				for (size_t l = 0; l < dim; ++l){
					tens4_out[i][j][k][l] = tens4a[i][j][k][l]
					                      + tens4b[i][j][k][l];
				}
}


////////////////////////////////////////////////////////////////////////////////
// Subtraction of tensors
////////////////////////////////////////////////////////////////////////////////

//	subtracts a fourth order tensor:
//	tens4_out = tens4a - tens4b
template <std::size_t TDim>
void
Tens4Subtract(MathTensor4<TDim, TDim, TDim, TDim>& tens4_out,
		const MathTensor4<TDim, TDim, TDim, TDim>& tens4a,
		const MathTensor4<TDim, TDim, TDim, TDim>& tens4b)
{
	static constexpr size_t dim = TDim;

	for (size_t i = 0; i < dim; ++i)
		for (size_t j = 0; j < dim; ++j)
			for (size_t k = 0; k < dim; ++k)
				for (size_t l = 0; l < dim; ++l){
					tens4_out[i][j][k][l] = tens4a[i][j][k][l]
					                      - tens4b[i][j][k][l];
				}
}


////////////////////////////////////////////////////////////////////////////////
// Transpose of a tensor
////////////////////////////////////////////////////////////////////////////////


//	transposes a fourth order tensor:
//	tens4_out = (tens4)^T
template <std::size_t TDim>
void
TransTens4(MathTensor4<TDim, TDim, TDim, TDim>& tens4_out,
		const MathTensor4<TDim, TDim, TDim, TDim>& tens4)
{
	static constexpr size_t dim = TDim;

	for (size_t i = 0; i < dim; ++i)
		for (size_t j = 0; j < dim; ++j)
			for (size_t k = 0; k < dim; ++k)
				for (size_t l = 0; l < dim; ++l){
					tens4_out[k][l][i][j] = tens4[i][j][k][l];
				}
}


////////////////////////////////////////////////////////////////////////////////
// Invert a rank-4-tensor
////////////////////////////////////////////////////////////////////////////////

//	This function inverts a tensor of order 4 by inverting the associated matrix
//	tens4_out = (tens4)^{-1}
template <std::size_t TDim>
void
InvertTensor4(MathTensor4<TDim, TDim, TDim, TDim>& tens4_out,
		const MathTensor4<TDim, TDim, TDim, TDim>& tens4)
{
	static constexpr size_t dim = TDim;
	static constexpr size_t dimSQ = dim * dim;

	DenseMatrix< FixedArray2<number, dimSQ, dimSQ> > mat;

//	get associated matrix from tens4
	Tens4ToMat(mat, tens4);

//	we now create a matrix, where we store the inverse matrix
	typename block_traits<DenseMatrix< FixedArray2<number, dimSQ, dimSQ> > >::inverse_type inv;

//	get the inverse
	if(!GetInverse(inv, mat))
		UG_THROW("Could not compute inverse.");

//	get associated tensor of rank 4 from the inverse matrix
	MatToTens4(tens4_out, inv);
}


/*
 * This function solves
 *
 * 	A : X = RHS
 *
 * X, RHS: tensors of second order
 * A: tensor of fourth order
 *
 * The inversion of the 4th order tensors (R^dim x R^dim x R^dim x R^dim) therein is done by associating
 * a R^(dim^2) x R^(dim^2) Matrix to the 4th order tensor and by associating a R^(dim^2)
 * vector to the second order tensors.
 * Hence, we can solve this 'associated' matrix-vector system (A_mat : x_vec = rhs_vec)
 */
template <std::size_t TDim>
void
SolveTensorMatrixEquation(MathMatrix<TDim, TDim>& X,
		const MathTensor4<TDim, TDim, TDim, TDim>& A,
		const MathMatrix<TDim, TDim>& rhs)
{
	static constexpr size_t dim = TDim;
	static constexpr size_t dimSQ = dim * dim;

	DenseMatrix< FixedArray2<number, dimSQ, dimSQ> > A_mat;
	DenseVector< FixedArray1<number, dimSQ> > rhs_vec, x_vec;

//	transform MathTensors of order 4 resp. 2 to DenseMatrices resp. DenseVectors:
	Tens4ToMat(A_mat,A);
	Tens2ToVec(rhs_vec,rhs);

	typename block_traits<DenseMatrix< FixedArray2<number, dimSQ, dimSQ> > >::inverse_type A_matI;

//	get the inverse
	if(!GetInverse(A_matI, A_mat))
		UG_THROW("Could not compute inverse.");

//	x_vec = A_matI * rhs_vec
	MatMult(x_vec, 1.0, A_matI, rhs_vec);
//	get MathMatrix X from DenseVector x:
	VecToTens2(X,x_vec);
}


////////////////////////////////////////////////////////////////////////////////
// Contractions of rank-4-tensors
////////////////////////////////////////////////////////////////////////////////

//	this function computes the contraction of a 4th order tensor by a second order tensor
//	tens2_out = tens4 : tens2
template <std::size_t TDim>
void
Tens4Contract(MathMatrix<TDim, TDim>& tens2_out,
		const MathTensor4<TDim, TDim, TDim, TDim>& tens4,
		const MathMatrix<TDim, TDim>& tens2)
{
	static constexpr size_t dim = TDim;

	for(size_t i = 0; i < dim; ++i)
		for(size_t j = 0; j < dim; ++j)
		{
			tens2_out[i][j] = 0.0;

			for(size_t k = 0; k < dim; ++k)
				for(size_t l = 0; l < dim; ++l)
				{
					tens2_out[i][j] += tens4[i][j][k][l] * tens2[k][l];
				}
		}
}

//	this function computes the contraction of a 4th order tensor by a fourth order tensor
//	tens4_out = tens4a : tens4b
template <std::size_t TDim>
void
Tens4Contract(MathTensor4<TDim, TDim, TDim, TDim>& tens4_out,
		const MathTensor4<TDim, TDim, TDim, TDim>& tens4a,
		const MathTensor4<TDim, TDim, TDim, TDim>& tens4b)
{
	static constexpr size_t dim = TDim;

	for(size_t i = 0; i < dim; ++i)
		for(size_t j = 0; j < dim; ++j)
			for(size_t k = 0; k < dim; ++k)
				for(size_t l = 0; l < dim; ++l)
				{
					tens4_out[i][j][k][l] = 0.0;

					for(size_t m = 0; m < dim; ++m)
						for(size_t n = 0; n < dim; ++n)
						{
							tens4_out[i][j][k][l] += tens4a[i][j][m][n] * tens4b[m][n][k][l];
						}
				}
}

//	this function computes
//	tens4_out = tens4a : tens4b : tens4c
template <std::size_t TDim>
void
Tens4Contract(MathTensor4<TDim, TDim, TDim, TDim>& tens4_out,
		const MathTensor4<TDim, TDim, TDim, TDim>& tens4a,
		const MathTensor4<TDim, TDim, TDim, TDim>& tens4b,
		const MathTensor4<TDim, TDim, TDim, TDim>& tens4c)
{
	static constexpr size_t dim = TDim;

	MathTensor4<dim, dim, dim, dim> help;

	for(size_t i = 0; i < dim; ++i)
		for(size_t j = 0; j < dim; ++j)
			for(size_t k = 0; k < dim; ++k)
				for(size_t l = 0; l < dim; ++l)
				{
					tens4_out[i][j][k][l] = 0.0;

					for(size_t m = 0; m < dim; ++m)
						for(size_t n = 0; n < dim; ++n)
						{
							help[m][n][k][l] = 0.0;

							for(size_t r = 0; r < dim; ++r)
								for(size_t s = 0; s < dim; ++s)
								{
									help[m][n][k][l] += tens4b[m][n][r][s] * tens4c[r][s][k][l];
								}

							tens4_out[i][j][k][l] += tens4a[i][j][m][n] * help[m][n][k][l];
						}
				}
}

////////////////////////////////////////////////////////////////////////////////
// Some tensor func
////////////////////////////////////////////////////////////////////////////////

template <std::size_t TDim>
void
Tens4Zero(MathTensor4<TDim, TDim, TDim, TDim>& tensOut)
{
	static constexpr size_t dim = TDim;

	for(size_t i = 0; i < dim; ++i)
		for(size_t j = 0; j < dim; ++j)
			for(size_t k = 0; k < dim; ++k)
				for(size_t l = 0; l < dim; ++l)
					tensOut[i][j][k][l] = 0.0;

}

//	this function computes the 4th order identity tensor
template <std::size_t TDim>
void
Tens4Identity(MathTensor4<TDim, TDim, TDim, TDim>& Ident)
{
	static constexpr size_t dim = TDim;

	for(size_t i = 0; i < dim; ++i)
		for(size_t j = 0; j < dim; ++j)
			for(size_t k = 0; k < dim; ++k)
				for(size_t l = 0; l < dim; ++l)
				{
					Ident[i][j][k][l] = 0.0;
					if ((i==k) && (j==l)) Ident[i][j][k][l] = 1.0;
				}
}

//	this function computes the symmetric 4th order identity tensor
template <std::size_t TDim>
void
Tens4IdentitySym(MathTensor4<TDim, TDim, TDim, TDim>& Ident)
{
	static constexpr size_t dim = TDim;

	for(size_t i = 0; i < dim; ++i)
		for(size_t j = 0; j < dim; ++j)
			for(size_t k = 0; k < dim; ++k)
				for(size_t l = 0; l < dim; ++l)
				{
					Ident[i][j][k][l] = 0.0;
					if ((i==k) && (j==l)) Ident[i][j][k][l] += 0.5;
					if ((i==l) && (j==k)) Ident[i][j][k][l] += 0.5;
				}
}


} // end of namespace

#endif /* __H__UG__COMMON__MATH_TENSOR_FUNCTIONS_COMMON_IMPL__ */
