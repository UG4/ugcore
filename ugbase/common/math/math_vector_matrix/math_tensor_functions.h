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

#ifndef __H__UG__COMMON__MATH_TENSOR_FUNCTIONS__
#define __H__UG__COMMON__MATH_TENSOR_FUNCTIONS__

#include "math_tensor.h"

#include "lib_algebra/small_algebra/small_matrix/densevector.h"
#include "lib_algebra/small_algebra/small_matrix/densematrix.h"

namespace ug{

/// \addtogroup math_tensor
/// \{

/// transformation of a tensor of rank 2 (R^dim x R^dim) into a vector (R^(dim^2))
template <std::size_t TDim, std::size_t TDimSQ>
void
Tens2ToVec(DenseVector<FixedArray1<number, TDimSQ> > &vec, const MathMatrix<TDim, TDim>& tens2);

///	transformation of a vector (R^(dim^2)) into a tensor of rank 2 (R^dim x R^dim)
template <std::size_t TDim, std::size_t TDimSQ>
void
VecToTens2(MathMatrix<TDim, TDim>& tens2, const DenseVector<FixedArray1<number, TDimSQ> > &vec);

///	transformation of a tensor of rank 4 (R^dim x R^dim x R^dim x R^dim)
///	into a matrix (R^(dim^2) x R^(dim^2))
template <std::size_t TDim, std::size_t TDimSQ>
void
Tens4ToMat(DenseMatrix<FixedArray2<number, TDimSQ, TDimSQ> > &mat, const MathTensor4<TDim, TDim, TDim, TDim>& tens4);

/// transformation of a matrix (R^(dim^2) x R^(dim^2))
///	into a tensor of rank 4 (R^dim x R^dim x R^dim x R^dim)
template <std::size_t TDim, std::size_t TDimSQ>
void
MatToTens4(MathTensor4<TDim, TDim, TDim, TDim>& tens4, const DenseMatrixInverse<FixedArray2<number, TDimSQ, TDimSQ> > &mat);



///	adds a fourth order tensor:
///	tens4_out = tens4a + tens4b
template <std::size_t TDim>
void
Tens4Add(MathTensor4<TDim, TDim, TDim, TDim>& tens4_out,
		const MathTensor4<TDim, TDim, TDim, TDim>& tens4a,
		const MathTensor4<TDim, TDim, TDim, TDim>& tens4b);

///	subtracts a fourth order tensor:
///	tens4_out = tens4a - tens4b
template <std::size_t TDim>
void
Tens4Subtract(MathTensor4<TDim, TDim, TDim, TDim>& tens4_out,
		const MathTensor4<TDim, TDim, TDim, TDim>& tens4a,
		const MathTensor4<TDim, TDim, TDim, TDim>& tens4b);



///	transposes a fourth order tensor:
///	tens4_out = (tens4)^T
template <std::size_t TDim>
void
TransTens4(MathTensor4<TDim, TDim, TDim, TDim>& tens4_out,
		const MathTensor4<TDim, TDim, TDim, TDim>& tens4);



///	This function inverts a tensor of order 4 by inverting the associated matrix
///	tens4_out = (tens4)^{-1}
template <std::size_t TDim>
void
InvertTensor4(MathTensor4<TDim, TDim, TDim, TDim>& tens4_out,
		const MathTensor4<TDim, TDim, TDim, TDim>& tens4);

/// This function solves
///	A : X = RHS
///	X, RHS: tensors of second order
///	A: tensor of fourth order
template <std::size_t TDim>
void
SolveTensorMatrixEquation(MathMatrix<TDim, TDim>& X,
		const MathTensor4<TDim, TDim, TDim, TDim>& A,
		const MathMatrix<TDim, TDim>& rhs);



///	this function computes the contraction of a 4th order tensor by a second order tensor
///	tens2_out = tens4 : tens2
template <std::size_t TDim>
void
Tens4Contract(MathMatrix<TDim, TDim>& tens2_out,
		const MathTensor4<TDim, TDim, TDim, TDim>& tens4,
		const MathMatrix<TDim, TDim>& tens2);

///	this function computes the contraction of a 4th order tensor by a fourth order tensor
///	tens4_out = tens4a : tens4b
template <std::size_t TDim>
void
Tens4Contract(MathTensor4<TDim, TDim, TDim, TDim>& tens4_out,
		const MathTensor4<TDim, TDim, TDim, TDim>& tens4a,
		const MathTensor4<TDim, TDim, TDim, TDim>& tens4b);

///	this function computes
///	tens4_out = tens4a : tens4b : tens4c
template <std::size_t TDim>
void
Tens4Contract(MathTensor4<TDim, TDim, TDim, TDim>& tens4_out,
		const MathTensor4<TDim, TDim, TDim, TDim>& tens4a,
		const MathTensor4<TDim, TDim, TDim, TDim>& tens4b,
		const MathTensor4<TDim, TDim, TDim, TDim>& tens4c);



template <std::size_t TDim>
void
Tens4Zero(MathTensor4<TDim, TDim, TDim, TDim>& tensOut);

///	this function computes the 4th order identity tensor
template <std::size_t TDim>
void
Tens4Identity(MathTensor4<TDim, TDim, TDim, TDim>& Ident);

///	this function computes the symmetric 4th order identity tensor
template <std::size_t TDim>
void
Tens4IdentitySym(MathTensor4<TDim, TDim, TDim, TDim>& Ident);


// end group math_tensor
/// \}

} //end of namespace

//	include a general, but not very fast implementation of the declared methods above.
#include "math_tensor_functions_common_impl.hpp"

#endif /* __H__UG__COMMON__MATH_TENSOR_FUNCTIONS__ */
