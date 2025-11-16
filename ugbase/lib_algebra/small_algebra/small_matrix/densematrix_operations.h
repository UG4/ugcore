/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__COMMON__DENSEMATRIX_OPERATIONS_H__
#define __H__UG__COMMON__DENSEMATRIX_OPERATIONS_H__

#include "densematrix.h"
#include "densevector.h"

#include "../../common/operations.h"

namespace ug{

/// \addtogroup small_algebra
/// \{

template<typename T>
struct matrix_algebra_type_traits;

template<typename T>
struct matrix_algebra_type_traits<DenseMatrix<T> >
{
	static constexpr int type = MATRIX_USE_GLOBAL_FUNCTIONS;
};

//! calculates dest = beta1 * A1 * w1;
template<typename vector_t, typename matrix_t>
inline void MatMult(DenseVector<vector_t> &dest,
		const number &beta1, const DenseMatrix<matrix_t> &A1, const DenseVector<vector_t> &w1)
{
	for(size_t r = 0; r < dest.size(); ++r)
	{
		MatMult(dest[r], beta1, A1(r,0), w1[0]);
		for(size_t c = 1; c < w1.size(); ++c)
			MatMultAdd(dest[r], 1.0, dest[r], beta1, A1(r,c), w1[c]);
	}
}

//! calculates dest = alpha1*v1 + beta1 * A1 *w1;
template<typename vector_t, typename matrix_t>
inline void MatMultAdd(DenseVector<vector_t> &dest,
		const number &alpha1, const DenseVector<vector_t> &v1,
		const number &beta1, const DenseMatrix<matrix_t> &A1, const DenseVector<vector_t> &w1)
{
	for(size_t r = 0; r < dest.size(); ++r)
	{
		VecScaleAssign(dest[r], alpha1, v1[r]);
		for(size_t c = 0; c < w1.size(); ++c)
			MatMultAdd(dest[r], 1.0, dest[r], beta1, A1(r,c), w1[c]);
	}
}

//! calculates dest = alpha1*v1 + beta1 * A1^T *w1;
template<typename vector_t, typename matrix_t>
inline void MatMultTransposedAdd(DenseVector<vector_t> &dest,
		const number &alpha1, const DenseVector<vector_t> &v1,
		const number &beta1, const DenseMatrix<matrix_t> &A1, const DenseVector<vector_t> &w1)
{
	for(size_t r = 0; r < dest.size(); ++r)
	{
		VecScaleAssign(dest[r], alpha1, v1[r]);
		for(size_t c = 0; c < w1.size(); ++c)
			MatMultTransposedAdd(dest[r], 1.0, dest[r], beta1, A1(c,r), w1[c]);
	}
}

// end group small_algebra
/// \}

}

#endif // __H__UG__COMMON__DENSEMATRIX_OPERATIONS_H__
