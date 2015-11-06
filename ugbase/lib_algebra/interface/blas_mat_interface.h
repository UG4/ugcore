/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef BLAS_MAT_INTERFACE_H_
#define BLAS_MAT_INTERFACE_H_



//! calculates dest = beta1 * A1 * w1;
template<typename vector_t, typename matrix_t>
inline bool MatMultScale(vector_t &dest,	const number &beta1, const matrix_t &A1, const vector_t &w1);

//! calculates dest = A1 * w1;
template<typename vector_t, typename matrix_t>
inline bool MatMult(vector_t &dest,	const matrix_t &A1, const vector_t &w1);


//! calculates dest = alpha1*v1 + beta1 * A1 *w1;
template<typename vector_t, typename matrix_t>
inline bool MatMultAdd(vector_t &dest,
		const number &alpha1, const vector_t &v1,
		const number &beta1, const matrix_t &A1, const vector_t &w1	);

//! calculates dest += A1*w1
template<typename vector_t, typename matrix_t>
inline bool MatMultAppend(vector_t &dest,
		const matrix_t &A1, const vector_t &w1	);


//! calculates dest = beta1 * A1^T *w1;
template<typename vector_t, typename matrix_t>
inline bool MatMultScaleTransposed(vector_t &dest,
		const number &beta1, const matrix_t &A1, const vector_t &w1	);

//! calculates dest = A1^T *w1;
template<typename vector_t, typename matrix_t>
inline bool MatMultScaleTransposed(vector_t &dest,
		const matrix_t &A1, const vector_t &w1	);


//! calculates dest = alpha1*v1 + beta1 * A1^T *w1;
template<typename vector_t, typename matrix_t>
inline bool MatMultTransposedScaledAdd(vector_t &dest,
		const number &alpha1, const vector_t &v1,
		const number &beta1, const matrix_t &A1, const vector_t &w1	)
{
	return mat_operations_class<vector_t, matrix_t, matrix_algebra_type_traits<matrix_t>::type>
			::MatMultTransposedAdd(dest, alpha1, v1, beta1, A1, w1);
}

// row functions

//! calculates dest += beta1 * A1[row] *w1;
template<typename vector_t, typename matrix_t>
inline bool MatMultAddRow(typename vector_t::value_type &dest,
		const number &alpha1, const vector_t &v1,
		const number &beta1, const matrix_t &A1, const vector_t &w1	)
{
	return mat_operations_class<vector_t, matrix_t, matrix_algebra_type_traits<matrix_t>::type>
			::MatMultAdd(dest, alpha1, v1, beta1, A1, w1);
}

#endif /* BLAS_MAT_INTERFACE_H_ */
