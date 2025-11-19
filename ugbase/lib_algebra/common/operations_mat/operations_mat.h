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

#ifndef __H__UG__CPU_ALGEBRA__OPERATIONS_MAT__
#define __H__UG__CPU_ALGEBRA__OPERATIONS_MAT__

#include "common/common.h"
#include "matrix_algebra_types.h"

// operations for matrices
//-----------------------------------------------------------------------------
// these functions execute matrix operations

#include "matrix_use_operators.h"

#include "matrix_use_global_functions.h"

#include "matrix_use_row_functions.h"

#include "matrix_use_member_functions.h"

namespace ug
{

// functions transforming MatMult into the right operation, depending on matrix_algebra_type_traits.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//! calculates dest = beta1 * A1;
template<typename vector_t, typename matrix_t>
inline bool MatMult(vector_t &dest,
		const number &beta1, const matrix_t &A1, const vector_t &w1)
{
	return mat_operations_class<vector_t, matrix_t, matrix_algebra_type_traits<matrix_t>::type>
		::MatMult(dest, beta1, A1, w1);
}

//! calculates dest = alpha1*v1 + beta1 * A1 *w1;
template<typename vector_t, typename matrix_t>
inline bool MatMultAdd(vector_t &dest,
		const number &alpha1, const vector_t &v1,
		const number &beta1, const matrix_t &A1, const vector_t &w1	)
{
	return mat_operations_class<vector_t, matrix_t, matrix_algebra_type_traits<matrix_t>::type>
			::MatMultAdd(dest, alpha1, v1, beta1, A1, w1);
}

//! calculates dest = alpha1*v1 + alpha2*v2 + beta1 * A1 *w1;
template<typename vector_t, typename matrix_t>
inline bool MatMultAdd(vector_t &dest,
		const number &alpha1, const vector_t &v1,
		const number &alpha2, const vector_t &v2,
		const number &beta1, const matrix_t &A1, const vector_t &w1)
{
	return mat_operations_class<vector_t, matrix_t, matrix_algebra_type_traits<matrix_t>::type>
			::MatMultAdd(dest, alpha1, v1, alpha2, v2, beta1, A1, w1);
}

//! calculates dest = beta1 * A1 *w1 + beta2 * A2*w2;
template<typename vector_t, typename matrix_t>
inline bool MatMultAdd(vector_t &dest,
		const number &beta1, const matrix_t &A1, const vector_t &w1,
		const number &beta2, const matrix_t &A2, const vector_t &w2)
{
	return mat_operations_class<vector_t, matrix_t, matrix_algebra_type_traits<matrix_t>::type>
			::MatMultAdd(dest, beta1, A1, w1, A2, w2);
}

//! calculates dest = alpha1*v1 + beta1 * A1 *w1 + beta2 * A2*w2;
template<typename vector_t, typename matrix_t>
inline bool MatMultAdd(vector_t &dest,
		const number &alpha1, const vector_t &v1,
		const number &beta1, const matrix_t &A1, const vector_t &w1,
		const number &beta2, const matrix_t &A2, const vector_t &w2)
{
	return mat_operations_class<vector_t, matrix_t, matrix_algebra_type_traits<matrix_t>::type>
			::MatMultAdd(dest, alpha1, v1, beta1, A1, w1, beta2, A2, w2);
}


//! calculates dest = alpha1*v1 + beta1 * A1 *w1;
template<typename vector_t, typename matrix_t>
inline bool MatMultTransposed(vector_t &dest,
		const number &beta1, const matrix_t &A1, const vector_t &w1	)
{
	return mat_operations_class<vector_t, matrix_t, matrix_algebra_type_traits<matrix_t>::type>
			::MatMultTransposed(dest, beta1, A1, w1);
}


//! calculates dest = alpha1*v1 + beta1 * A1 *w1;
template<typename vector_t, typename matrix_t>
inline bool MatMultTransposedAdd(vector_t &dest,
		const number &alpha1, const vector_t &v1,
		const number &beta1, const matrix_t &A1, const vector_t &w1	)
{
	return mat_operations_class<vector_t, matrix_t, matrix_algebra_type_traits<matrix_t>::type>
			::MatMultTransposedAdd(dest, alpha1, v1, beta1, A1, w1);
}
/*
//! calculates dest = beta1*A1*w1 + alpha1*v1, and norm = norm^2(dest)
template<typename vector_t, typename matrix_t>
inline void MatMultAddNorm(vector_t &dest,
		const number &beta1, const TE_MAT_RowApplicable<matrix_t> &A1, const vector_t &w1,
		const number &alpha1, const vector_t &v1,
		double &norm)
{
	if(OPERATIONS_VERBOSE) cout << "dest: " << dest << endl <<
			"\t= " << beta1 << " * " << A1.cast() << "\t * " << w1 << endl <<
			"\t+ " << alpha1 << " * " << v1 << endl;

	norm=0;
	for(size_t i=0; i<dest.size(); i++)
	{
		VecScaleAssign(dest[i], alpha1, v1[i]);
		A1.cast().mat_mult_add_row(i, dest[i], beta1, w1);
		VecNormSquaredAdd(dest[i], norm);
	}
}
*/



} // namespace ug
#endif