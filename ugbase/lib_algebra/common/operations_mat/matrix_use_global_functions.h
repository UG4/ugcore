/*
 * Copyright (c) 2011-2013:  G-CSC, Goethe University Frankfurt
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

// MATRIX_USE_GLOBAL_FUNCTIONS
///////////////////////////////////////////////////////////////////////////////////////////////////

#include "../operations_vec.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
//! calculates dest = alpha1*v1 + beta1 * A1 *w1;
template<typename vector_t, typename matrix_t>
inline bool MatMultAddDirect(vector_t &dest,
		const number &beta1, const matrix_t &A1, const vector_t &w1,
		const number &alpha1, const vector_t &v1)
{
	VecScaleAssign(dest, alpha1, v1);
	return MatMultAddDirect(dest, beta1, A1.cast(), w1);
}

//! calculates dest = alpha1*v1 + beta1 * A1 *w1 + alpha2*v2;
template<typename vector_t, typename matrix_t>
inline bool MatMultAddDirect(vector_t &dest,
		const number &beta1, const matrix_t &A1, const vector_t &w1,
		const number &alpha1, const vector_t &v1,
		const number &alpha2, const vector_t &v2)
{
	bool b = MatMultDirect(dest, beta1, A1, w1);
	VecScaleAdd(dest, 1.0, dest, alpha1, v1, alpha2, v2);
	return b;
}


//! calculates dest = beta1 * A1 *w1 + beta2 * A2*w2;
template<typename vector_t, typename matrix_t>
inline bool MatMultAddDirect(vector_t &dest,
		const number &beta1, const matrix_t &A1, const vector_t &w1,
		const number &beta2, const matrix_t &A2, const vector_t &w2)
{
	MatMultDirect(dest, beta1, A1.cast(), w1);
	return MatMultAddDirect(dest, beta2, A2, w2, 1.0, dest);
}


//! calculates dest = beta1 * A1 *w1 + beta2 * A2*w2 + alpha1*v1;
template<typename vector_t, typename matrix_t>
inline bool MatMultAddDirect(vector_t &dest,
		const number &beta1, const matrix_t &A1, const vector_t &w1,
		const number &beta2, const matrix_t &A2, const vector_t &w2,
		const number &alpha1, const vector_t &v1)
{
	bool b = MatMultAddDirect(dest, beta1, A1, w1, beta2, A2, w2);
	VecScaleAdd(dest, 1.0, dest, alpha1, v1);
	return b;
}


template<typename vector_t, typename matrix_t>
struct mat_operations_class<vector_t, matrix_t, MATRIX_USE_GLOBAL_FUNCTIONS>
{
	//! calculates dest = beta1 * A1;
	static inline bool MatMult(vector_t &dest,
			const number &beta1, const matrix_t &A1, const vector_t &w1)
	{
		return MatMultDirect(dest, beta1, A1, w1);
	}

	//! calculates dest = alpha1*v1 + beta1 * A1 *w1;
	static inline bool MatMultAdd(vector_t &dest,
			const number &alpha1, const vector_t &v1,
			const number &beta1, const matrix_t &A1, const vector_t &w1)
	{
		return MatMultAddDirect(dest, alpha1, v1, beta1, A1, w1);
	}

	//! calculates dest = alpha1*v1 + beta1 * A1 *w1 + alpha2*v2;
	static inline bool MatMultAdd(vector_t &dest,
			const number &alpha1, const vector_t &v1,
			const number &alpha2, const vector_t &v2,
			const number &beta1, const matrix_t &A1, const vector_t &w1)
	{
		return MatMultAddDirect(dest, alpha1, v1, alpha2, v2, beta1, A1, w1);
	}

	//! calculates dest = beta1 * A1 *w1 + beta2 * A2*w2;
	static inline bool MatMultAdd(vector_t &dest,
			const number &beta1, const matrix_t &A1, const vector_t &w1,
			const number &beta2, const matrix_t &A2, const vector_t &w2)
	{
		return MatMultAddDirect(dest, beta1, A1, w1, A2, w2);
	}

	//! calculates dest = beta1 * A1 *w1 + beta2 * A2*w2 + alpha1*v1;
	static inline bool MatMultAdd(vector_t &dest,
			const number &alpha1, const vector_t &v1,
			const number &beta1, const matrix_t &A1, const vector_t &w1,
			const number &beta2, const matrix_t &A2, const vector_t &w2)
	{
		return MatMultAddDirect(dest, alpha1, v1, beta1, A1, w1, beta2, A2, w2);
	}

	//! calculates dest = beta1 * A1^T *w1;
	static inline bool MatMultTransposed(vector_t &dest,
			const number &beta1, const matrix_t &A1, const vector_t &w1)
	{
		return MatMultTransposed(dest, beta1, A1, w1);
	}

	//! calculates dest = alpha1*v1 + beta1 * A1^T *w1;
	static inline bool MatMultTransposedAdd(vector_t &dest,
			const number &alpha1, const vector_t &v1,
			const number &beta1, const matrix_t &A1, const vector_t &w1)
	{
		return MatMultTransposedAddDirect(dest, alpha1, w1, beta1, A1, w1);
	}
};
}
