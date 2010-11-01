/*
 * operations_mat.h
 *
 *  Created on: 29.09.2010
 *      Author: mrupp
 */

#ifndef __H__UG__CPU_ALGEBRA__OPERATIONS_MAT__
#define __H__UG__CPU_ALGEBRA__OPERATIONS_MAT__

namespace ug
{

// operations for matrices
//-----------------------------------------------------------------------------
// these functions execute matrix operations


/*
suppose you have your algebra with your vector class yourvector and your matrix class yourmatrix.
when your matrix class is rowwise multiplicable, and you want to do so, set
template<>
class matrix_algebra_type_traits<yourclass>
{
	static const matrix_algebra_type type = MATRIX_USE_ROW_FUNCTIONS;
};
and add a function
inline void yourmatrix::mat_mult_add_row(size_t row, yourvector::value_type &dest, double beta1, const yourvector &v) const;
to your class.

If you cannot or dont want to use rowwise multiplication, set
template<>
class matrix_algebra_type_traits<yourclass>
{
	static const matrix_algebra_type type = MATRIX_USE_GLOBAL_FUNCTIONS;
};

and implement (at least)
inline void MatMultAddDirect(yourvector &dest, const number &beta1, const yourmatrix &A1, const yourvector &w1);
inline void MatMultDirect(yourvector &dest, const number &beta1, const yourmatrix &A1, const yourvector &w1,
							number &alpha1, const yourvector &v1);
Other functions like MatMultAdd(dest, beta1, A1, w1, beta2, A2, w2, alpha1, v1) are then implemented through your functions.
(see class mat_operations_class<vector_t, matrix_t, MATRIX_USE_GLOBAL_FUNCTIONS> for all functions).

*/


enum matrix_algebra_type
{
	MATRIX_USE_ROW_FUNCTIONS,
	MATRIX_USE_GLOBAL_FUNCTIONS,
	MATRIX_USE_OPERATORS,
	MATRIX_USE_MEMBER_FUNCTIONS
};

template<typename vector_t, typename matrix_t, matrix_algebra_type type>
class mat_operations_class;

template<typename T>
struct matrix_algebra_type_traits
{
	static const matrix_algebra_type type = MATRIX_USE_OPERATORS;
};

// MATRIX_USE_OPERATORS: elementary types like double, float, or template operatored classes: use operators
///////////////////////////////////////////////////////////////////////////////////////////////////

template<typename vector_t, typename matrix_t>
struct mat_operations_class<vector_t, matrix_t, MATRIX_USE_OPERATORS>
{
	//! calculates dest = beta1 * A1 *w1;
	static inline void MatMult(vector_t &dest, const number &beta1, const matrix_t &A1, const vector_t &w1)
	{
		dest = beta1 * A1 * w1;
	}

	//! calculates dest = alpha1*v1 + beta1 * A1 *w1;
	static inline void MatMultAdd(vector_t &dest,
			const number &alpha1, const vector_t &v1,
			const number &beta1, const matrix_t &A1, const vector_t &w1)
	{
		dest = alpha1*v1 + beta1 * A1 * w1;
	}
	static inline void MatMultAdd(vector_t &dest,
			const number &alpha1, const vector_t &v1,
			const number &alpha2, const vector_t &v2,
			const number &beta1, const matrix_t &A1, const vector_t &w1)
	{
		dest = alpha1*v1 + alpha2*v2 + beta1 * A1 * w1;
	}
	static inline void MatMultAdd(vector_t &dest,
			const number &beta1, const matrix_t &A1, const vector_t &w1,
			const number &beta2, const matrix_t &A2, const vector_t &w2)
	{
		dest = beta1 * A1 * w1 + beta2 * A2 * w2;
	}
	static inline void MatMultAdd(vector_t &dest,
			const number &alpha1, const vector_t &v1,
			const number &beta1, const matrix_t &A1, const vector_t &w1,
			const number &beta2, const matrix_t &A2, const vector_t &w2)
	{
		dest = alpha1*v1 + beta1 * A1 * w1 + beta2 * A2 * w2;
	}

	//! calculates dest = alpha1*v1 + beta1 * A1 *w1;
	static inline void MatMultTransposedAdd(vector_t &dest,
			const number &alpha1, const vector_t &v1,
			const number &beta1, const matrix_t &A1, const vector_t &w1)
	{
		dest = alpha1*v1 + beta1 * A1 * w1;
	}
};

// MATRIX_USE_GLOBAL_FUNCTIONS
///////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//! calculates dest = alpha1*v1 + beta1 * A1 *w1;
template<typename vector_t, typename matrix_t>
inline void MatMultAddDirect(vector_t &dest,
		const number &beta1, const matrix_t &A1, const vector_t &w1,
		const number &alpha1, const vector_t &v1)
{
	VecScaleAssign(dest, alpha1, v1);
	MatMultAddDirect(dest, beta1, A1.cast(), w1);
}

//! calculates dest = alpha1*v1 + beta1 * A1 *w1 + alpha2*v2;
template<typename vector_t, typename matrix_t>
inline void MatMultAddDirect(vector_t &dest,
		const number &beta1, const matrix_t &A1, const vector_t &w1,
		const number &alpha1, const vector_t &v1,
		const number &alpha2, const vector_t &v2)
{
	MatMultDirect(dest, beta1, A1, w1);
	VecScaleAdd(dest, 1.0, dest, alpha1, v1, alpha2, v2);
}


//! calculates dest = beta1 * A1 *w1 + beta2 * A2*w2;
template<typename vector_t, typename matrix_t>
inline void MatMultAddDirect(vector_t &dest,
		const number &beta1, const matrix_t &A1, const vector_t &w1,
		const number &beta2, const matrix_t &A2, const vector_t &w2)
{
	MatMultDirect(dest, beta1, A1.cast(), w1);
	MatMultAddDirect(dest, beta2, A2, w2, 1.0, dest);
}


//! calculates dest = beta1 * A1 *w1 + beta2 * A2*w2 + alpha1*v1;
template<typename vector_t, typename matrix_t>
inline void MatMultAddDirect(vector_t &dest,
		const number &beta1, const matrix_t &A1, const vector_t &w1,
		const number &beta2, const matrix_t &A2, const vector_t &w2,
		const number &alpha1, const vector_t &v1)
{
	MatMultAddDirect(dest, beta1, A1, w1, beta2, A2, w2);
	VecScaleAdd(dest, 1.0, dest, alpha1, v1);
}


template<typename vector_t, typename matrix_t>
struct mat_operations_class<vector_t, matrix_t, MATRIX_USE_GLOBAL_FUNCTIONS>
{
	//! calculates dest = beta1 * A1;
	static inline void MatMult(vector_t &dest,
			const number &beta1, const matrix_t &A1, const vector_t &w1)
	{
		MatMultDirect(dest, beta1, A1, w1);
	}

	//! calculates dest = alpha1*v1 + beta1 * A1 *w1;
	static inline void MatMultAdd(vector_t &dest,
			const number &beta1, const matrix_t &A1, const vector_t &w1,
			const number &alpha1, const vector_t &v1)
	{
		MatMultAddDirect(dest, beta1, A1, w1, alpha1, w1);
	}

	//! calculates dest = alpha1*v1 + beta1 * A1 *w1 + alpha2*v2;
	static inline void MatMultAdd(vector_t &dest,
			const number &beta1, const matrix_t &A1, const vector_t &w1,
			const number &alpha1, const vector_t &v1,
			const number &alpha2, const vector_t &v2)
	{
		MatMultAddDirect(dest, beta1, A1, w1, alpha1, w1, alpha2, v2);
	}

	//! calculates dest = beta1 * A1 *w1 + beta2 * A2*w2;
	static inline void MatMultAdd(vector_t &dest,
			const number &beta1, const matrix_t &A1, const vector_t &w1,
			const number &beta2, const matrix_t &A2, const vector_t &w2)
	{
		MatMultAddDirect(dest, beta1, A1, w1, A2, w2);
	}

	//! calculates dest = beta1 * A1 *w1 + beta2 * A2*w2 + alpha1*v1;
	static inline void MatMultAdd(vector_t &dest,
			const number &beta1, const matrix_t &A1, const vector_t &w1,
			const number &beta2, const matrix_t &A2, const vector_t &w2,
			const number &alpha1, const vector_t &v1)
	{
		MatMultAddDirect(dest, beta1, A1, w1, beta2, A2, w2, alpha1, v1);
	}

	//! calculates dest = alpha1*v1 + beta1 * A1^T *w1;
	static inline void MatMultTransposedAdd(vector_t &dest,
			const number &alpha1, const vector_t &v1,
			const number &beta1, const matrix_t &A1, const vector_t &w1)
	{
		MatMultTransposedAddDirect(dest, beta1, A1, w1, alpha1, w1);
	}
};


// MATRIX_USE_ROW_FUNCTIONS
///////////////////////////////////////////////////////////////////////////////////////////////////

template<typename vector_t, typename matrix_t>
struct mat_operations_class<vector_t, matrix_t, MATRIX_USE_ROW_FUNCTIONS>
{
	//! calculates dest = beta1 * A1;
	static inline void MatMult(vector_t &dest,
			const number &beta1, const matrix_t &A1, const vector_t &w1)
	{
		for(size_t i=0; i<dest.size(); i++)
		{
			dest[i] = 0.0;
			A1.mat_mult_add_row(i, dest[i], beta1, w1);
		}
	}

	//! calculates dest = alpha1*v1 + beta1 * A1 *w1;
	static inline void MatMultAdd(vector_t &dest,
			const number &alpha1, const vector_t &v1,
			const number &beta1, const matrix_t &A1, const vector_t &w1)
	{
		for(size_t i=0; i<dest.size(); i++)
		{
			VecScaleAssign(dest[i], alpha1, v1[i]);
			A1.mat_mult_add_row(i, dest[i], beta1, w1);
		}
	}

	//! calculates dest = alpha1*v1 + alpha2*v2 + beta1 * A1 *w1;
	static inline void MatMultAdd(vector_t &dest,
			const number &alpha1, const vector_t &v1,
			const number &alpha2, const vector_t &v2,
			const number &beta1, const matrix_t &A1, const vector_t &w1)
	{
		for(size_t i=0; i<dest.size(); i++)
		{
			VecScaleAdd(dest[i], alpha1, v1[i], alpha2, v2[i]);
			A1.cast().mat_mult_add_row(i, dest[i], beta1, w1);
		}
	}


	//! calculates dest = beta1 * A1 *w1 + beta2 * A2*w2;
	static inline void MatMultAdd(vector_t &dest,
			const number &beta1, const matrix_t &A1, const vector_t &w1,
			const number &beta2, const matrix_t &A2, const vector_t &w2)
	{
		for(size_t i=0; i<dest.size(); i++)
		{
			dest[i] = 0.0;
			A1.cast().mat_mult_add_row(i, dest[i], beta1, w1);
			A2.cast().mat_mult_add_row(i, dest[i], beta2, w2);
		}
	}


	//! calculates dest = beta1 * A1 *w1 + beta2 * A2*w2 + alpha1*v1;
	static inline void MatMultAdd(vector_t &dest,
			const number &alpha1, const vector_t &v1,
			const number &beta1, const matrix_t &A1, const vector_t &w1,
			const number &beta2, const matrix_t &A2, const vector_t &w2)
	{
		for(size_t i=0; i<dest.size(); i++)
		{
			VecScaleAssign(dest[i], alpha1, v1[i]);
			A1.cast().mat_mult_add_row(i, dest[i], beta1, w1);
			A2.cast().mat_mult_add_row(i, dest[i], beta2, w2);
		}
	}

	//! calculates dest = alpha1*v1 + beta1 * A1^T *w1;
	static inline void MatMultTransposedAdd(vector_t &dest,
			const number &alpha1, const vector_t &v1,
			const number &beta1, const matrix_t &A1, const vector_t &w1)
	{
		MatMultTransposedAddDirect(dest, beta1, A1, w1, alpha1, w1);
	}
};

// MATRIX_USE_MEMBER_FUNCTIONS
///////////////////////////////////////////////////////////////////////////////////////////////////

template<typename vector_t, typename matrix_t>
struct mat_operations_class<vector_t, matrix_t, MATRIX_USE_MEMBER_FUNCTIONS>
{
	//! calculates dest = beta1 * A1;
	static inline void MatMult(vector_t &dest,
			const number &beta1, const matrix_t &A1, const vector_t &w1)
	{
		dest.mat_mult(beta1, A1, w1);
	}

	//! calculates dest = alpha1*v1 + beta1 * A1 *w1;
	static inline void MatMultAdd(vector_t &dest,
			const number &alpha1, const vector_t &v1,
			const number &beta1, const matrix_t &A1, const vector_t &w1)
	{
		dest.mat_mult_add(alpha1, w1, beta1, A1, w1);
	}

	//! calculates dest = alpha1*v1 + beta1 * A1 *w1 + alpha2*v2;
	static inline void MatMultAdd(vector_t &dest,
			const number &alpha1, const vector_t &v1,
			const number &alpha2, const vector_t &v2,
			const number &beta1, const matrix_t &A1, const vector_t &w1)
	{
		dest.mat_mult_add(alpha1, w1, alpha2, v2, beta1, A1, w1);
	}

	//! calculates dest = beta1 * A1 *w1 + beta2 * A2*w2;
	static inline void MatMultAdd(vector_t &dest,
			const number &beta1, const matrix_t &A1, const vector_t &w1,
			const number &beta2, const matrix_t &A2, const vector_t &w2)
	{
		dest.mat_mult_add(beta1, A1, w1, A2, w2);
	}

	//! calculates dest = alpha1*v1 + beta1 * A1 *w1 + beta2 * A2*w2;
	static inline void MatMultAdd(vector_t &dest,
			const number &alpha1, const vector_t &v1,
			const number &beta1, const matrix_t &A1, const vector_t &w1,
			const number &beta2, const matrix_t &A2, const vector_t &w2)
	{
		dest.mat_mult_add(alpha1, v1, beta1, A1, w1, beta2, A2, w2);
	}


	//! calculates dest = alpha1*v1 + beta1 * A1^T *w1;
	static inline void MatMultTransposedAdd(vector_t &dest,
			const number &alpha1, const vector_t &v1,
			const number &beta1, const matrix_t &A1, const vector_t &w1)
	{
		dest.mat_mult_tranposed_add(dest, beta1, A1, w1, alpha1, w1);
	}
};

// functions transforming MatMult into the right operation, depending on matrix_algebra_type_traits.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//! calculates dest = beta1 * A1;
template<typename vector_t, typename matrix_t>
inline void MatMult(vector_t &dest,
		const number &beta1, const matrix_t &A1, const vector_t &w1)
{
	mat_operations_class<vector_t, matrix_t, matrix_algebra_type_traits<matrix_t>::type>
		::MatMult(dest, beta1, A1, w1);
}

//! calculates dest = alpha1*v1 + beta1 * A1 *w1;
template<typename vector_t, typename matrix_t>
inline void MatMultAdd(vector_t &dest,
		const number &alpha1, const vector_t &v1,
		const number &beta1, const matrix_t &A1, const vector_t &w1	)
{
	mat_operations_class<vector_t, matrix_t, matrix_algebra_type_traits<matrix_t>::type>
			::MatMultAdd(dest, alpha1, v1, beta1, A1, w1);
}

//! calculates dest = alpha1*v1 + beta1 * A1 *w1 + alpha2*v2;
template<typename vector_t, typename matrix_t>
inline void MatMultAdd(vector_t &dest,
		const number &alpha1, const vector_t &v1,
		const number &alpha2, const vector_t &v2,
		const number &beta1, const matrix_t &A1, const vector_t &w1)
{
	mat_operations_class<vector_t, matrix_t, matrix_algebra_type_traits<matrix_t>::type>
			::MatMultAdd(dest, alpha1, v1, alpha2, v2, beta1, A1, w1);
}

//! calculates dest = beta1 * A1 *w1 + beta2 * A2*w2;
template<typename vector_t, typename matrix_t>
inline void MatMultAdd(vector_t &dest,
		const number &beta1, const matrix_t &A1, const vector_t &w1,
		const number &beta2, const matrix_t &A2, const vector_t &w2)
{
	mat_operations_class<vector_t, matrix_t, matrix_algebra_type_traits<matrix_t>::type>
			::MatMultAdd(dest, beta1, A1, w1, A2, w2);
}

//! calculates dest = beta1 * A1 *w1 + beta2 * A2*w2 + alpha1*v1;
template<typename vector_t, typename matrix_t>
inline void MatMultAdd(vector_t &dest,
		const number &alpha1, const vector_t &v1,
		const number &beta1, const matrix_t &A1, const vector_t &w1,
		const number &beta2, const matrix_t &A2, const vector_t &w2)
{
	mat_operations_class<vector_t, matrix_t, matrix_algebra_type_traits<matrix_t>::type>
			::MatMultAdd(dest, alpha1, v1, beta1, A1, w1, beta2, A2, w2);
}



//! calculates dest = alpha1*v1 + beta1 * A1 *w1;
template<typename vector_t, typename matrix_t>
inline void MatMultTransposedAdd(vector_t &dest,
		const number &alpha1, const vector_t &v1,
		const number &beta1, const matrix_t &A1, const vector_t &w1	)
{
	mat_operations_class<vector_t, matrix_t, matrix_algebra_type_traits<matrix_t>::type>
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
#endif /* __H__UG__CPU_ALGEBRA__OPERATIONS_MAT__ */
