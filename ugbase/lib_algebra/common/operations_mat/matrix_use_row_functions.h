#include "../operations_vec.h"
namespace ug
{
// MATRIX_USE_ROW_FUNCTIONS
///////////////////////////////////////////////////////////////////////////////////////////////////

template<typename vector_t, typename matrix_t>
struct mat_operations_class<vector_t, matrix_t, MATRIX_USE_ROW_FUNCTIONS>
{
	//! calculates dest = beta1 * A1;
	static inline bool MatMult(vector_t &dest,
			const number &beta1, const matrix_t &A1, const vector_t &w1)
	{
		for(size_t i=0; i<dest.size(); i++)
		{
			dest[i] = 0.0;
			A1.mat_mult_add_row(i, dest[i], beta1, w1);
		}
		return true;
	}

	//! calculates dest = alpha1*v1 + beta1 * A1 *w1;
	static inline bool MatMultAdd(vector_t &dest,
			const number &alpha1, const vector_t &v1,
			const number &beta1, const matrix_t &A1, const vector_t &w1)
	{
		for(size_t i=0; i<dest.size(); i++)
		{
			VecScaleAssign(dest[i], alpha1, v1[i]);
			A1.mat_mult_add_row(i, dest[i], beta1, w1);
		}
		return true;
	}

	//! calculates dest = alpha1*v1 + alpha2*v2 + beta1 * A1 *w1;
	static inline bool MatMultAdd(vector_t &dest,
			const number &alpha1, const vector_t &v1,
			const number &alpha2, const vector_t &v2,
			const number &beta1, const matrix_t &A1, const vector_t &w1)
	{
		for(size_t i=0; i<dest.size(); i++)
		{
			VecScaleAdd(dest[i], alpha1, v1[i], alpha2, v2[i]);
			A1.cast().mat_mult_add_row(i, dest[i], beta1, w1);
		}
		return true;
	}


	//! calculates dest = beta1 * A1 *w1 + beta2 * A2*w2;
	static inline bool MatMultAdd(vector_t &dest,
			const number &beta1, const matrix_t &A1, const vector_t &w1,
			const number &beta2, const matrix_t &A2, const vector_t &w2)
	{
		for(size_t i=0; i<dest.size(); i++)
		{
			dest[i] = 0.0;
			A1.cast().mat_mult_add_row(i, dest[i], beta1, w1);
			A2.cast().mat_mult_add_row(i, dest[i], beta2, w2);
		}
		return true;
	}


	//! calculates dest = beta1 * A1 *w1 + beta2 * A2*w2 + alpha1*v1;
	static inline bool MatMultAdd(vector_t &dest,
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
		return true;
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
		return MatMultTransposedAddDirect(dest, beta1, A1, w1, alpha1, w1);
	}
};
}
