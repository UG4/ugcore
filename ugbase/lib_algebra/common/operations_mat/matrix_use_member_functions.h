
// MATRIX_USE_MEMBER_FUNCTIONS
///////////////////////////////////////////////////////////////////////////////////////////////////

namespace ug
{

template<typename vector_t, typename matrix_t>
struct mat_operations_class<vector_t, matrix_t, MATRIX_USE_MEMBER_FUNCTIONS>
{
	//! calculates dest = beta1 * A1;
	static inline bool MatMult(vector_t &dest,
			const number &beta1, const matrix_t &A1, const vector_t &w1)
	{
		return dest.mat_mult(beta1, A1, w1);
	}

	//! calculates dest = alpha1*v1 + beta1 * A1 *w1;
	static inline bool MatMultAdd(vector_t &dest,
			const number &alpha1, const vector_t &v1,
			const number &beta1, const matrix_t &A1, const vector_t &w1)
	{
		return dest.mat_mult_add(alpha1, w1, beta1, A1, w1);
	}

	//! calculates dest = alpha1*v1 + beta1 * A1 *w1 + alpha2*v2;
	static inline bool MatMultAdd(vector_t &dest,
			const number &alpha1, const vector_t &v1,
			const number &alpha2, const vector_t &v2,
			const number &beta1, const matrix_t &A1, const vector_t &w1)
	{
		return dest.mat_mult_add(alpha1, w1, alpha2, v2, beta1, A1, w1);
	}

	//! calculates dest = beta1 * A1 *w1 + beta2 * A2*w2;
	static inline bool MatMultAdd(vector_t &dest,
			const number &beta1, const matrix_t &A1, const vector_t &w1,
			const number &beta2, const matrix_t &A2, const vector_t &w2)
	{
		return dest.mat_mult_add(beta1, A1, w1, A2, w2);
	}

	//! calculates dest = alpha1*v1 + beta1 * A1 *w1 + beta2 * A2*w2;
	static inline bool MatMultAdd(vector_t &dest,
			const number &alpha1, const vector_t &v1,
			const number &beta1, const matrix_t &A1, const vector_t &w1,
			const number &beta2, const matrix_t &A2, const vector_t &w2)
	{
		return dest.mat_mult_add(alpha1, v1, beta1, A1, w1, beta2, A2, w2);
	}

	//! calculates dest = beta1 * A1^T *w1;
	static inline bool MatMultTransposed(vector_t &dest,
			const number &beta1, const matrix_t &A1, const vector_t &w1)
	{
		return dest.mat_mult_transposed(beta1, A1, w1);
	}

	//! calculates dest = alpha1*v1 + beta1 * A1^T *w1;
	static inline bool MatMultTransposedAdd(vector_t &dest,
			const number &alpha1, const vector_t &v1,
			const number &beta1, const matrix_t &A1, const vector_t &w1)
	{
		return dest.mat_mult_tranposed_add(alpha1, w1, beta1, A1, w1);
	}
};

}
