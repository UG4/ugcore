
// MATRIX_USE_OPERATORS: elementary types like double, float, or template operatored classes: use operators
///////////////////////////////////////////////////////////////////////////////////////////////////

namespace ug{

template<typename vector_t, typename matrix_t>
struct mat_operations_class<vector_t, matrix_t, MATRIX_USE_OPERATORS>
{
	//! calculates dest = beta1 * A1 *w1;
	static inline bool MatMult(vector_t &dest, const number &beta1, const matrix_t &A1, const vector_t &w1)
	{
		dest = beta1 * A1 * w1;
		return true;
	}

	//! calculates dest = alpha1*v1 + beta1 * A1 *w1;
	static inline bool MatMultAdd(vector_t &dest,
			const number &alpha1, const vector_t &v1,
			const number &beta1, const matrix_t &A1, const vector_t &w1)
	{
		dest = alpha1*v1 + beta1 * A1 * w1;
		return true;
	}
	static inline bool MatMultAdd(vector_t &dest,
			const number &alpha1, const vector_t &v1,
			const number &alpha2, const vector_t &v2,
			const number &beta1, const matrix_t &A1, const vector_t &w1)
	{
		dest = alpha1*v1 + alpha2*v2 + beta1 * A1 * w1;
		return true;
	}
	static inline bool MatMultAdd(vector_t &dest,
			const number &beta1, const matrix_t &A1, const vector_t &w1,
			const number &beta2, const matrix_t &A2, const vector_t &w2)
	{
		dest = beta1 * A1 * w1 + beta2 * A2 * w2;
		return true;
	}
	static inline bool MatMultAdd(vector_t &dest,
			const number &alpha1, const vector_t &v1,
			const number &beta1, const matrix_t &A1, const vector_t &w1,
			const number &beta2, const matrix_t &A2, const vector_t &w2)
	{
		dest = alpha1*v1 + beta1 * A1 * w1 + beta2 * A2 * w2;
		return true;
	}

	//! calculates dest = beta1 * A1^T *w1;
	static inline bool MatMultTransposed(vector_t &dest,
			const number &beta1, const matrix_t &A1, const vector_t &w1)
	{
		dest = beta1 * A1 * w1;
		return true;
	}

	//! calculates dest = alpha1*v1 + beta1 * A1^T *w1;
	static inline bool MatMultTransposedAdd(vector_t &dest,
			const number &alpha1, const vector_t &v1,
			const number &beta1, const matrix_t &A1, const vector_t &w1)
	{
		dest = alpha1*v1 + beta1 * A1 * w1;
		return true;
	}
};

}
