
#ifndef __H__UG__CPU_ALGEBRA__CORE_SMOOTHERS__
#define __H__UG__CPU_ALGEBRA__CORE_SMOOTHERS__
////////////////////////////////////////////////////////////////////////////////////////////////

namespace ug
{

/// \addtogroup lib_algebra
///	@{

/////////////////////////////////////////////////////////////////////////////////////////////
//	gs_step_LL
/**
 * \brief Performs a forward gauss-seidel-step, that is, solve on the lower left of A.
 * \param A Matrix \f$A = D - L - U\f$
 * \param x will be \f$x = (D-L)^{-1}b \f$
 * \param b the vector b.
 * \sa gs_step_UR, sgs_step
 */
template<typename Matrix_type, typename Vector_type>
bool gs_step_LL(const Matrix_type &A, Vector_type &x, const Vector_type &b)
{
	// gs LL has preconditioning matrix
	// (D-L)^{-1}

	typename Vector_type::value_type s;

	for(size_t i=0; i < x.size(); i++)
	{
		s = b[i];

		for(typename Matrix_type::const_row_iterator it = A.begin_row(i); it != A.end_row(i) && it.index() < i; ++it)
			// s -= it.value() * x[it.index()];
			MatMultAdd(s, 1.0, s, -1.0, it.value(), x[it.index()]);

		// x[i] = s/A(i,i)
		InverseMatMult(x[i], 1.0, A(i,i), s);
	}

	return true;
}

/////////////////////////////////////////////////////////////////////////////////////////////
//	gs_step_UR
/**
 * \brief Performs a backward gauss-seidel-step, that is, solve on the upper right of A.
 * \param A Matrix \f$A = D - L - U\f$
 * \param x will be \f$x = (D-U)^{-1}b \f$
 * \param b the vector b.
 * \sa gs_step_LL, sgs_step
 */
template<typename Matrix_type, typename Vector_type>
bool gs_step_UR(const Matrix_type &A, Vector_type &x, const Vector_type &b)
{
	// gs UR has preconditioning matrix
	// (D-U)^{-1}
	typename Vector_type::value_type s;

	if(x.size() == 0) return true;
	size_t i = x.size()-1;
	do
	{
		s = b[i];
		typename Matrix_type::const_row_iterator diag = A.get_connection(i, i);
		
		typename Matrix_type::const_row_iterator it = diag; ++it;
		for(; it != A.end_row(i); ++it)
			// s -= it.value() * x[it.index()];
			MatMultAdd(s, 1.0, s, -1.0, it.value(), x[it.index()]);

		// x[i] = s/A(i,i)
		InverseMatMult(x[i], 1.0, diag.value(), s);
	} while(i-- != 0);

	return true;
}

/////////////////////////////////////////////////////////////////////////////////////////////
//	sgs_step
/**
 * \brief Performs a symmetric gauss-seidel step
 * \param A Matrix \f$A = D - L - R\f$
 * \param x will be \f$x = (D-U)^{-1} D (D-L)^{-1} b \f$
 * \param b the vector b.
 * \sa gs_step_LL, gs_step_LL
 */
template<typename Matrix_type, typename Vector_type>
bool sgs_step(const Matrix_type &A, Vector_type &x, const Vector_type &b)
{
	// sgs has preconditioning matrix
	// (D-U)^{-1} D (D-L)^{-1}

	// x1 = (D-L)^{-1} b
	gs_step_LL(A, x, b);

	// x2 = D x1
	for(size_t i = 0; i<x.size(); i++)
		MatMult(x[i], 1.0, A(i, i), x[i]);

	// x3 = (D-U)^{-1} x2
	gs_step_UR(A, x, x);

	return true;
}

/////////////////////////////////////////////////////////////////////////////////////////////
//	diag_step
/**
 * \brief Performs a jacobi-step
 * \param A Matrix \f$A = D - L - R\f$
 * \param x will be \f$x = D^{-1} b \f$
 * \param b the vector b.
 * \param damp the damping factor
  */
template<typename Matrix_type, typename Vector_type>
bool diag_step(const Matrix_type& A, Vector_type& x, const Vector_type& b, number damp)
{
	//exit(3);
	UG_ASSERT(x.size() == b.size() && x.size() == A.num_rows(), x << ", " << b << " and " << A << " need to have same size.");

	for(size_t i=0; i < x.size(); i++)
		// x[i] = damp * b[i]/A(i,i)
		InverseMatMult(x[i], damp, A(i,i), b[i]);

	return true;
}

/// @}
}
#endif // __H__UG__CPU_ALGEBRA__CORE_SMOOTHERS__
