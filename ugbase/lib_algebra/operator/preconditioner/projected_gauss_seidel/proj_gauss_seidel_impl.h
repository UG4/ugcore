/*
 * proj_gauss_seidel_impl.h
 *
 *  Created on: 10.10.2013
 *      Author: raphaelprohl
 */

#ifndef __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_GAUSS_SEIDEL__PROJ_GAUSS_SEIDEL_IMPL__
#define __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_GAUSS_SEIDEL__PROJ_GAUSS_SEIDEL_IMPL__

#include "proj_gauss_seidel.h"

namespace ug{

///	commmon GaussSeidel-step-calls for a single index 'i'
template<typename Matrix_type, typename Vector_type>
void forward_gs_step(Vector_type& c, const Matrix_type& A, const Vector_type& d,
		const size_t i, const number relaxFactor)
{
	typename Vector_type::value_type s = d[i];

	for(typename Matrix_type::const_row_iterator it = A.begin_row(i);
			it != A.end_row(i) && it.index() < i; ++it)
		// s -= it.value() * x[it.index()];
		MatMultAdd(s, 1.0, s, -1.0, it.value(), c[it.index()]);

	//	c[i] = relaxFactor * s / A(i,i)
	InverseMatMult(c[i], relaxFactor, A(i,i), s);
}

template<typename Matrix_type, typename Vector_type>
void backward_gs_step(Vector_type& c, const Matrix_type& A, const Vector_type& d,
		const size_t i, const number relaxFactor)
{
	typename Vector_type::value_type s = d[i];

	typename Matrix_type::const_row_iterator diag = A.get_connection(i, i);
	typename Matrix_type::const_row_iterator it = diag; ++it;

	for(; it != A.end_row(i); ++it)
		// s -= it.value() * x[it.index()];
		MatMultAdd(s, 1.0, s, -1.0, it.value(), c[it.index()]);

	// c[i] = relaxFactor * s/A(i,i)
	InverseMatMult(c[i], relaxFactor, diag.value(), s);
}


template <typename TAlgebra>
void
ProjGaussSeidel<TAlgebra>::
step(const matrix_type& A, vector_type& c, const vector_type& d, const number relax)
{

	for(size_t i = 0; i < c.size(); i++)
	{
		forward_gs_step(c, A, d, i, relax);

		this->project_correction(c, i);
	}

}

template <typename TAlgebra>
void
ProjBackwardGaussSeidel<TAlgebra>::
step(const matrix_type& A, vector_type& c, const vector_type& d, const number relax)
{

	if(c.size() == 0) return;
	size_t i = c.size()-1;
	do
	{
		backward_gs_step(c, A, d, i, relax);

		this->project_correction(c, i);

	} while(i-- != 0);

}

template <typename TAlgebra>
void
ProjSymmetricGaussSeidel<TAlgebra>::
step(const matrix_type& A, vector_type& c, const vector_type& d, const number relax)
{

	for(size_t i = 0; i < c.size(); i++)
	{
		//	1. perform a forward GaussSeidel step
		//	c1 = (D-L)^{-1} d
		forward_gs_step(c, A, d, i, relax);

		//	2. c2 = D c1
		MatMult(c[i], 1.0, A(i, i), c[i]);

		//	3. perform a backward GaussSeidel step
		//	c3 = (D-U)^{-1} c2
		backward_gs_step(c, A, c, i, relax);

		this->project_correction(c, i);
	}
}

} // end namespace ug

#endif /* __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_GAUSS_SEIDEL__PROJ_GAUSS_SEIDEL_IMPL__ */
