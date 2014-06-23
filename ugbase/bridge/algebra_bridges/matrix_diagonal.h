/*
 * matrix_diagonal.h
 *
 *  Created on: 23.06.2014
 *      Author: mrupp
 */

#ifndef MATRIX_DIAGONAL_H_
#define MATRIX_DIAGONAL_H_

#include "lib_algebra/operator/interface/linear_iterator.h"

/**
 * D = MatrixDiagonal(mat) creates a LinearOperator which acts like D = diag(mat)
 */
namespace ug{
template <typename M, typename X, typename Y = X>
class MatrixDiagonal :	public virtual ILinearOperator<X,Y>,
						public M
{
	public:
	// 	Domain space
		typedef X domain_function_type;

	// 	Range space
		typedef Y codomain_function_type;

	// 	Matrix type
		typedef M matrix_type;
		typedef MatrixOperator<M, X, Y> mo_type;

		SmartPtr<mo_type> m_mo;

	public:
		MatrixDiagonal(SmartPtr<mo_type> mo)  {
			m_mo = mo;
		}

	// 	Init Operator J(u)
		virtual void init(const X& u)
		{
			m_mo->init(u);
		}

	// 	Init Operator L
		virtual void init() { m_mo->init(); }

	// 	Apply Operator f = L*u (e.g. d = J(u)*c in iterative scheme)
		virtual void apply(Y& f, const X& u)
		{
			matrix_type &A = m_mo->get_matrix();
			for(size_t i=0; i<A.num_rows(); i++)
				f[i] = A(i,i)*u[i];
		}

	// 	Apply Operator, i.e. f = f - L*u;
		virtual void apply_sub(Y& f, const X& u)
		{
			matrix_type &A = m_mo->get_matrix();
			for(size_t i=0; i<A.num_rows(); i++)
				f[i] -= A(i,i)*u[i];
		}
};


/**
 * Dinv = MatrixDiagonalInverse(mat) creates a LinearOperator which acts like Dinv = diag(mat)^{-1}
 */
template <typename M, typename X, typename Y = X>
class MatrixDiagonalInverse :	public virtual ILinearOperator<X,Y>,
						public M
{
	public:
	// 	Domain space
		typedef X domain_function_type;

	// 	Range space
		typedef Y codomain_function_type;

	// 	Matrix type
		typedef M matrix_type;
		typedef MatrixOperator<M, X, Y> mo_type;

		SmartPtr<mo_type> m_mo;

	public:
		MatrixDiagonalInverse(SmartPtr<mo_type> mo)  {
			m_mo = mo;
		}

	// 	Init Operator J(u)
		virtual void init(const X& u)
		{
			m_mo->init(u);
		}

	// 	Init Operator L
		virtual void init() { m_mo->init(); }

	// 	Apply Operator f = L*u (e.g. d = J(u)*c in iterative scheme)
		virtual void apply(Y& f, const X& u)
		{
			matrix_type &A = m_mo->get_matrix();
			for(size_t i=0; i<A.num_rows(); i++)
				InverseMatMult(f[i], 1.0, A(i,i), u[i]);
		}

	// 	Apply Operator, i.e. f = f - L*u;
		virtual void apply_sub(Y& f, const X& u)
		{
			matrix_type &A = m_mo->get_matrix();
			typename X::value_type t;
			for(size_t i=0; i<A.num_rows(); i++)
			{
				InverseMatMult(t, 1.0, A(i,i), u[i]);
				f[i] -= t;
			}

		}
};

}
#endif /* MATRIX_DIAGONAL_H_ */
