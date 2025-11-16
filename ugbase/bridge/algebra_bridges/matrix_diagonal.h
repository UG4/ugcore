/*
 * Copyright (c) 2014-2015:  G-CSC, Goethe University Frankfurt
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
		using domain_function_type = X;

	// 	Range space
		using codomain_function_type = Y;

	// 	Matrix type
		using matrix_type = M;
		using mo_type = MatrixOperator<M, X, Y>;

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
		using domain_function_type = X;

	// 	Range space
		using codomain_function_type = Y;

	// 	Matrix type
		using matrix_type = M;
		using mo_type = MatrixOperator<M, X, Y>;

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
