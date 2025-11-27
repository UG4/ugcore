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

#ifndef __H__UG__CPU_ALGEBRA__CORE_SMOOTHERS__
#define __H__UG__CPU_ALGEBRA__CORE_SMOOTHERS__
////////////////////////////////////////////////////////////////////////////////////////////////

namespace ug
{

/// \addtogroup lib_algebra
///	@{

/////////////////////////////////////////////////////////////////////////////////////////////
///		Gauss-Seidel-Iterations
/**
 * Here, iteration schemes of gauss-seidel-type are described for solving the
 * linear equation
 *
 * 		\f$ A * x = b.			A \in R^{nxn}, x \in R^n, b \in R^n \f$.
 *
 * 	Most of the common linear iteration-methods base on the decomposition of A into
 * 	its diagonal (D) and strict-upper(-U) and strict-lower part (-L),
 *
 * 		\f$ A = D - L - U \f$.
 *
 * 	Among others, W. Hackbusch ('Iterative Loesung grosser Gleichungssysteme'),
 * 	distinguishes three different forms for describing a linear iteration scheme.
 * 	The general 'first normal-form' of a linear iteration scheme takes the form
 *
 * 		\f$ x^{m+1} = M * x^m + N * b \f$,
 *
 * 	with some Matrices \f$ M \f$ and \f$ N \in R^{nxn} \f$. m denotes the iteration index.
 * 	The general 'second normal-form' of a linear iteration scheme takes the form
 *
 * 		\f$ x^{m+1} = x^m - N * (A * x^m - b) \f$.
 *
 * 	Those linear iteration schemes, which can be represented by the second normal-form
 * 	are the linear, consistent iteration schemes.
 *
 * 	Introducing the correction \f$ c{m+1} := x^{m+1} - x^m \f$ and the defect
 * 	\f$ d^m := b - A * x^m \f$ the second normal-form can be rewritten as
 *
 * 		\f$ c = N * d \f$.
 *
 *	Below, methods for the (forward) Gauss-Seidel step, the backward Gauss-Seidel step and the symmetric
 *	Gauss-Seidel step are implemented ('gs_step_LL', 'gs_step_UR' resp. 'sgs_step').
 *	The matrices of the second normal-form are
 *
 *		\f$ N = (D-L)^{-1}\f$ 				for the (forward) Gauss-Seidel step,
 *		\f$ N = (D-U)^{-1}\f$ 				for the backward Gauss-Seidel step and
 *		\f$ N = (D-u)^{-1}D(D-U)^{-1} \f$	for the symmetric Gauss-Seidel step.
 *
 *	References:
 * <ul>
 * <li> W. Hackbusch. Iterative Loesung grosser Gleichungssysteme
 * </ul>
 *
 * \sa gs_step_LL, gs_step_UR, sgs_step
 */


/////////////////////////////////////////////////////////////////////////////////////////////
//	gs_step_LL
/** \brief Performs a forward gauss-seidel-step, that is, solve on the lower left of A.
 * Using gs_step_LL within a preconditioner-scheme leads to the fact that we get the correction
 * by successive inserting the already computed values of c in c = N * d, with c being the correction
 * and d being the defect. N denotes the matrix of the second normal-form of a linear iteration
 * scheme.
 *
* \param A Matrix \f$A = D - L - U\f$
* \param c Vector. \f$ c = N * d = (D-L)^{-1} * d \f$
* \param d Vector d.
* \sa gs_step_UR, sgs_step
*/
template<typename Matrix_type, typename Vector_type>
void gs_step_LL(const Matrix_type &A, Vector_type &c, const Vector_type &d, const number relaxFactor)
{
	// gs LL has preconditioning matrix N = (D-L)^{-1}

	using matrix_block = typename Matrix_type::value_type;
	using const_row_it = typename Matrix_type::const_row_iterator;
	typename Vector_type::value_type s;

	const size_t sz = c.size();
	for (size_t i = 0; i < sz; ++i)
	{
		s = d[i];

		//	loop over all lower left matrix entries.
		//	Note: Here the corrections c, which have already been computed in previous loops (wrt. i),
		//	are taken to compute the i-th correction. For example the correction of the second row
		//	is computed by s[2] = (d[2] - A[2][1] * c[1]); and c[2] = s[2]/A[2][2];
		const const_row_it rowEnd = A.end_row(i);
		const_row_it it = A.begin_row(i);
		for(; it != rowEnd && it.index() < i; ++it)
			// s -= it.value() * c[it.index()];
			MatMultAdd(s, 1.0, s, -1.0, it.value(), c[it.index()]);

		// c[i] = relaxFactor * s/A(i,i)
		const matrix_block& A_ii = it.index() == i ? it.value() : matrix_block(0);
		InverseMatMult(c[i], relaxFactor, A_ii, s);
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////
//	gs_step_UR
/**
 * \brief Performs a backward gauss-seidel-step, that is, solve on the upper right of A.
 * Using gs_step_UR within a preconditioner-scheme leads to the fact that we get the correction
 * by successive inserting the already computed values of c in c = N * d, with c being the correction
 * and d being the defect. N denotes the matrix of the second normal-form of a linear iteration
 * scheme.
 *
 * \param A Matrix \f$A = D - L - U\f$
 * \param c will be \f$c = N * d = (D-U)^{-1} * d \f$
 * \param d the vector d.
 * \sa gs_step_LL, sgs_step
 */
template<typename Matrix_type, typename Vector_type>
void gs_step_UR(const Matrix_type &A, Vector_type &c, const Vector_type &d, const number relaxFactor)
{
	// gs UR has preconditioning matrix N = (D-U)^{-1}

	typename Vector_type::value_type s;

	if(c.size() == 0) return;
	size_t i = c.size()-1;
	do
	{
		s = d[i];
		typename Matrix_type::const_row_iterator diag = A.get_connection(i, i);
		
		typename Matrix_type::const_row_iterator it = diag; ++it;
		for(; it != A.end_row(i); ++it)
			// s -= it.value() * x[it.index()];
			MatMultAdd(s, 1.0, s, -1.0, it.value(), c[it.index()]);

		// c[i] = relaxFactor * s/A(i,i)
		InverseMatMult(c[i], relaxFactor, diag.value(), s);
	} while(i-- != 0);

}

/////////////////////////////////////////////////////////////////////////////////////////////
//	sgs_step
/**
 * \brief Performs a symmetric gauss-seidel step.
 * Using sgs_step within a preconditioner-scheme leads to the fact that we get the correction
 * by successive inserting the already computed values of c in c = N * d, with c being the correction
 * and d being the defect. N denotes the matrix of the second normal-form of a linear iteration
 * scheme.
 *
 * \param A Matrix \f$A = D - L - R\f$
 * \param c will be \f$c = N * d = (D-U)^{-1} D (D-L)^{-1} d \f$
 * \param d the vector d.
 * \sa gs_step_LL, gs_step_LL
 */
template<typename Matrix_type, typename Vector_type>
void sgs_step(const Matrix_type &A, Vector_type &c, const Vector_type &d, const number relaxFactor)
{
	// sgs has preconditioning matrix N = (D-U)^{-1} D (D-L)^{-1}

	// c1 = (D-L)^{-1} d
	gs_step_LL(A, c, d, relaxFactor);

	// c2 = D c1
	typename Vector_type::value_type s;
	for(size_t i = 0; i<c.size(); i++)
	{
		s=c[i];
		MatMult(c[i], 1.0, A(i, i), s);
	}

	// c3 = (D-U)^{-1} c2
	gs_step_UR(A, c, c, relaxFactor);
}

/////////////////////////////////////////////////////////////////////////////////////////////
//	diag_step
/***
 * \brief Performs a jacobi-step
 * \param A Matrix \f$A = D - L - R\f$
 * \param c will be \f$c = N * d = D^{-1} d \f$
 * \param d the vector d.
 * \param damp the damping factor
  */
template<typename Matrix_type, typename Vector_type>
void diag_step(const Matrix_type& A, Vector_type& c, const Vector_type& d, number damp)
{
	//exit(3);
	UG_ASSERT(c.size() == d.size() && c.size() == A.num_rows(), c << ", " << d <<
			" and " << A << " need to have same size.");

	for(size_t i=0; i < c.size(); i++)
		// c[i] = damp * d[i]/A(i,i)
		InverseMatMult(c[i], damp, A(i,i), d[i]);
}


/// @}
}
#endif