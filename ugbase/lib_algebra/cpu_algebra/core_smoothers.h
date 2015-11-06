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
