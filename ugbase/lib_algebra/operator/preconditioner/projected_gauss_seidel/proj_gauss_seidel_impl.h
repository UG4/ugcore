/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Raphael Prohl
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


template <typename TDomain, typename TAlgebra>
void
ProjGaussSeidel<TDomain,TAlgebra>::
step(const matrix_type& A, vector_type& c, const vector_type& d, const number relax)
{
	for(size_t i = 0; i < c.size(); i++)
	{
		forward_gs_step(c, A, d, i, relax);

		//	project correction on the subspace defined by the obstacle constraints
		this->project_correction(c[i], i);
	}
}

template <typename TDomain, typename TAlgebra>
void
ProjBackwardGaussSeidel<TDomain,TAlgebra>::
step(const matrix_type& A, vector_type& c, const vector_type& d, const number relax)
{
	if(c.size() == 0) return;
	size_t i = c.size()-1;
	do
	{
		backward_gs_step(c, A, d, i, relax);

		//	project correction on the subspace defined by the obstacle constraints
		this->project_correction(c[i], i);

	} while(i-- != 0);
}

template <typename TDomain, typename TAlgebra>
void
ProjSymmetricGaussSeidel<TDomain,TAlgebra>::
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

		//	project correction on the subspace defined by the obstacle constraints
		this->project_correction(c[i], i);
	}
}

} // end namespace ug

#endif /* __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_GAUSS_SEIDEL__PROJ_GAUSS_SEIDEL_IMPL__ */
