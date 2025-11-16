/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
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

#ifndef __H__LIB_ALGEBRA__OPERATOR__INTERFACE__MATRIX_OPERATOR__
#define __H__LIB_ALGEBRA__OPERATOR__INTERFACE__MATRIX_OPERATOR__

#include "linear_operator.h"
#include "lib_algebra/common/operations_mat/matrix_algebra_types.h"

namespace ug{


///////////////////////////////////////////////////////////////////////////////
// Matrix based linear operator
///////////////////////////////////////////////////////////////////////////////

template <typename M, typename X, typename Y = X>
class MatrixOperator :	public virtual ILinearOperator<X,Y>,
						public M
{
	public:
	// 	Domain space
		using domain_function_type = X;

	// 	Range space
		using codomain_function_type = Y;

	// 	Matrix type
		using matrix_type = M;

	public:
	// 	Init Operator J(u)
		virtual void init(const X& u) {}

	// 	Init Operator L
		virtual void init() {}

	// 	Apply Operator f = L*u (e.g. d = J(u)*c in iterative scheme)
		virtual void apply(Y& f, const X& u) {matrix_type::apply(f,u);}

	// 	Apply Operator, i.e. f = f - L*u;
		virtual void apply_sub(Y& f, const X& u) {matrix_type::matmul_minus(f,u);}

	// 	Access to matrix
		virtual M& get_matrix() {return *this;};
};

template<typename M, typename X, typename Y>
struct matrix_algebra_type_traits<MatrixOperator<M, X, Y> >
{
	enum
	{
		type = matrix_algebra_type_traits<M>::type
	};
};

} // end namespace ug
#endif /* __H__LIB_ALGEBRA__OPERATOR__INTERFACE__MATRIX_OPERATOR__ */
