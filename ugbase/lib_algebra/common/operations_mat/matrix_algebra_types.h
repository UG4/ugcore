/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef MATRIX_ALGEBRA_TYPES_H_
#define MATRIX_ALGEBRA_TYPES_H_

namespace ug{

/*
suppose you have your algebra with your vector class yourvector and your matrix class yourmatrix.
when your matrix class is rowwise multiplicable, and you want to do so, set
template<>
class matrix_algebra_type_traits<yourclass>
{
	static constexpr matrix_algebra_type type = MATRIX_USE_ROW_FUNCTIONS;
};
and add a function
inline void yourmatrix::mat_mult_add_row(size_t row, yourvector::value_type &dest, double beta1, const yourvector &v) const;
to your class.

If you cannot or don't want to use rowwise multiplication, set
template<>
class matrix_algebra_type_traits<yourclass>
{
	static constexpr matrix_algebra_type type = MATRIX_USE_GLOBAL_FUNCTIONS;
};

and implement (at least)
inline void MatMultAddDirect(yourvector &dest, const number &beta1, const yourmatrix &A1, const yourvector &w1);
inline void MatMultDirect(yourvector &dest, const number &beta1, const yourmatrix &A1, const yourvector &w1,
							number &alpha1, const yourvector &v1);
Other functions like MatMultAdd(dest, beta1, A1, w1, beta2, A2, w2, alpha1, v1) are then implemented through your functions.
(see class mat_operations_class<vector_t, matrix_t, MATRIX_USE_GLOBAL_FUNCTIONS> for all functions).

*/

enum matrix_algebra_type
{
	MATRIX_USE_ROW_FUNCTIONS,
	MATRIX_USE_GLOBAL_FUNCTIONS,
	MATRIX_USE_OPERATORS,
	MATRIX_USE_MEMBER_FUNCTIONS
};

template<typename vector_t, typename matrix_t, int type>
struct mat_operations_class;

template<typename T>
struct matrix_algebra_type_traits
{
	static constexpr int type = MATRIX_USE_OPERATORS;
};


}

#endif