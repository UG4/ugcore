/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__LAGRANGEP1__LAGRANGEP1__
#define __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__LAGRANGEP1__LAGRANGEP1__

#include "../local_finite_element_provider.h"
// #include "../local_dof_set.h"
#include "lagrange_local_dof.h"

namespace ug {

/// \ingroup lib_discretization_local_shape_function_set
/// @{

/// Lagrange shape functions of first order
/**
 * This class provides Lagrangian Shape Functions of order 1 for a
 * Reference Element.
 * \tparam 	TRefElem		Reference Element Type
 */
template <typename TRefElem>
class LagrangeP1
	: public LagrangeLDS<TRefElem>,
	  public BaseLSFS<LagrangeP1<TRefElem>, TRefElem::dim>
{
	///	base class
		using base_type = BaseLSFS<LagrangeP1, TRefElem::dim>;

	public:
	///	Shape type
		using shape_type = typename base_type::shape_type;

	///	Gradient type
		using grad_type = typename base_type::grad_type;

	public:
	///	Reference Element type
		using reference_element_type = TRefElem;

	///	Dimension, where shape functions are defined
		static constexpr int dim = TRefElem::dim;

	///	Order of Shape functions
		static constexpr size_t order = 1;

	/// Number of shape functions
		static constexpr size_t nsh = TRefElem::numCorners;

	public:
	///	constructor
		LagrangeP1() : LagrangeLDS<TRefElem>(1) {}

	///	\copydoc ug::LocalShapeFunctionSet::continuous()
		bool continuous() const {return true;}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		size_t num_sh() const { return nsh;}

	///	\copydoc ug::LocalShapeFunctionSet::position()
		bool position(size_t i, MathVector<dim>& pos) const;

	///	\copydoc ug::LocalShapeFunctionSet::shape()
		shape_type shape(size_t i, const MathVector<dim>& x) const;

	///	\copydoc ug::LocalShapeFunctionSet::grad()
		void grad(grad_type& value, size_t i, const MathVector<dim>& x) const;
};

/// @}

} //namespace ug

#endif