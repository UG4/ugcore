/*
 * lagrangep1.h
 *
 *  Created on: 17.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__LAGRANGEP1__LAGRANGEP1__
#define __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__LAGRANGEP1__LAGRANGEP1__

#include "../local_finite_element_provider.h"
#include "../local_dof_set.h"
#include "lagrange_local_dof.h"

namespace ug{

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
	  public BaseLocalShapeFunctionSet<LagrangeP1<TRefElem>, TRefElem::dim>
{
	///	base class
		typedef BaseLocalShapeFunctionSet<LagrangeP1<TRefElem>, TRefElem::dim> base_type;

	public:
	///	Shape type
		typedef typename base_type::shape_type shape_type;

	///	Gradient type
		typedef typename base_type::grad_type grad_type;

	public:
	///	Reference Element type
		typedef TRefElem reference_element_type;

	///	Dimension, where shape functions are defined
		static const int dim = TRefElem::dim;

	///	Order of Shape functions
		static const size_t order = 1;

	/// Number of shape functions
		static const size_t nsh = TRefElem::numCorners;

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

#endif /* __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__LAGRANGEP1__LAGRANGEP1__ */
