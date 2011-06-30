/*
 * lagrangep1.h
 *
 *  Created on: 17.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISCRETIZATION__LOCAL_SHAPE_FUNCTION_SET__LAGRANGEP1__LAGRANGEP1__
#define __H__UG__LIB_DISCRETIZATION__LOCAL_SHAPE_FUNCTION_SET__LAGRANGEP1__LAGRANGEP1__

#include "../local_shape_function_set.h"
#include "../local_dof.h"

namespace ug{

/// \ingroup lib_discretization_local_shape_function_set
/// @{

/// Lagrange shape functions of first order
/**
 * This class provides Lagrangian Shape Functions of order 1 for a
 * Reference Element.
 * \tparam 	TRefElem		Reference Element Type
 */
template <typename TRefElem, int TOrder>
class LagrangeP1
{
	public:
	///	Reference Element type
		typedef TRefElem reference_element_type;

	///	Dimension, where shape functions are defined
		static const int dim = TRefElem::dim;

	///	Domain position type
		typedef MathVector<dim> position_type;

	///	Shape type
		typedef number shape_type;

	///	Gradient type
		typedef MathVector<dim> grad_type;

	///	Order of Shape functions
		static const size_t order = 1;

	/// Number of shape functions
		static const size_t nsh = TRefElem::num_corners;

	public:
	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		size_t num_sh() const { return nsh;}

	///	\copydoc ug::LocalShapeFunctionSet::position()
		bool position(size_t i, position_type& pos) const;

	///	\copydoc ug::LocalShapeFunctionSet::shape()
		shape_type shape(size_t i, const position_type& x) const;

	///	\copydoc ug::LocalShapeFunctionSet::shapes()
		void shapes(shape_type* sOut, const position_type& x) const;

	///	\copydoc ug::LocalShapeFunctionSet::grad()
		grad_type grad(size_t i, const position_type& x) const;

	///	\copydoc ug::LocalShapeFunctionSet::grads()
		void grads(grad_type* gOut, const position_type& x) const;

	///	\copydoc ug::LocalShapeFunctionSet::has_sh_on()
		bool has_sh_on(int d) const {return d == 0;}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh(int,size_t)
		size_t num_sh(int d, size_t id) const
			{if(d==0)return 1; else return 0;}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh(ReferenceObjectID)
		size_t num_sh(ReferenceObjectID type) const
			{if(type==ROID_VERTEX)return 1; else return 0;}

	protected:
	///	evaluate gradient of i'th shape function at position x
		void evaluate_grad(size_t i, const position_type& x,
		                   grad_type& value) const;
};

/// @}

} //namespace ug

// include implementation
#include "lagrangep1_impl.h"

#endif /* __H__UG__LIB_DISCRETIZATION__LOCAL_SHAPE_FUNCTION_SET__LAGRANGEP1__LAGRANGEP1__ */
