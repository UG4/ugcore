/*
 * lagrangep1.h
 *
 *  Created on: 17.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISCRETIZATION__LOCAL_SHAPE_FUNCTION_SET__LAGRANGEP1__LAGRANGEP1__
#define __H__UG__LIB_DISCRETIZATION__LOCAL_SHAPE_FUNCTION_SET__LAGRANGEP1__LAGRANGEP1__

#include "../local_shape_function_set.h"
#include "../local_dof_pattern.h"

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
	: public ug::LocalShapeFunctionSet<TRefElem>
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
	///	constructor
		LagrangeP1(){m_ElementDoFPattern.set_num_dofs(ROID_VERTEX, 1);}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		virtual size_t num_sh() const { return nsh;}

	///	\copydoc ug::LocalShapeFunctionSet::position()
		virtual bool position(size_t i, position_type& pos) const;

	///	\copydoc ug::LocalShapeFunctionSet::shape()
		virtual shape_type shape(size_t i, const position_type& x) const;

	///	\copydoc ug::LocalShapeFunctionSet::shapes()
		virtual void shapes(shape_type* sOut, const position_type& x) const;

	///	\copydoc ug::LocalShapeFunctionSet::grad()
		virtual grad_type grad(size_t i, const position_type& x) const;

	///	\copydoc ug::LocalShapeFunctionSet::grads()
		virtual void grads(grad_type* gOut, const position_type& x) const;

		const LocalDoFPattern<TRefElem>& local_dof_pattern() const
		{
			return m_ElementDoFPattern;
		}

	protected:
	///	evaluate gradient of i'th shape function at position x
		void evaluate_grad(size_t i, const position_type& x,
		                   grad_type& value) const;

	private:
		LocalDoFPattern<TRefElem> m_ElementDoFPattern;
};

/// @}

} //namespace ug

// include implementation
#include "lagrangep1_impl.h"

#endif /* __H__UG__LIB_DISCRETIZATION__LOCAL_SHAPE_FUNCTION_SET__LAGRANGEP1__LAGRANGEP1__ */
