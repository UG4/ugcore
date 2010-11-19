/*
 * trialspace.h
 *
 *  Created on: 12.05.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__LOCAL_SHAPE_FUCNTION_SET__
#define __H__LIBDISCRETIZATION__LOCAL_SHAPE_FUCNTION_SET__

// extern libraries
#include <cassert>
#include <map>

// other ug4 modules
#include "common/math/ugmath.h"
#include "lib_grid/lg_base.h"

// library intern headers
#include "../reference_element/reference_element.h"
#include "./local_dof_pattern.h"

namespace ug {

// Exception
struct UG_ERROR_InvalidShapeFunctionIndex
{
	UG_ERROR_InvalidShapeFunctionIndex(size_t i)
		: index(i)
	{}
	size_t index;
};

// Doxygen group
////////////////////////////////////////////////////////////////////////
/**
 * \brief provides local shape function sets.
 *
 * The Local Shape Function Set section contains the shape functions that
 * can be used in discretization.
 *
 * \defgroup lib_discretization_local_shape_function_set Local Shape Function Sets
 * \ingroup lib_discretization
 */

/// \ingroup lib_discretization_local_shape_function_set
/// @{

// LocalShapeFunctionSet
/** base class for local shape functions
 *
 * This class is a base class for the supply of shape functions
 * on finite elements. It can evaluate all trial functions for a
 * given Geometric Object. We call a 'finite element' a geometric
 * entity, that contains degrees of freedoms (DoF). DoFs can be
 * located in the inner or at the boundary of the geometric object.
 * (E.g in vertices, edges on a triangle) In contrast a geometric
 * object is a grid entity. It can contain DoFs, induced by the
 * finite element it belongs to. We distinguish between shape functions
 * on different element types by a Reference Element template
 * argument.
 * \tparam 	TRefElem	Reference Element Type
 */
template <typename TRefElem>
class LocalShapeFunctionSet
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

	public:
	///	Number of DoFs (shapes) on finite element
		virtual size_t num_sh() const = 0;

	///	local position of DoF i
	/**
	 * This function returns the local position of a DoF if possible.
	 * \param[in] 	i		number of DoF
	 * \param[out]	pos		Position of DoF
	 * \retval		true 	if position exists
	 * \retval		false 	if no meaningful position available
	 */
		virtual bool position(size_t i, position_type& pos) const = 0;

	/// evaluates the shape function
	/**
	 * This function returns the value of Shape Function i at
	 * an element-local evaluation point.
	 * \param[in]	i		number of DoF
	 * \param[in]	x		Position on reference element (evaluation point)
	 * \return shape function value at point
	 */
		virtual shape_type shape(size_t i, const position_type& x) const = 0;

	/// returns all shape functions evaluated at a point
	/**
	 * This function returns the values of all Shape Functions at
	 * an element-local evaluation point in an array.
	 * \param[out]	sOut	Vector of Shapes
	 * \param[in]	x		Position on reference element (evaluation point)
	 */
		virtual void shapes(shape_type* sOut, const position_type& x) const = 0;

	/// evaluates the gradient of the shape function
	/** This function returns the gradient of Shape Function i at
	 * an element-local evaluation point.
	 * \param[in]	i		number of DoF
	 * \param[in]	x		Position on reference element (evaluation point)
	 * \return gradient at point
	 */
		virtual grad_type grad(size_t i, const position_type& x) const = 0;

	/// returns all gradients evaluated at a point
	/**
	 * This function returns the gradients of all Shape Functions at
	 * an element-local evaluation point in an array.
	 * \param[out]	gOut	Vector of gradients
	 * \param[in]	x		Position on reference element (evaluation point)
	 */
		virtual void grads(grad_type* gOut, const position_type& x) const = 0;

	///	returns the local dof pattern for this element type
		virtual const LocalDoFPattern<TRefElem>& local_dof_pattern() const = 0;

	///	virtual destructor
		virtual ~LocalShapeFunctionSet()
		{};
};

/// @}

/** Traits for local shape function sets
 *
 * Every set has to supply the following info:
 *
 * - corresponding LocalDoFPattern as local_dof_pattern_type
 *
 */
class local_shape_function_set_traits {};

} // namespace ug

#endif /* __H__LIBDISCRETIZATION__LOCAL_SHAPE_FUCNTION_SET__ */
