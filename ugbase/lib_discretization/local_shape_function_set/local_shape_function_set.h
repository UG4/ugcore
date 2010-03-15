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
#include "lib_grid/lib_grid.h"

// library intern headers
#include "lib_discretization/reference_element/reference_elements.h"
#include "lib_discretization/local_shape_function_set/local_dof_pattern.h"

namespace ug {

/*
 * Trial Function Interface
 *
 * This class can evaluate all trial functions for a given Geometric Object.
 * The use of terms is as follows:
 * We call a 'finite element' a geometric entity, that contains dofs. Specially dofs
 * can be located at the boundary of the geometric object, e.g. in vertices, edges on a triangle.
 * In contrast a geometric object is a grid entity. It can contain dofs, induced by the finite element it
 * belongs to.
 */
template <typename TRefElem>
class LocalShapeFunctionSet
{
	public:
	typedef TRefElem reference_element_type;

	// dimension, where shape functions are defined (i.e. reference element dimension)
	static const std::size_t dim = TRefElem::dim;

	// domain position type
	typedef MathVector<dim> position_type;

	// evaluation type
	typedef number shape_value_type;

	// gradient value type
	typedef MathVector<dim> grad_value_type;

public:
	// number of dofs on finite element
	virtual uint num_shape_functions() const = 0;

	// local position of dof, returns true if exists, returns false if no meaningful position possible
	virtual bool position_of_dof(int nrShapeFct, position_type& value) const = 0;

	// evaluate shape function
	virtual bool evaluate(int nrShapeFct, const position_type& locPos, shape_value_type& value) const = 0;

	// evaluate gradient
	virtual bool evaluate_grad(int nrShapeFct, const position_type& locPos, grad_value_type& value) const = 0;

	virtual const LocalDoFPattern<TRefElem>& local_dof_pattern() const = 0;

	// virtual destructor
	virtual ~LocalShapeFunctionSet()
	{};
};

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
