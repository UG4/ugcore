/*
 * p1conform.h
 *
 *  Created on: 17.02.2010
 *      Author: andreasvogel
 */

#ifndef P1CONFORM_H_
#define P1CONFORM_H_

#include "../local_shape_function_set.h"
#include "../local_dof_pattern.h"

namespace ug{

template <typename TRefElem>
class P1conform : public ug::LocalShapeFunctionSet<TRefElem>{
	protected:
		typedef TRefElem reference_element_type;

		// dimension, where shape functions are defined (i.e. reference element dimension)
		static const std::size_t dim = TRefElem::dim;

		// domain position type
		typedef MathVector<dim> position_type;

		// evaluation type
		typedef number shape_value_type;

		// gradient value type
		typedef MathVector<dim> grad_value_type;

		// number of shape functions
		static const std::size_t nsh = TRefElem::num_corners;

	public:
		uint num_shape_functions() const { return nsh;	}

		bool position_of_dof(int nrShapeFct, position_type& value) const;

		bool evaluate(int nrShapeFct, const position_type& locPos, shape_value_type& value) const;

		bool evaluate_grad(int nrShapeFct, const position_type& locPos, grad_value_type& value) const;

		const LocalDoFPattern<TRefElem>& local_dof_pattern() const;

		static P1conform<TRefElem>& inst();

	private:
		P1conform();

		static const uint _order = 1;
		LocalDoFPattern<TRefElem> m_ElementDoFPattern;

};

} //namespace ug

#include "p1conform_impl.h"

#endif /* P1CONFORM_H_ */
