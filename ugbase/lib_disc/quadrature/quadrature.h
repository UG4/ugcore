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

#ifndef __H__UG__LIB_DISC__QUADRATURE__
#define __H__UG__LIB_DISC__QUADRATURE__

#include "common/common.h"
#include "common/math/ugmath.h"

namespace ug{

// Doxygen group
////////////////////////////////////////////////////////////////////////
/**
 * \brief supply of quadrature rules.
 *
 * The Quadrature Rule section provides the user with several quadrature
 * rules for all reference elements.
 *
 * \defgroup lib_discretization_quadrature_rules Quadrature Rules
 * \ingroup lib_discretization
 */

/// \addtogroup lib_discretization_quadrature_rules
/// @{

/// provides quadrature rule for a Reference Dimension
/**
 * A Quadrature Rule provides for a given Reference Element integration points
 * and weights. An Integral over the Reference Element T is approximated by
 * \f[
 * 		\int\limits_T f(\mathbf{x}) \; d\mathbf{x} \approx \sum_{i=0}^{n-1}
 * 			f(\mathbf{x}_{i}) \cdot w_i
 * \f]
 * with the \f$n\f$ integration points \f$\mathbf{x}_i\f$ and weights
 * \f$ w_i \f$.
 *
 * \tparam 		TDim 		Dimension of Reference Element
 */
template <int TDim>
class QuadratureRule{
	public:
	///	Dimension of Reference Element
		static constexpr int dim = TDim;

	/// Position Type in Reference Element Space
	using position_type = MathVector<dim>;

	///	Type of weights
	using weight_type = number;

	public:
	///	constructor
		QuadratureRule() :
			m_pvPoint(nullptr), m_pvWeight(nullptr),
			m_numPoints(0), m_order(0)
		{}

	///	destructor
		virtual ~QuadratureRule() {}

	///	number of integration points
		inline size_t size() const {return m_numPoints;}

	///	returns i'th integration point
		inline const MathVector<dim>& point(size_t i) const
		{
			UG_ASSERT(i < size(), "Wrong index");
			return m_pvPoint[i];
		}

	///	returns all positions in an array of size()
		inline const MathVector<dim>* points() const {return m_pvPoint;}

	///	return the i'th weight
		inline number weight(size_t i) const
		{
			UG_ASSERT(i < size(), "Wrong index");
			return m_pvWeight[i];
		}

	/// returns all weights in an array of size()
		inline const number* weights() const	{return m_pvWeight;}

	///	returns the order
		inline size_t order() const {return m_order;}

	protected:
		const MathVector<dim>* m_pvPoint;	///< Integration points
		const number* m_pvWeight; 			///< Weights
		size_t m_numPoints;					///< number of points
		int m_order;						///< Order of rule
};

/// @}

} // namespace ug

#endif /* __H__UG__LIB_DISC__QUADRATURE__ */
