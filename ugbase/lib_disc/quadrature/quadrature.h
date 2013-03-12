/*
 * quadrature.h
 *
 *  Created on: 15.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__QUADRATURE__
#define __H__UG__LIB_DISC__QUADRATURE__

#include "../reference_element/reference_element.h"

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
		static const int dim = TDim;

	/// Position Type in Reference Element Space
		typedef MathVector<dim> position_type;

	///	Type of weights
		typedef number weight_type;

	public:
	///	number of integration points
		inline size_t size() const {return m_numPoints;}

	///	returns i'th integration point
		inline const position_type& point(size_t i) const
		{
			UG_ASSERT(i < size(), "Wrong index");
			return m_pvPoint[i];
		}

	///	returns all positions in an array of size()
		inline const position_type* points() const {return m_pvPoint;}

	///	return the i'th weight
		inline weight_type weight(size_t i) const
		{
			UG_ASSERT(i < size(), "Wrong index");
			return m_pvWeight[i];
		}

	/// returns all weights in an array of size()
		inline const weight_type* weights() const	{return m_pvWeight;}

	///	returns the order
		inline size_t order() const {return m_order;}

	protected:
		const position_type* m_pvPoint;	///< Integration points
		const weight_type* m_pvWeight; 	///< Weights
		size_t m_numPoints;				///< number of points
		int m_order;					///< Order of rule
};

/// @}

} // namespace ug

#endif /* __H__UG__LIB_DISC__QUADRATURE__ */
