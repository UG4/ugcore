/*
 * quadrature.h
 *
 *  Created on: 15.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__QUADRATURE__
#define __H__LIBDISCRETIZATION__QUADRATURE__

#include "../reference_element/reference_element.h"

namespace ug{

/// Exception thrown when quadrature rule not found
struct UG_ERROR_QuadratureRuleNotRegistered
{
		UG_ERROR_QuadratureRuleNotRegistered(size_t order_)
			: order(order_)
		{}

		size_t order;
};

// predeclaration
template <typename TRefElem>
class QuadratureRuleProvider;

// registering function
template <typename TRefElem>
bool RegisterQuadratureRule(QuadratureRuleProvider<TRefElem>& factory);

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

/// provides quadrature rule for a given Reference Element
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
 * \tparam 		TRefElem 		Reference Element Type
 */
template <typename TRefElem>
class QuadratureRule{
	public:
	///	Reference Element Type
		typedef TRefElem reference_element_type;

	///	Dimension of Reference Element
		static const int dim = TRefElem::dim;

	/// Position Type in Reference Element Space
		typedef MathVector<dim> position_type;

	///	Type of weights
		typedef number weight_type;

	public:
	///	number of integration points
		inline size_t size() const {return m_num_points;}

	///	returns i'th integration point
		inline const position_type& point(size_t i) const
		{
			UG_ASSERT(i < m_num_points, "Wrong index");
			return m_points[i];
		}

	///	returns all positions in an array of size()
		inline position_type* points() const {return m_points;}

	///	return the i'th weight
		inline weight_type weight(size_t i) const
		{
			UG_ASSERT(i < m_num_points, "Wrong index");
			return m_weights[i];
		}

	/// returns all weights in an array of size()
		inline weight_type* weights() const	{return m_weights;}

	///	returns the order
		inline size_t order() const {return m_order;}

	protected:
		position_type* m_points;	///< Integration points
		weight_type* m_weights;	 	///< Weights
		size_t m_num_points;		///< number of points
		int m_order;				///< Order of rule
};

/// provides quadrature rules for an element type
/**
 * This class serves as a provider for quadrature rules. It is templated for a
 * reference element type.
 * \tparam 	TRefElem	Reference Element Type
 */
template <typename TRefElem>
class QuadratureRuleProvider{
	private:
	///	private constructor performing standard registering
		QuadratureRuleProvider()
		{
		//	register standard rules
			RegisterQuadratureRule<TRefElem>(*this);
		}

	//	disallow copy
		QuadratureRuleProvider(const QuadratureRuleProvider&);
		QuadratureRuleProvider& operator=(const QuadratureRuleProvider&);

	//	provide rule
		static const QuadratureRule<TRefElem>& get_quad_rule(size_t order)
		{
		//	check if order or higerh order registered
			if(order >= m_vRule.size())
				throw(UG_ERROR_QuadratureRuleNotRegistered(order));

		//	look for rule of order or next higher one
			if(m_vRule[order] == 0)
			{
				for(size_t i = order + 1; i < m_vRule.size(); ++i)
				{
				//	return higher order than requested
					if(m_vRule[i] != 0)
						return *m_vRule[i];
				}
				throw(UG_ERROR_QuadratureRuleNotRegistered(order));
			}

		//	return correct order
			return *m_vRule[order];
		}

	///	singleton provider
		static QuadratureRuleProvider<TRefElem>& instance()
		{
			static QuadratureRuleProvider<TRefElem> inst;
			return inst;
		}

	private:
	///	Vector, holding all registered rules
		static std::vector<const QuadratureRule<TRefElem>*> m_vRule;
		static bool m_initialized;

	public:
	///	register rule at this provider
	/**
	 * This function registers a quadrature rule at the Provider. If there is
	 * already a rule registered for the order, the rule is overwritten.
	 */
		static bool register_rule(const QuadratureRule<TRefElem>& rule)
		{
		//	get order of rule to register
			size_t order = rule.order();

		//	resize vector if needed
			if(m_vRule.size() <= order) m_vRule.resize(order+1, 0);

		//	set or override rule
			m_vRule[order] = &rule;

		//	we're done
			return true;
		}

	///	gets quadrature rule of requested order
	/**
	 * This function returns the next quadrature rule of order >= 'order'
	 * that is registered to this Provider. If no rule is found an
	 * Exception is thrown.
	 * \param[in]	order		Order of requested quadrature rule
	 */
		inline static const QuadratureRule<TRefElem>& get_rule(size_t order)
		{
			return instance().get_quad_rule(order);
		}
};


/// gauss quadrature
/**
 * Providing gauss quadrature for an reference element
 * \tparam 		TRefElem		Reference Element Type
 */
template <typename TRefElem>
class GaussQuadrature
	: public QuadratureRule<TRefElem>
{
	public:
		GaussQuadrature(int order);
		~GaussQuadrature();
	protected:
		inline bool allocate_memory(std::size_t n);
};

/// @}

} // namespace ug

#include "quadrature_impl.h"

#endif /* __H__LIBDISCRETIZATION__QUADRATURE__ */
