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

// predeclaration
template <int TDim>
class QuadratureRuleProvider;

// registering function
template <typename TRefElem>
bool RegisterQuadratureRule(QuadratureRuleProvider<TRefElem::dim>& factory);
template <> bool RegisterQuadratureRule<ReferenceVertex>(QuadratureRuleProvider<ReferenceVertex::dim>& factory);
template <> bool RegisterQuadratureRule<ReferenceEdge>(QuadratureRuleProvider<ReferenceEdge::dim>& factory);
template <> bool RegisterQuadratureRule<ReferenceTriangle>(QuadratureRuleProvider<ReferenceTriangle::dim>& factory);
template <> bool RegisterQuadratureRule<ReferenceQuadrilateral>(QuadratureRuleProvider<ReferenceQuadrilateral::dim>& factory);
template <> bool RegisterQuadratureRule<ReferenceTetrahedron>(QuadratureRuleProvider<ReferenceTetrahedron::dim>& factory);
template <> bool RegisterQuadratureRule<ReferencePrism>(QuadratureRuleProvider<ReferencePrism::dim>& factory);
template <> bool RegisterQuadratureRule<ReferencePyramid>(QuadratureRuleProvider<ReferencePyramid::dim>& factory);
template <> bool RegisterQuadratureRule<ReferenceHexahedron>(QuadratureRuleProvider<ReferenceHexahedron::dim>& factory);

// registering function
template <int dim>
bool RegisterQuadratureRuleDim(QuadratureRuleProvider<dim>& factory);

// implementation 0d
template <>
inline bool RegisterQuadratureRuleDim(QuadratureRuleProvider<0>& factory)
{
	bool bRet = true;
	bRet &= RegisterQuadratureRule<ReferenceVertex>(factory);
	return bRet;
}
// implementation 1d
template <>
inline bool RegisterQuadratureRuleDim(QuadratureRuleProvider<1>& factory)
{
	bool bRet = true;
	bRet &= RegisterQuadratureRule<ReferenceEdge>(factory);
	return bRet;
}
// implementation 2d
template <>
inline bool RegisterQuadratureRuleDim(QuadratureRuleProvider<2>& factory)
{
	bool bRet = true;
	bRet &= RegisterQuadratureRule<ReferenceTriangle>(factory);
	bRet &= RegisterQuadratureRule<ReferenceQuadrilateral>(factory);
	return bRet;
}
// implementation 3d
template <>
inline bool RegisterQuadratureRuleDim(QuadratureRuleProvider<3>& factory)
{
	bool bRet = true;
	bRet &= RegisterQuadratureRule<ReferenceTetrahedron>(factory);
	bRet &= RegisterQuadratureRule<ReferencePrism>(factory);
	bRet &= RegisterQuadratureRule<ReferencePyramid>(factory);
	bRet &= RegisterQuadratureRule<ReferenceHexahedron>(factory);
	return bRet;
}

/// provides quadrature rules for a reference dimension
/**
 * This class serves as a provider for quadrature rules. It is templated for a
 * reference element dimension.
 *
 * \tparam 	TDim	Reference Element Dimension
 */
template <int TDim>
class QuadratureRuleProvider
{
	public:
	///	dimension of reference element
		static const int dim = TDim;

	private:
	///	private constructor performing standard registering
		QuadratureRuleProvider()
		{
			m_vRule.resize(NUM_REFERENCE_OBJECTS);

		//	register standard rules
			RegisterQuadratureRuleDim<dim>(*this);
		}

	//	disallow copy
		QuadratureRuleProvider(const QuadratureRuleProvider&);
		QuadratureRuleProvider& operator=(const QuadratureRuleProvider&);

	//	provide rule
		static const QuadratureRule<dim>& get_quad_rule(ReferenceObjectID roid,
		                                                size_t order)
		{
		//	check if order or higerh order registered
			if(order >= m_vRule[roid].size())
				UG_THROW("QuadratureRuleProvider: Quadrature Rule not found for "
						<<roid<<" (dim="<<dim<<") and order "<<order);

		//	look for rule of order or next higher one
			if(m_vRule[roid][order] == NULL)
			{
				for(size_t i = order + 1; i < m_vRule[roid].size(); ++i)
				{
				//	return higher order than requested
					if(m_vRule[roid][i] != NULL) return *m_vRule[roid][i];
				}
				UG_THROW("QuadratureRuleProvider: Quadrature Rule not found for "
						<<roid<<" (dim="<<dim<<") and order "<<order);
			}

		//	return correct order
			return *m_vRule[roid][order];
		}

	///	singleton provider
		static QuadratureRuleProvider<dim>& instance()
		{
			static QuadratureRuleProvider<dim> inst;
			return inst;
		}

	private:
	///	Vector, holding all registered rules
		static std::vector<std::vector<const QuadratureRule<TDim>*> > m_vRule;

	public:
	///	register rule at this provider
	/**
	 * This function registers a quadrature rule at the Provider. If there is
	 * already a rule registered for the order, the rule is overwritten.
	 */
		static void register_rule(ReferenceObjectID roid,
		                          const QuadratureRule<dim>& rule)
		{
			m_vRule.resize(NUM_REFERENCE_OBJECTS);

		//	get order of rule to register
			size_t order = rule.order();

		//	resize vector if needed
			if(m_vRule[roid].size() <= order) m_vRule[roid].resize(order+1, NULL);

		//	set or override rule
			m_vRule[roid][order] = &rule;
		}

	///	register rule at this provider
	/**
	 * This function registers a quadrature rule at the Provider. If there is
	 * already a rule registered for the order, the rule is overwritten.
	 */
		template <typename TRefElem>
		static void register_rule(const QuadratureRule<dim>& rule)
		{
		//	check that dimension is correct
			if(TRefElem::dim != dim)
				UG_THROW("QuadratureRuleProvider: registering by reference"
						" element, but at provider of different dimension.");

		//	get reference object id
			ReferenceObjectID roid = TRefElem::REFERENCE_OBJECT_ID;

		//	forward request
			register_rule(roid, rule);
		}

	///	gets quadrature rule of requested order
	/**
	 * This function returns the next quadrature rule of order >= 'order'
	 * that is registered to this Provider. If no rule is found an
	 * Exception is thrown.
	 *
	 * \param[in]	order		Order of requested quadrature rule
	 * \tparam		TRefElem	Reference element type
	 */
		template <typename TRefElem>
		inline static const QuadratureRule<dim>& get_rule(size_t order)
		{
		//	check that dimension is correct
			if(TRefElem::dim != dim)
				UG_THROW("QuadratureRuleProvider: requesting by reference"
						" element, but at provider of different dimension.");

		//	get reference object id
			ReferenceObjectID roid = TRefElem::REFERENCE_OBJECT_ID;

		//	forward request
			return instance().get_quad_rule(roid, order);
		}

	///	gets quadrature rule of requested order
	/**
	 * This function returns the next quadrature rule of order >= 'order'
	 * that is registered to this Provider. If no rule is found an
	 * Exception is thrown.
	 *
	 * \param[in]	roid		Reference Object id
	 * \param[in]	order		Order of requested quadrature rule
	 */
		inline static const QuadratureRule<dim>& get_rule(ReferenceObjectID roid,
		                                                  size_t order)
		{
		//	forward request
			return instance().get_quad_rule(roid, order);
		}

};

// Init static member
template <int dim>
std::vector<std::vector<const QuadratureRule<dim>*> > QuadratureRuleProvider<dim>::m_vRule
	= std::vector<std::vector<const QuadratureRule<dim>*> >();

/// flexible order gauss quadrature
/**
 * Providing gauss quadrature for an reference element
 * \tparam 		TRefElem		Reference Element Type
 */
template <typename TRefElem>
class FlexGaussQuadrature
	: public QuadratureRule<TRefElem::dim>
{
	public:
	///	Constructor
		FlexGaussQuadrature(int order);

	///	Destructor
		~FlexGaussQuadrature() {}
};

/// @}

} // namespace ug

#endif /* __H__UG__LIB_DISC__QUADRATURE__ */
