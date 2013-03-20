/*
 * quadrature_provider.h
 *
 *  Created on: 15.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__QUADRATURE_PROVIDER__
#define __H__UG__LIB_DISC__QUADRATURE_PROVIDER__

#include "lib_grid/grid/geometric_base_objects.h"
#include "quadrature.h"

namespace ug{

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

	///	types of quadratures
		enum QuadratureType {
			BEST = 0,
			GAUSS,
			GAUSS_LEGENDRE,
			NEWTON_COTES,
			NUM_QUADRATURE_TYPES // always last
		};

	private:
	///	private constructor performing standard registering
		QuadratureRuleProvider();

	//	disallow copy
		QuadratureRuleProvider(const QuadratureRuleProvider&);
		QuadratureRuleProvider& operator=(const QuadratureRuleProvider&);

	//	provide rule
		static const QuadratureRule<TDim>& get_quad_rule(ReferenceObjectID roid,
		                                                size_t order,
		                                                QuadratureType type);

	///	singleton provider
		static QuadratureRuleProvider<dim>& instance()
		{
			static QuadratureRuleProvider<dim> inst;
			return inst;
		}

	public:
	///	destructor
		~QuadratureRuleProvider();

	protected:
	///	Vector, holding all registered rules
		static std::vector<const QuadratureRule<TDim>*> m_vRule[NUM_QUADRATURE_TYPES][NUM_REFERENCE_OBJECTS];

	///	creates rule at this provider
	/**
	 * This function registers a quadrature rule at the Provider.
	 */
		static void create_rule(ReferenceObjectID roid,
		                        size_t order,
                                QuadratureType type);

		static const QuadratureRule<TDim>* create_gauss_rule(ReferenceObjectID roid, size_t order);
		static const QuadratureRule<TDim>* create_newton_cotes_rule(ReferenceObjectID roid, size_t order);
		static const QuadratureRule<TDim>* create_gauss_legendre_rule(ReferenceObjectID roid, size_t order);
		static const QuadratureRule<TDim>* create_gauss_jacobi10_rule(size_t order);
		static const QuadratureRule<TDim>* create_gauss_jacobi20_rule(size_t order);

	public:
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
		static const QuadratureRule<TDim>& get_rule(size_t order,
		                                            QuadratureType type = BEST);

	///	gets quadrature rule of requested order
	/**
	 * This function returns the next quadrature rule of order >= 'order'
	 * that is registered to this Provider. If no rule is found an
	 * Exception is thrown.
	 *
	 * \param[in]	roid		Reference Object id
	 * \param[in]	order		Order of requested quadrature rule
	 */
		static const QuadratureRule<TDim>& get_rule(ReferenceObjectID roid,
		                                            size_t order,
		                                            QuadratureType type = BEST);
};

// Init static member
template <int dim>
std::vector<const QuadratureRule<dim>*> QuadratureRuleProvider<dim>::m_vRule[NUM_QUADRATURE_TYPES][NUM_REFERENCE_OBJECTS];

/// writes the Identifier to the output stream
template <int TDim>
std::ostream& operator<<(std::ostream& out,	const typename QuadratureRuleProvider<TDim>::QuadratureType& v);

/// returns Identifier from string
template <int TDim>
typename QuadratureRuleProvider<TDim>::QuadratureType GetQuadratureType(const std::string& name);

/// @}

} // namespace ug

#include "quadrature_provider_impl.h"

#endif /* __H__UG__LIB_DISC__QUADRATURE_PROVIDER__ */
