/*
 * quadrature_provider.h
 *
 *  Created on: 15.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__QUADRATURE_PROVIDER__
#define __H__UG__LIB_DISC__QUADRATURE_PROVIDER__

#include "lib_grid/grid/grid_base_objects.h"
#include "quadrature.h"


///	types of quadratures
enum QuadType {
	BEST = 0,
	GAUSS,
	GAUSS_LEGENDRE,
	NEWTON_COTES,
	NUM_QUADRATURE_TYPES // always last
};

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

	private:
	///	private constructor performing standard registering
		QuadratureRuleProvider();

	//	disallow copy
		QuadratureRuleProvider(const QuadratureRuleProvider&);
		QuadratureRuleProvider& operator=(const QuadratureRuleProvider&);

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

	///	provide rule, try to create it if not already present
		static const QuadratureRule<TDim>&
		get_quad_rule(ReferenceObjectID roid, size_t order, QuadType type);

	///	creates rule at this provider
		static void create_rule(ReferenceObjectID roid, size_t order, QuadType type);

	///	rule creation, returns NULL if unavailable
	/// \{
		static const QuadratureRule<TDim>* create_gauss_rule(ReferenceObjectID roid, size_t order);
		static const QuadratureRule<TDim>* create_newton_cotes_rule(ReferenceObjectID roid, size_t order);
		static const QuadratureRule<TDim>* create_gauss_legendre_rule(ReferenceObjectID roid, size_t order);
	/// \}

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
		static const QuadratureRule<TDim>&
		get(size_t order, QuadType type = BEST);

	///	gets quadrature rule of requested order
	/**
	 * This function returns the next quadrature rule of order >= 'order'
	 * that is registered to this Provider. If no rule is found an
	 * Exception is thrown.
	 *
	 * \param[in]	roid		Reference Object id
	 * \param[in]	order		Order of requested quadrature rule
	 */
		static const QuadratureRule<TDim>&
		get(ReferenceObjectID roid, size_t order, QuadType type = BEST);
};

// Init static member
template <int dim>
std::vector<const QuadratureRule<dim>*> QuadratureRuleProvider<dim>::m_vRule[NUM_QUADRATURE_TYPES][NUM_REFERENCE_OBJECTS];

/// writes the Identifier to the output stream
std::ostream& operator<<(std::ostream& out,	const QuadType& v);

/// returns Identifier from string
QuadType GetQuadratureType(const std::string& name);

/// @}

} // namespace ug

#include "quadrature_provider_impl.h"

#endif /* __H__UG__LIB_DISC__QUADRATURE_PROVIDER__ */
