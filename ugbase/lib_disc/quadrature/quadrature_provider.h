/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
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
		static constexpr int dim = TDim;

	private:
	///	private constructor performing standard registering
		QuadratureRuleProvider();

	//	disallow copy
		QuadratureRuleProvider(const QuadratureRuleProvider&) = delete;
		QuadratureRuleProvider& operator = (const QuadratureRuleProvider&) = delete;

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
		static std::vector<const QuadratureRule<TDim>*> m_vRule[QuadType::NUM_QUADRATURE_TYPES][NUM_REFERENCE_OBJECTS];

	///	provide rule, try to create it if not already present
		static const QuadratureRule<TDim>&
		get_quad_rule(ReferenceObjectID roid, size_t order, QuadType type);

	///	creates rule at this provider
		static void create_rule(ReferenceObjectID roid, size_t order, QuadType type);

	///	rule creation, returns nullptr if unavailable
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
		get(size_t order, QuadType type = QuadType::BEST);

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
		get(ReferenceObjectID roid, size_t order, QuadType type = QuadType::BEST);
};

// Init static member
template <int dim>
std::vector<const QuadratureRule<dim>*> QuadratureRuleProvider<dim>::m_vRule[QuadType::NUM_QUADRATURE_TYPES][NUM_REFERENCE_OBJECTS];

/// writes the Identifier to the output stream
std::ostream& operator << (std::ostream& out,	const QuadType& v);

/// returns Identifier from string
QuadType GetQuadratureType(const std::string& name);

/// @}

} // namespace ug

#include "quadrature_provider_impl.h"

#endif