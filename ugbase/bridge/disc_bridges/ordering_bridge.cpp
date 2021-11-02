/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
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

// extern headers
#include <iostream>
#include <sstream>
#include <string>

// include bridge
#include "bridge/bridge.h"
#include "bridge/util.h"
#include "bridge/util_domain_dependent.h"
#include "bridge/util_domain_algebra_dependent.h"

// lib_disc includes
#include "lib_disc/domain.h"

// ordering algorithms
#include "lib_disc/ordering_strategies/algorithms/ordering_algorithms.cpp"

using namespace std;

namespace ug{
namespace bridge{
namespace Ordering{

/**
 * \defgroup ordering_bridge Ordering Bridge
 * \ingroup disc_bridge
 * \{
 */

/**
 * Class exporting the functionality. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */
struct Functionality
{

/**
 * Function called for the registration of Domain and Algebra dependent parts.
 * All Functions and Classes depending on both Domain and Algebra
 * are to be placed here when registering. The method is called for all
 * available Domain and Algebra types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TDomain, typename TAlgebra>
static void DomainAlgebra(Registry& reg, string grp)
{
	string suffix = GetDomainAlgebraSuffix<TDomain,TAlgebra>();
	string tag = GetDomainAlgebraTag<TDomain,TAlgebra>();

//	typedefs for this algebra
	typedef typename TAlgebra::vector_type vector_type;
	typedef typename TAlgebra::matrix_type matrix_type;
	typedef std::vector<size_t> ordering_container_type;

//	Lexicographic ordering
	{
		typedef LexOrdering<TDomain, typename TAlgebra::matrix_type, ordering_container_type> T;
		typedef IOrderingAlgorithm<typename TAlgebra::matrix_type, ordering_container_type> TBase;
		string name = string("LexOrdering").append(suffix);
		reg.add_class_<T, TBase>(name, grp, "LexOrdering")
			.add_constructor()
			.add_method("set_approximation_space", &T::set_approximation_space)
			.add_method("set_direction", &T::set_direction)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "LexOrdering", tag);
	}
}

/**
 * Function called for the registration of Domain dependent parts.
 * All Functions and Classes depending on the Domain
 * are to be placed here when registering. The method is called for all
 * available Domain types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TDomain>
static void Domain(Registry& reg, string grp)
{
	string suffix = GetDomainSuffix<TDomain>();
	string tag = GetDomainTag<TDomain>();
	typedef ApproximationSpace<TDomain> approximation_space_type;

//	group string
	grp.append("/ApproximationSpace");

//	Order Cuthill-McKee
	{
		reg.add_function("OrderCuthillMcKee", static_cast<void (*)(approximation_space_type&, bool)>(&OrderCuthillMcKee), grp);
	}

//	Order lexicographically
	{
		reg.add_function("OrderLex", static_cast<void (*)(approximation_space_type&, const char*)>(&OrderLex<TDomain>), grp);
	}
//	Order in downwind direction
	{
		reg.add_function("OrderDownwind", static_cast<void (*)(approximation_space_type&, SmartPtr<UserData<MathVector<TDomain::dim>, TDomain::dim> >)> (&ug::OrderDownwind<TDomain>), grp);
		reg.add_function("OrderDownwind", static_cast<void (*)(approximation_space_type&, const std::vector<number>&)>(&ug::OrderDownwind<TDomain>), grp);

        reg.add_function("OrderDownwind", static_cast<void (*)(approximation_space_type&, SmartPtr<UserData<MathVector<TDomain::dim>, TDomain::dim> >, number)> (&ug::OrderDownwind<TDomain>), grp);
        reg.add_function("OrderDownwind", static_cast<void (*)(approximation_space_type&, const std::vector<number>&, number)>(&ug::OrderDownwind<TDomain>), grp);

#ifdef UG_FOR_LUA
		reg.add_function("OrderDownwind", static_cast<void (*)(approximation_space_type&, const char*)>(&ug::OrderDownwind<TDomain>), grp);
        reg.add_function("OrderDownwind", static_cast<void (*)(approximation_space_type&, const char*, number)>(&ug::OrderDownwind<TDomain>), grp);
#endif
	}

}

/**
 * Function called for the registration of Dimension dependent parts.
 * All Functions and Classes depending on the Dimension
 * are to be placed here when registering. The method is called for all
 * available Dimension types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <int dim>
static void Dimension(Registry& reg, string grp)
{
	string suffix = GetDimensionSuffix<dim>();
	string tag = GetDimensionTag<dim>();

}

/**
 * Function called for the registration of Algebra dependent parts.
 * All Functions and Classes depending on Algebra
 * are to be placed here when registering. The method is called for all
 * available Algebra types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TAlgebra>
static void Algebra(Registry& reg, string grp)
{
	string suffix = GetAlgebraSuffix<TAlgebra>();
	string tag = GetAlgebraTag<TAlgebra>();

}

/**
 * Function called for the registration of Domain and Algebra independent parts.
 * All Functions and Classes not depending on Domain and Algebra
 * are to be placed here when registering.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
static void Common(Registry& reg, string grp)
{}

}; // end Functionality

// end group ordering_bridge
/// \}

}// namespace Ordering

/// \addtogroup ordering_bridge
void RegisterBridge_Ordering(Registry& reg, string grp)
{
	grp.append("/Discretization");
	typedef Ordering::Functionality Functionality;

	try{
//		RegisterCommon<Functionality>(reg,grp);
//		RegisterDimensionDependent<Functionality>(reg,grp);
		RegisterDomainDependent<Functionality>(reg,grp);
//		RegisterAlgebraDependent<Functionality>(reg,grp);
		RegisterDomainAlgebraDependent<Functionality>(reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

}//	end of namespace bridge
}//	end of namespace ug
