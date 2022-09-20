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
#include "bridge/util_domain_algebra_dependent.h"

// lib_disc includes
#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/function_spaces/approximation_space.h"

#include "lib_disc/spatial_disc/disc_util/fv_output.h"

using namespace std;

namespace ug{
namespace bridge{
namespace FiniteVolume{

/**
 * \defgroup finitvolume_bridge Finite Volume Bridge
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
template <typename TDomain, typename TAlgebra, typename TRegistry=Registry>
static void DomainAlgebra(Registry& reg, string grp)
{
	string suffix = GetDomainAlgebraSuffix<TDomain,TAlgebra>();
	string tag = GetDomainAlgebraTag<TDomain,TAlgebra>();

}

template <template <class, int> class TFVGeom, typename TDomain, typename TRegistry=Registry>
static void DomainFVGeom(TRegistry& reg, string grp, string append)
{
	string suffix = GetDomainSuffix<TDomain>();
	string tag = GetDomainTag<TDomain>();

	std::string name = "";
	{
		name = std::string("CreateSubControlVolumeFaceDomain"); name.append(append);
		reg.add_function(name, &CreateSubControlVolumeFaceDomain<TFVGeom, TDomain>, grp);
		name = std::string("CreateSubControlVolumeDomain"); name.append(append);
		reg.add_function(name, &CreateSubControlVolumeDomain<TFVGeom, TDomain>, grp);
		name = std::string("CreateControlVolumeDomain");name.append(append);
		reg.add_function(name, &CreateControlVolumeDomain<TFVGeom, TDomain>, grp);
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
template <typename TDomain, typename TRegistry=Registry>
static void Domain(TRegistry& reg, string grp)
{
	string suffix = GetDomainSuffix<TDomain>();
	string tag = GetDomainTag<TDomain>();

//	CreateSubControlVolumeFaceDomain
	{
		DomainFVGeom<FV1Geometry, TDomain, TRegistry>(reg, grp, "_FV1");
		DomainFVGeom<HFV1Geometry, TDomain, TRegistry>(reg, grp, "_HFV1");
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
template <int dim, typename TRegistry=Registry>
static void Dimension(TRegistry& reg, string grp)
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
template <typename TAlgebra, typename TRegistry>
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
{
}

}; // end Functionality

// end group finitvolume_bridge
/// \}

}// end FiniteVolume

/// \addtogroup finitvolume_bridge
template <typename TRegistry=Registry>
void RegisterBridge_FiniteVolume_(TRegistry& reg, string grp)
{
	grp.append("/Discretization");
	typedef FiniteVolume::Functionality Functionality;

	try{
		//RegisterCommon<Functionality,TRegistry>(reg,grp);
		//RegisterDimensionDependent<Functionality,TRegistry>(reg,grp);
		RegisterDomainDependent<Functionality,TRegistry>(reg,grp);
		//RegisterAlgebraDependent<Functionality,TRegistry>(reg,grp);
		//RegisterDomainAlgebraDependent<Functionality,TRegistry>(reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

}// namespace bridge

UG_REGISTRY_DEFINE(RegisterBridge_FiniteVolume);

}// namespace ug
