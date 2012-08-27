/*
 * finite_volume_bridge.cpp
 *
 *  Created on: 24.08.2012
 *      Author: andreasvogel
 */

// extern headers
#include <iostream>
#include <sstream>
#include <string>

// include bridge
#include "bridge/bridge.h"
#include "bridge/util.h"

// lib_disc includes
#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/function_spaces/approximation_space.h"

#include "lib_disc/spatial_disc/disc_util/finite_volume_output.h"

using namespace std;

namespace ug{
namespace bridge{
namespace FiniteVolume{

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

}

template <template <class, int> class TFVGeom, typename TDomain>
static void DomainFVGeom(Registry& reg, string grp, string append)
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
template <typename TDomain>
static void Domain(Registry& reg, string grp)
{
	string suffix = GetDomainSuffix<TDomain>();
	string tag = GetDomainTag<TDomain>();

//	CreateSubControlVolumeFaceDomain
	{
		DomainFVGeom<FV1Geometry, TDomain>(reg, grp, "_FV1");
		DomainFVGeom<HFV1Geometry, TDomain>(reg, grp, "_HFV1");
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
{
}

}; // end Functionality
}// end FiniteVolume

void RegisterBridge_FiniteVolume(Registry& reg, string grp)
{
	grp.append("/Discretization");
	typedef FiniteVolume::Functionality Functionality;

	try{
		RegisterCommon<Functionality>(reg,grp);
		RegisterDimensionDependent<Functionality>(reg,grp);
		RegisterDomainDependent<Functionality>(reg,grp);
		RegisterAlgebraDependent<Functionality>(reg,grp);
		RegisterDomainAlgebraDependent<Functionality>(reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

}// namespace bridge
}// namespace ug
