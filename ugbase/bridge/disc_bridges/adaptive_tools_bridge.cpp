/*
 * adaptive_tools_bridge.cpp
 *
 *  Created on: 06.03.2012
 *      Author: andreasvogel
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
#include "lib_disc/function_spaces/error_indicator.h"
#include "lib_disc/function_spaces/level_transfer.h"
#include "lib_disc/function_spaces/local_transfer.h"

using namespace std;

namespace ug{
namespace bridge{
namespace AdaptiveTools{

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

//	MarkForAdaption_GradientIndicator
	{
		string grp("ug4/Refinement/");
		reg.add_function("MarkForAdaption_GradientIndicator",
						 &MarkForAdaption_GradientIndicator<TDomain, SurfaceDoFDistribution, TAlgebra>, grp);
	}

//	Prolongate
	{
		reg.add_function("Prolongate",
						 &Prolongate<TDomain, SurfaceDoFDistribution, TAlgebra>, grp);
	}

//	Restrict
	{
		reg.add_function("Restrict",
						 &Restrict<TDomain, SurfaceDoFDistribution, TAlgebra>, grp);
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

//	ILocalTransferAlgebra
	{
		typedef ILocalTransferAlgebra<TAlgebra> T;
		typedef ILocalTransfer TBase;
		string name = string("ILocalTransferAlgebra").append(suffix);
		reg.add_class_<T, TBase>(name)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ILocalTransferAlgebra", tag);
	}

//	P1LocalTransfer
	{
		typedef P1LocalTransfer<TAlgebra> T;
		typedef ILocalTransferAlgebra<TAlgebra> TBase;
		string name = string("P1LocalTransfer").append(suffix);
		reg.add_class_<T, TBase>(name)
			.template add_constructor<void (*)(size_t)>("fct")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "P1LocalTransfer", tag);
	}

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
//	ILocalTransfer
	{
		reg.add_class_<ILocalTransfer>("ILocalTransfer");
	}
}

}; // end Functionality
}// end AdaptiveTools

void RegisterBridge_AdaptiveTools(Registry& reg, string grp)
{
	grp.append("/Discretization");
	typedef AdaptiveTools::Functionality Functionality;

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
