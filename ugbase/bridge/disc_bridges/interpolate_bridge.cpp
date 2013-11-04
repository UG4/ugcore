/*
 * interpolate_bridge.cpp
 *
 *  Created on: 22.09.2010
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
#include "lib_disc/function_spaces/interpolate.h"

using namespace std;

namespace ug{
namespace bridge{
namespace Interpolate{

/**
 * \defgroup interpolate_bridge Interpolation Bridge
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

//	typedef
	static const int dim = TDomain::dim;
	typedef ug::GridFunction<TDomain, TAlgebra> TFct;

//	Interpolate
	{
		reg.add_function("Interpolate", static_cast<void (*)(SmartPtr<UserData<number, dim> >, SmartPtr<TFct>, const char*, const char*, number)>(&ug::Interpolate<TFct>),grp, "Integral", "Data#GridFunction#Component#Subsets#Time");
		reg.add_function("Interpolate", static_cast<void (*)(SmartPtr<UserData<number, dim> >, SmartPtr<TFct>, const char*, number)>(&ug::Interpolate<TFct>),grp, "Integral", "Data#GridFunction#Component#Time");
		reg.add_function("Interpolate", static_cast<void (*)(SmartPtr<UserData<number, dim> >, SmartPtr<TFct>, const char*, const char*)>(&ug::Interpolate<TFct>), grp, "Integral", "Data#GridFunction#Component#Subsets");
		reg.add_function("Interpolate", static_cast<void (*)(SmartPtr<UserData<number, dim> >, SmartPtr<TFct>, const char*)>(&ug::Interpolate<TFct>),grp, "Integral", "Data#GridFunction#Component");

		reg.add_function("Interpolate", static_cast<void (*)(number, SmartPtr<TFct>, const char*, const char*, number)>(&ug::Interpolate<TFct>),grp, "Integral", "ConstantValue#GridFunction#Component#Subsets#Time");
		reg.add_function("Interpolate", static_cast<void (*)(number, SmartPtr<TFct>, const char*, number)>(&ug::Interpolate<TFct>),grp, "Integral", "ConstantValue#GridFunction#Component#Time");
		reg.add_function("Interpolate", static_cast<void (*)(number, SmartPtr<TFct>, const char*, const char*)>(&ug::Interpolate<TFct>), grp, "Integral", "ConstantValue#GridFunction#Component#Subsets");
		reg.add_function("Interpolate", static_cast<void (*)(number, SmartPtr<TFct>, const char*)>(&ug::Interpolate<TFct>),grp, "Integral", "ConstantValue#GridFunction#Component");

		#ifdef UG_FOR_LUA
		reg.add_function("Interpolate", static_cast<void (*)(const char*, SmartPtr<TFct>, const char*, const char*, number)>(&ug::Interpolate<TFct>),grp, "Integral", "LuaFunction#GridFunction#Component#Subsets#Time");
		reg.add_function("Interpolate", static_cast<void (*)(const char*, SmartPtr<TFct>, const char*, number)>(&ug::Interpolate<TFct>),grp, "Integral", "LuaFunction#GridFunction#Component#Time");
		reg.add_function("Interpolate", static_cast<void (*)(const char*, SmartPtr<TFct>, const char*, const char*)>(&ug::Interpolate<TFct>), grp, "Integral", "LuaFunction#GridFunction#Component#Subsets");
		reg.add_function("Interpolate", static_cast<void (*)(const char*, SmartPtr<TFct>, const char*)>(&ug::Interpolate<TFct>),grp, "Integral", "LuaFunction#GridFunction#Component");
		#endif
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

// end group interpolate_bridge
/// \}

}// namespace Interpolate

/// \addtogroup interpolate_bridge
void RegisterBridge_Interpolate(Registry& reg, string grp)
{
	grp.append("/Discretization");
	typedef Interpolate::Functionality Functionality;

	try{
//		RegisterCommon<Functionality>(reg,grp);
//		RegisterDimensionDependent<Functionality>(reg,grp);
//		RegisterDomainDependent<Functionality>(reg,grp);
//		RegisterAlgebraDependent<Functionality>(reg,grp);
		RegisterDomainAlgebraDependent<Functionality>(reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

}//	end of namespace bridge
}//	end of namespace ug
