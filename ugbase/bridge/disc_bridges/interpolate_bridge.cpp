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
#include "lib_disc/function_spaces/interpolate.h"
#include "lib_disc/function_spaces/interpolate_inner.h"

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
		// function-string
		reg.add_function("Interpolate", static_cast<void (*)(const char*, SmartPtr<TFct>, const char*, const char*, number)>(&ug::Interpolate<TFct>),grp, "Integral", "LuaFunction#GridFunction#Component#Subsets#Time");
		reg.add_function("Interpolate", static_cast<void (*)(const char*, SmartPtr<TFct>, const char*, number)>(&ug::Interpolate<TFct>),grp, "Integral", "LuaFunction#GridFunction#Component#Time");
		reg.add_function("Interpolate", static_cast<void (*)(const char*, SmartPtr<TFct>, const char*, const char*)>(&ug::Interpolate<TFct>), grp, "Integral", "LuaFunction#GridFunction#Component#Subsets");
		reg.add_function("Interpolate", static_cast<void (*)(const char*, SmartPtr<TFct>, const char*)>(&ug::Interpolate<TFct>),grp, "Integral", "LuaFunction#GridFunction#Component");

		// function-handle
		reg.add_function("Interpolate", static_cast<void (*)(LuaFunctionHandle, SmartPtr<TFct>, const char*, const char*, number)>(&ug::Interpolate<TFct>),grp, "Integral", "LuaFunction#GridFunction#Component#Subsets#Time");
		reg.add_function("Interpolate", static_cast<void (*)(LuaFunctionHandle, SmartPtr<TFct>, const char*, number)>(&ug::Interpolate<TFct>),grp, "Integral", "LuaFunction#GridFunction#Component#Time");
		reg.add_function("Interpolate", static_cast<void (*)(LuaFunctionHandle, SmartPtr<TFct>, const char*, const char*)>(&ug::Interpolate<TFct>), grp, "Integral", "LuaFunction#GridFunction#Component#Subsets");
		reg.add_function("Interpolate", static_cast<void (*)(LuaFunctionHandle, SmartPtr<TFct>, const char*)>(&ug::Interpolate<TFct>),grp, "Integral", "LuaFunction#GridFunction#Component");
		#endif
	}

	//IterpolateDiff
	{

		reg.add_function("Interpolate", static_cast<void (*)(SmartPtr<UserData<number, dim> >, SmartPtr<TFct>, const char*, const MathVector<dim>&)>(&ug::Interpolate<TFct>),grp, "Integral", "Data#GridFunction#Component#MoveVector");
	#ifdef UG_FOR_LUA
		//reg.add_function("Interpolate", static_cast<void (*)(LuaFunctionHandle, SmartPtr<TFct>, const char*,const SmartPtr<CplUserData<MathVector<dim>, dim> >)>(&ug::Interpolate<TFct>),grp, "Integral", "LuaFunction#GridFunction#Component#DiffVector");
	#endif
	}

	// InterpolateInner
	{
		reg.add_function("InterpolateInner", static_cast<void (*)(SmartPtr<UserData<number, dim> >, SmartPtr<TFct>, const char*, const char*, number)>(&ug::InterpolateInner<TFct>),grp, "Integral", "Data#GridFunction#Component#Subsets#Time");
		reg.add_function("InterpolateInner", static_cast<void (*)(SmartPtr<UserData<number, dim> >, SmartPtr<TFct>, const char*, number)>(&ug::InterpolateInner<TFct>),grp, "Integral", "Data#GridFunction#Component#Time");
		reg.add_function("InterpolateInner", static_cast<void (*)(SmartPtr<UserData<number, dim> >, SmartPtr<TFct>, const char*, const char*)>(&ug::InterpolateInner<TFct>), grp, "Integral", "Data#GridFunction#Component#Subsets");
		reg.add_function("InterpolateInner", static_cast<void (*)(SmartPtr<UserData<number, dim> >, SmartPtr<TFct>, const char*)>(&ug::InterpolateInner<TFct>),grp, "Integral", "Data#GridFunction#Component");

		reg.add_function("InterpolateInner", static_cast<void (*)(number, SmartPtr<TFct>, const char*, const char*, number)>(&ug::InterpolateInner<TFct>),grp, "Integral", "ConstantValue#GridFunction#Component#Subsets#Time");
		reg.add_function("InterpolateInner", static_cast<void (*)(number, SmartPtr<TFct>, const char*, number)>(&ug::InterpolateInner<TFct>),grp, "Integral", "ConstantValue#GridFunction#Component#Time");
		reg.add_function("InterpolateInner", static_cast<void (*)(number, SmartPtr<TFct>, const char*, const char*)>(&ug::InterpolateInner<TFct>), grp, "Integral", "ConstantValue#GridFunction#Component#Subsets");
		reg.add_function("InterpolateInner", static_cast<void (*)(number, SmartPtr<TFct>, const char*)>(&ug::InterpolateInner<TFct>),grp, "Integral", "ConstantValue#GridFunction#Component");

		#ifdef UG_FOR_LUA
		reg.add_function("InterpolateInner", static_cast<void (*)(const char*, SmartPtr<TFct>, const char*, const char*, number)>(&ug::InterpolateInner<TFct>),grp, "Integral", "LuaFunction#GridFunction#Component#Subsets#Time");
		reg.add_function("InterpolateInner", static_cast<void (*)(const char*, SmartPtr<TFct>, const char*, number)>(&ug::InterpolateInner<TFct>),grp, "Integral", "LuaFunction#GridFunction#Component#Time");
		reg.add_function("InterpolateInner", static_cast<void (*)(const char*, SmartPtr<TFct>, const char*, const char*)>(&ug::InterpolateInner<TFct>), grp, "Integral", "LuaFunction#GridFunction#Component#Subsets");
		reg.add_function("InterpolateInner", static_cast<void (*)(const char*, SmartPtr<TFct>, const char*)>(&ug::InterpolateInner<TFct>),grp, "Integral", "LuaFunction#GridFunction#Component");
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
