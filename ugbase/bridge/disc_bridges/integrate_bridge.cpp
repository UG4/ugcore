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
#include "lib_disc/function_spaces/integrate.h"
#include "lib_disc/function_spaces/integrate_flux.h"

#include "lib_disc/quadrature/quad_test.h"

using namespace std;

namespace ug{
namespace bridge{
namespace Integrate{

/**
 * \defgroup integrate_bridge Integration Bridge
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
 * @param grp				group for sorting of functionality
 */
template <typename TDomain, typename TAlgebra>
static void DomainAlgebra(Registry& reg, string grp)
{
	static constexpr int dim = TDomain::dim;
	using TFct = GridFunction<TDomain, TAlgebra>;

//	Integral
	{
		reg.add_function("Integral", static_cast<number (*)(SmartPtr<UserData<number,dim> >, SmartPtr<TFct>, const char*, number, int)>(&Integral<TFct>), grp, "Integral", "Data#GridFunction#Subsets#Time#QuadOrder");
		reg.add_function("Integral", static_cast<number (*)(SmartPtr<UserData<number,dim> >, SmartPtr<TFct>, const char*, number)>(&Integral<TFct>), grp, "Integral", "Data#GridFunction#Subsets#Time");
		reg.add_function("Integral", static_cast<number (*)(SmartPtr<UserData<number,dim> >, SmartPtr<TFct>, number)>(&Integral<TFct>), grp, "Integral", "Data#GridFunction#Time");
		reg.add_function("Integral", static_cast<number (*)(SmartPtr<UserData<number,dim> >, SmartPtr<TFct>, const char*)>(&Integral<TFct>), grp, "Integral", "Data#GridFunction#Subsets");
		reg.add_function("Integral", static_cast<number (*)(SmartPtr<UserData<number,dim> >, SmartPtr<TFct>)>(&Integral<TFct>), grp, "Integral", "Data#GridFunction");

		reg.add_function("Integral", static_cast<number (*)(number, SmartPtr<TFct>, const char*, number, int)>(&Integral<TFct>), grp, "Integral", "ConstantValue#GridFunction#Subsets#Time#QuadOrder");
		reg.add_function("Integral", static_cast<number (*)(number, SmartPtr<TFct>, const char*, number)>(&Integral<TFct>), grp, "Integral", "ConstantValue#GridFunction#Subsets#Time");
		reg.add_function("Integral", static_cast<number (*)(number, SmartPtr<TFct>, number)>(&Integral<TFct>), grp, "Integral", "ConstantValue#GridFunction#Time");
		reg.add_function("Integral", static_cast<number (*)(number, SmartPtr<TFct>, const char*)>(&Integral<TFct>), grp, "Integral", "ConstantValue#GridFunction#Subsets");
		reg.add_function("Integral", static_cast<number (*)(number, SmartPtr<TFct>)>(&Integral<TFct>), grp, "Integral", "ConstantValue#GridFunction");

#ifdef UG_FOR_LUA
		reg.add_function("Integral", static_cast<number (*)(const char*, SmartPtr<TFct>, const char*, number, int, std::string)>(&Integral<TFct>), grp, "Integral", "LuaFunction#GridFunction#Subsets#Time#QuadOrder#QuadType");
		reg.add_function("Integral", static_cast<number (*)(const char*, SmartPtr<TFct>, const char*, number, int)>(&Integral<TFct>), grp, "Integral", "LuaFunction#GridFunction#Subsets#Time#QuadOrder");
		reg.add_function("Integral", static_cast<number (*)(const char*, SmartPtr<TFct>, const char*, number)>(&Integral<TFct>), grp, "Integral", "LuaFunction#GridFunction#Subsets#Time");
		reg.add_function("Integral", static_cast<number (*)(const char*, SmartPtr<TFct>, number)>(&Integral<TFct>), grp, "Integral", "LuaFunction#GridFunction#Time");
		reg.add_function("Integral", static_cast<number (*)(const char*, SmartPtr<TFct>, const char*)>(&Integral<TFct>), grp, "Integral", "LuaFunction#GridFunction#Subsets");
		reg.add_function("Integral", static_cast<number (*)(const char*, SmartPtr<TFct>)>(&Integral<TFct>), grp, "Integral", "LuaFunction#GridFunction");
#endif

		reg.add_function("Integral",static_cast<number (*)(SmartPtr<TFct>, const char*, const char*, int)>(Integral<TFct>),grp, "Integral", "GridFunction#Component#Subsets#QuadOrder");
		reg.add_function("Integral",static_cast<number (*)(SmartPtr<TFct>, const char*, const char*)>(Integral<TFct>),grp, "Integral", "GridFunction#Component#Subsets");
		reg.add_function("Integral",static_cast<number (*)(SmartPtr<TFct>, const char*)>(Integral<TFct>),grp, "Integral", "GridFunction#Component");
	}

	{
		reg.add_function("Minimum",static_cast<number (*)(SmartPtr<TFct>, const char*, const char*)>(&Minimum<TFct>), grp);
	}
//	L2Error
	{
		reg.add_function("L2Error",static_cast<number (*)(SmartPtr<TFct>, const char*, SmartPtr<TFct>, const char*, int, const char*)>(&L2Error<TFct>), grp);
		reg.add_function("L2Error",static_cast<number (*)(SmartPtr<TFct>, const char*, SmartPtr<TFct>, const char*, int)>(&L2Error<TFct>), grp);
		reg.add_function("L2Error",static_cast<number (*)(SmartPtr<UserData<number, dim> >, SmartPtr<TFct>, const char*, number, int, const char*)>(&L2Error<TFct>), grp);
		reg.add_function("L2Error",static_cast<number (*)(SmartPtr<UserData<number, dim> >, SmartPtr<TFct>, const char*, number, int)>(&L2Error<TFct>), grp);
#ifdef UG_FOR_LUA
		reg.add_function("L2Error",static_cast<number (*)(const char*, SmartPtr<TFct>, const char*, number, int, const char*)>(&L2Error<TFct>), grp);
		reg.add_function("L2Error",static_cast<number (*)(const char*, SmartPtr<TFct>, const char*, number, int)>(&L2Error<TFct>), grp);
#endif
	}

//	H1Error
	{
		reg.add_function("H1Error",static_cast<number (*)(SmartPtr<TFct>, const char*, SmartPtr<TFct>, const char*, int, const char*)>(&H1Error<TFct>), grp);
		reg.add_function("H1Error",static_cast<number (*)(SmartPtr<TFct>, const char*, SmartPtr<TFct>, const char*, int)>(&H1Error<TFct>), grp);
		reg.add_function("H1Error",static_cast<number (*)(SmartPtr<UserData<number, dim> >, SmartPtr<UserData<MathVector<dim>, dim> >, SmartPtr<TFct>, const char*, number, int, const char*)>(&H1Error<TFct>), grp);
		reg.add_function("H1Error",static_cast<number (*)(SmartPtr<UserData<number, dim> >, SmartPtr<UserData<MathVector<dim>, dim> >, SmartPtr<TFct>, const char*, number, int)>(&H1Error<TFct>), grp);
#ifdef UG_FOR_LUA
		reg.add_function("H1Error",static_cast<number (*)(const char*, const char*, SmartPtr<TFct>, const char*, number, int, const char*)>(&H1Error<TFct>), grp);
		reg.add_function("H1Error",static_cast<number (*)(const char*, const char*, SmartPtr<TFct>, const char*, number, int)>(&H1Error<TFct>), grp);
#endif
	}

//	L2Norm
	{
		reg.add_function("L2Norm",static_cast<number (*)(SmartPtr<TFct>, const char*, int, const char*)>(&L2Norm<TFct>),grp);
		reg.add_function("L2Norm",static_cast<number (*)(SmartPtr<TFct>, const char*, int)>(&L2Norm<TFct>),grp);

		reg.add_function("H1SemiNorm",static_cast<number (*)(SmartPtr<TFct>, const char*, int)>(&H1SemiNorm<TFct>),grp);
		reg.add_function("H1SemiNorm",static_cast<number (*)(SmartPtr<TFct>, const char*, int, const char*)>(&H1SemiNorm<TFct>),grp);

		reg.add_function("H1Norm",static_cast<number (*)(SmartPtr<TFct>, const char*, int)>(&H1Norm<TFct>),grp);
		reg.add_function("H1Norm",static_cast<number (*)(SmartPtr<TFct>, const char*, int, const char*)>(&H1Norm<TFct>),grp);
	}

//	IntegrateNormalGradientOnManifold
	{
		reg.add_function("IntegrateNormalGradientOnManifold",static_cast<number (*)(TFct&, const char*, const char*, const char*)>(&IntegrateNormalGradientOnManifold<TFct>),grp, "Integral", "GridFunction#Component#BoundarySubset#InnerSubset");
	}

//	IntegralNormalComponentOnManifold
	{
		reg.add_function("IntegralNormalComponentOnManifold", static_cast<number (*)(SmartPtr<UserData<MathVector<dim>,dim> >, SmartPtr<TFct>, const char*, const char*, number, int)>(&IntegralNormalComponentOnManifold<TFct>), grp, "Integral", "Data#GridFunction#BoundarySubsets#InnerSubsets#Time#QuadOrder");
		reg.add_function("IntegralNormalComponentOnManifold", static_cast<number (*)(SmartPtr<UserData<MathVector<dim>,dim> >, SmartPtr<TFct>, const char*, const char*, number)>(&IntegralNormalComponentOnManifold<TFct>), grp, "Integral", "Data#GridFunction#BoundarySubsets#InnerSubsets#Time");
		reg.add_function("IntegralNormalComponentOnManifold", static_cast<number (*)(SmartPtr<UserData<MathVector<dim>,dim> >, SmartPtr<TFct>, const char*, number)>(&IntegralNormalComponentOnManifold<TFct>), grp, "Integral", "Data#GridFunction#BoundarySubsets#Time");
		reg.add_function("IntegralNormalComponentOnManifold", static_cast<number (*)(SmartPtr<UserData<MathVector<dim>,dim> >, SmartPtr<TFct>, const char*, const char*)>(&IntegralNormalComponentOnManifold<TFct>), grp, "Integral", "Data#GridFunction#BoundarySubsets#InnerSubsets");
		reg.add_function("IntegralNormalComponentOnManifold", static_cast<number (*)(SmartPtr<UserData<MathVector<dim>,dim> >, SmartPtr<TFct>, const char*)>(&IntegralNormalComponentOnManifold<TFct>), grp, "Integral", "Data#GridFunction#BoundarySubsets");
#ifdef UG_FOR_LUA
		reg.add_function("IntegralNormalComponentOnManifold", static_cast<number (*)(const char*, SmartPtr<TFct>, const char*, const char*, number, int)>(&IntegralNormalComponentOnManifold<TFct>), grp, "Integral", "LuaData#GridFunction#BoundarySubsets#InnerSubsets#Time#QuadOrder");
		reg.add_function("IntegralNormalComponentOnManifold", static_cast<number (*)(const char*, SmartPtr<TFct>, const char*, const char*, number)>(&IntegralNormalComponentOnManifold<TFct>), grp, "Integral", "LuaFunction#GridFunction#BoundarySubsets#InnerSubsets#Time");
		reg.add_function("IntegralNormalComponentOnManifold", static_cast<number (*)(const char*, SmartPtr<TFct>, const char*, number)>(&IntegralNormalComponentOnManifold<TFct>), grp, "Integral", "LuaFunction#GridFunction#BoundarySubsets#Time");
		reg.add_function("IntegralNormalComponentOnManifold", static_cast<number (*)(const char*, SmartPtr<TFct>, const char*, const char*)>(&IntegralNormalComponentOnManifold<TFct>), grp, "Integral", "LuaFunction#GridFunction#BoundarySubsets#InnerSubsets");
		reg.add_function("IntegralNormalComponentOnManifold", static_cast<number (*)(const char*, SmartPtr<TFct>, const char*)>(&IntegralNormalComponentOnManifold<TFct>), grp, "Integral", "LuaFunction#GridFunction#BoundarySubsets");
#endif
		{
			reg.add_function("IntegrateNormalComponentOnManifold", static_cast<number (*)(TFct&, const char*, const char*)>(&IntegrateNormalComponentOnManifold<TFct>),grp, "Integral", "GridFunction#Component#BoundarySubset");
	}
		}

//	IntegrateDiscFlux
	{
		reg.add_function("IntegrateDiscFlux", &IntegrateDiscFlux<TFct>, grp, "Integral");
	}

}

}; // end Functionality

// end group integrate_bridge
/// \}

}// end Integrate

/// \addtogroup integrate_bridge
void RegisterBridge_Integrate(Registry& reg, string grp)
{
	grp.append("/Discretization");
	using Functionality = Integrate::Functionality;

	{
		reg.add_function("TestQuadRule", &ug::TestQuadRule);
	}

	try{
		RegisterDomainAlgebraDependent<Functionality>(reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

}//	end of namespace bridge
}//	end of namespace ug
