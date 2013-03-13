/*
 * integrate_bridge.cpp
 *
 *  Created on: 21.05.2012
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
#include "lib_disc/dof_manager/surface_dof_distribution.h"
#include "lib_disc/function_spaces/integrate.h"
#include "lib_disc/function_spaces/integrate_flux.h"

#include "lib_disc/quadrature/quad_test.h"

using namespace std;

namespace ug{
namespace bridge{
namespace Integrate{

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
//	typedef
	static const int dim = TDomain::dim;
	typedef GridFunction<TDomain, TAlgebra> TFct;

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
	}

//	IntegrateDiscFlux
	{
		reg.add_function("IntegrateDiscFlux", &IntegrateDiscFlux<TFct>, grp, "Integral");
	}

}

}; // end Functionality
}// end Integrate

void RegisterBridge_Integrate(Registry& reg, string grp)
{
	grp.append("/Discretization");
	typedef Integrate::Functionality Functionality;


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
