/*
 * domain_bridge.cpp
 *
 *  Created on: 21.05.2012
 *      Author: andreasvogel
 */

// extern headers
#include <iostream>
#include <sstream>
#include <string>

// include bridge
#include "../bridge.h"
#include "registry/registry.h"

// lib_algebra includes
#include "lib_algebra/cpu_algebra_types.h"

// lib_disc includes
#include "lib_disc/domain.h"
#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/dof_manager/surface_dof_distribution.h"
#include "lib_disc/function_spaces/integrate.h"
#include "lib_disc/function_spaces/integrateDraft.h"

using namespace std;

namespace ug{
namespace bridge{

template <typename TDomain, typename TAlgebra>
static void Register__Algebra_Domain(Registry& reg, string parentGroup)
{
//	typedef
	static const int dim = TDomain::dim;
	typedef GridFunction<TDomain, SurfaceDoFDistribution, TAlgebra> TFct;

//	group string
	string grp = parentGroup; grp.append("");

//	Integral
	{
		reg.add_function("Integral", &Integral<TFct>, grp);
	}

//	L2Error
	{
		reg.add_function("L2Error", static_cast<number (*)(IPData<number, dim>&, TFct&, const char*, number)>(&L2Error<TFct>), grp);
		reg.add_function("L2Error",static_cast<number (*)(IPData<number, dim>&, TFct&, const char*, number, const char*)>(&L2Error<TFct>),grp);
#ifdef UG_FOR_LUA
		reg.add_function("L2Error", static_cast<number (*)(const char*, TFct&, const char*, number)>(&L2Error<TFct>), grp);
		reg.add_function("L2Error",static_cast<number (*)(const char*, TFct&, const char*, number, const char*)>(&L2Error<TFct>),grp);
#endif
	}

//	L2ErrorDraft
	{
		reg.add_function("L2ErrorDraft",static_cast<number (*)(IPData<number, dim>&, TFct&, const char*, number, int, const char*)>(&L2ErrorDraft<TFct>), grp);
#ifdef UG_FOR_LUA
		reg.add_function("L2ErrorDraft",static_cast<number (*)(const char*, TFct&, const char*, number, int, const char*)>(&L2ErrorDraft<TFct>), grp);
#endif
	}

//	L2Norm
	{
		typedef number (*fct_type)(TFct&, const char*, int, const char*);
		reg.add_function("L2Norm",static_cast<fct_type>(&L2Norm<TFct>),grp);
	}

//	StdFuncIntegral
	{
		typedef number (*fct_type)(TFct&, const char*, int, const char*);
		reg.add_function("StdFuncIntegral",static_cast<fct_type>(&StdFuncIntegral<TFct>),grp);
	}

//	IntegrateFluxOnBoundary
	{
		typedef number (*fct_type)(TFct&, const char*, const char*, const char*);
		reg.add_function("IntegrateFluxOnBoundary",static_cast<fct_type>(&IntegrateFluxOnBoundary<TFct>),grp,
		                 "Integral", "GridFunction#Component#BoundarySubset#InnerSubset");
	}



}


template <typename TAlgebra>
static void Register__Algebra(Registry& reg, string parentGroup)
{
//	get group string
	string grp = parentGroup; grp.append("/Discretization");

	try{
#ifdef UG_DIM_1
		Register__Algebra_Domain<Domain1d, TAlgebra>(reg, grp);
#endif
#ifdef UG_DIM_2
		Register__Algebra_Domain<Domain2d, TAlgebra>(reg, grp);
#endif
#ifdef UG_DIM_3
		Register__Algebra_Domain<Domain3d, TAlgebra>(reg, grp);
#endif
	}
	catch(UG_REGISTRY_ERROR_RegistrationFailed ex)
	{
		UG_LOG("### ERROR in RegisterIntegrate: "
				"Registration failed (using name " << ex.name << ").\n");
		UG_THROW("Registration failed.");
	}
}

void RegisterBridge_Integrate(Registry& reg, string parentGroup)
{
	try{
#ifdef UG_CPU_1
	Register__Algebra<CPUAlgebra>(reg, parentGroup);
#endif
#ifdef UG_CPU_2
	Register__Algebra<CPUBlockAlgebra<2> >(reg, parentGroup);
#endif
#ifdef UG_CPU_3
	Register__Algebra<CPUBlockAlgebra<3> >(reg, parentGroup);
#endif
#ifdef UG_CPU_4
	Register__Algebra<CPUBlockAlgebra<4> >(reg, parentGroup);
#endif
#ifdef UG_CPU_VAR
	Register__Algebra<CPUVariableBlockAlgebra >(reg, parentGroup);
#endif
	}
	catch(UG_REGISTRY_ERROR_RegistrationFailed ex)
	{
		UG_LOG("### ERROR in RegisterIntegrate: "
				"Registration failed (using name " << ex.name << ").\n");
		UG_THROW("Registration failed.");
	}
}

}//	end of namespace ug
}//	end of namespace interface
