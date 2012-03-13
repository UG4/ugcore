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
#include "../bridge.h"
#include "registry/registry.h"

// lib_algebra includes
#include "lib_algebra/cpu_algebra_types.h"

// lib_disc includes
#include "lib_disc/domain.h"
#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/function_spaces/approximation_space.h"
#include "lib_disc/function_spaces/error_indicator.h"
#include "lib_disc/function_spaces/level_transfer.h"
#include "lib_disc/function_spaces/local_transfer.h"


using namespace std;

namespace ug{
namespace bridge{

template <typename TDomain, typename TAlgebra>
static void Register__Algebra_Domain(Registry& reg, string parentGroup)
{
//	typedef
	typedef typename TAlgebra::vector_type vector_type;
	typedef typename TAlgebra::matrix_type matrix_type;
	typedef ApproximationSpace<TDomain> approximation_space_type;
	typedef GridFunction<TDomain, SurfaceDoFDistribution, TAlgebra> TFct;

//	suffix and tag
	string dimAlgSuffix = GetDomainSuffix<TDomain>();
	dimAlgSuffix.append(GetAlgebraSuffix<TAlgebra>());

	string dimAlgTag = GetDomainTag<TDomain>();
	dimAlgTag.append(GetAlgebraTag<TAlgebra>());

//	group string
	string approxGrp = parentGroup; approxGrp.append("/ApproximationSpace");
	string domDiscGrp = parentGroup; domDiscGrp.append("/SpatialDisc");

//	MarkForRefinement_GradientIndicator
	{
		string grp("ug4/Refinement/");
		reg.add_function("MarkForRefinement_GradientIndicator",
						 &MarkForRefinement_GradientIndicator<TDomain, SurfaceDoFDistribution, TAlgebra>, grp);
	}

//	Prolongate
	{
		string grp("ug4/");
		reg.add_function("Prolongate",
						 &Prolongate<TDomain, SurfaceDoFDistribution, TAlgebra>, grp);
	}
}


template <typename TAlgebra>
static void Register__Algebra(Registry& reg, string parentGroup)
{
//	get group string
	string grp = parentGroup; grp.append("/Discretization");

//	suffix and tag
	string algSuffix = GetAlgebraSuffix<TAlgebra>();
	string algTag = GetAlgebraTag<TAlgebra>();

//	ILocalTransferAlgebra
	{
		typedef ILocalTransferAlgebra<TAlgebra> T;
		typedef ILocalTransfer TBase;
		string name = string("ILocalTransferAlgebra").append(algSuffix);
		reg.add_class_<T, TBase>(name)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ILocalTransferAlgebra", algTag);
	}

//	P1LocalTransfer
	{
		typedef P1LocalTransfer<TAlgebra> T;
		typedef ILocalTransferAlgebra<TAlgebra> TBase;
		string name = string("P1LocalTransfer").append(algSuffix);
		reg.add_class_<T, TBase>(name)
			.template add_constructor<void (*)(size_t)>("fct")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "P1LocalTransfer", algTag);
	}

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
		UG_LOG("### ERROR in RegisterAdaptiveTools: "
				"Registration failed (using name " << ex.name << ").\n");
		UG_THROW_FATAL("Registration failed.");
	}
}

template <typename TDomain>
static void Register__Domain(Registry& reg, string parentGroup)
{
//	suffix and tag
	string dimSuffix = GetDomainSuffix<TDomain>();
	string dimTag = GetDomainTag<TDomain>();
}

bool RegisterAdaptiveTools(Registry& reg, string parentGroup)
{
//	ILocalTransfer
	{
		reg.add_class_<ILocalTransfer>("ILocalTransfer");
	}

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

#ifdef UG_DIM_1
	Register__Domain<Domain1d>(reg, parentGroup);
#endif
#ifdef UG_DIM_2
	Register__Domain<Domain2d>(reg, parentGroup);
#endif
#ifdef UG_DIM_3
	Register__Domain<Domain3d>(reg, parentGroup);
#endif
	}
	catch(UG_REGISTRY_ERROR_RegistrationFailed ex)
	{
		UG_LOG("### ERROR in RegisterAdaptiveTools: "
				"Registration failed (using name " << ex.name << ").\n");
		UG_THROW_FATAL("Registration failed.");
	}

	return true;
}

}//	end of namespace ug
}//	end of namespace interface
