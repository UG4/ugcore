/*
 * domain_bridge.cpp
 *
 *  Created on: 22.09.2010
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
#include "lib_algebra/operator/operator_util.h"
#include "lib_algebra/operator/interface/operator.h"
#include "lib_algebra/operator/interface/operator_inverse.h"

// lib_disc includes
#include "lib_disc/domain.h"
#include "lib_disc/spatial_disc/domain_disc.h"
#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/function_spaces/approximation_space.h"
#include "lib_disc/function_spaces/grid_function_util.h"
#include "lib_disc/function_spaces/interpolate.h"
#include "lib_disc/function_spaces/integrate.h"
#include "lib_disc/function_spaces/integrateDraft.h"
#include "lib_disc/dof_manager/cuthill_mckee.h"
#include "lib_disc/dof_manager/lexorder.h"

using namespace std;

namespace ug{
namespace bridge{

template <typename TDomain, typename TAlgebra>
static void Register__Algebra_Domain(Registry& reg, string parentGroup)
{
//	typedef
	static const int dim = TDomain::dim;
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

//	DomainDiscretization
	{
		typedef IDomainDiscretization<TAlgebra> TBase;
		typedef DomainDiscretization<TDomain, TAlgebra> T;
		string name = string("DomainDiscretization").append(dimAlgSuffix);
		reg.add_class_<T, TBase>(name, domDiscGrp)
			.template add_constructor<void (*)(SmartPtr<ApproximationSpace<TDomain> >)>("ApproximationSpace")
			.add_method("add", static_cast<void (T::*)(SmartPtr<IDomainConstraint<TDomain, TAlgebra> >)>(&T::add), "", "Post Process")
			.add_method("add", static_cast<void (T::*)(SmartPtr<IDomainElemDisc<TDomain> >)>(&T::add), "", "Element Discretization")
			.add_method("add", static_cast<void (T::*)(SmartPtr<IDiscretizationItem<TDomain, TAlgebra> >)>(&T::add), "", "DiscItem")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "DomainDiscretization", dimAlgTag);
	}

//	IDiscretizationItem
	{
		typedef IDiscretizationItem<TDomain, TAlgebra> T;
		string name = string("IDiscretizationItem").append(dimAlgSuffix);
		reg.add_class_<T>(name, domDiscGrp);
		reg.add_class_to_group(name, "IDiscretizationItem", dimAlgTag);
	}

//	GridFunction
	{
		string name = string("GridFunction").append(dimAlgSuffix);
		reg.add_class_<TFct, vector_type>(name, approxGrp)
			.template add_constructor<void (*)(SmartPtr<approximation_space_type>)>("ApproximationSpace")
			.add_method("assign", static_cast<void (TFct::*)(const vector_type&)>(&TFct::assign),
						"Success", "Vector")
			.add_method("clone", &TFct::clone)
			.add_method("add_transfer", &TFct::add_transfer)
			.add_method("remove_transfer", &TFct::remove_transfer)
			.add_method("clear_transfers", &TFct::clear_transfers)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "GridFunction", dimAlgTag);
	}

	string grp = parentGroup; grp.append("");

//	InterpolateFunction
	{
		reg.add_function("InterpolateFunction", static_cast<void (*)(IPData<number, dim>&, TFct&, const char*, number)>(&InterpolateFunction<TFct>), grp);
		reg.add_function("InterpolateFunction",static_cast<void (*)(IPData<number, dim>&, TFct&, const char*, number, const char*)>(&InterpolateFunction<TFct>),grp);
#ifdef UG_FOR_LUA
		reg.add_function("InterpolateFunction", static_cast<void (*)(const char*, TFct&, const char*, number)>(&InterpolateFunction<TFct>), grp);
		reg.add_function("InterpolateFunction",static_cast<void (*)(const char*, TFct&, const char*, number, const char*)>(&InterpolateFunction<TFct>),grp);
#endif
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
		UG_LOG("### ERROR in RegisterLibDisc_Domain: "
				"Registration failed (using name " << ex.name << ").\n");
		UG_THROW_FATAL("Registration failed.");
	}
}

template <typename TDomain>
static void Register__Domain(Registry& reg, string parentGroup)
{
//	typedef
//	static const int dim = TDomain::dim;
	typedef ApproximationSpace<TDomain> approximation_space_type;

//	group string
	string approxGrp = parentGroup; approxGrp.append("/ApproximationSpace");

//	suffix and tag
	string dimSuffix = GetDomainSuffix<TDomain>();
	string dimTag = GetDomainTag<TDomain>();


//  ApproximationSpace
	{
		typedef ApproximationSpace<TDomain> T;
		typedef IApproximationSpace TBase;
		string name = string("ApproximationSpace").append(dimSuffix);
		reg.add_class_<T, TBase>(name, approxGrp)
			.template add_constructor<void (*)(SmartPtr<TDomain>)>("Domain")
			.add_method("domain", static_cast<SmartPtr<TDomain> (T::*)()>(&T::domain))
			.add_method("surface_view", static_cast<ConstSmartPtr<SurfaceView> (T::*)() const>(&T::surface_view))
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ApproximationSpace", dimTag);
	}

//	Order Cuthill-McKee
	{
		reg.add_function("OrderCuthillMcKee", static_cast<void (*)(approximation_space_type&, bool)>(&OrderCuthillMcKee), approxGrp);
	}

//	Order lexicographically
	{
		reg.add_function("OrderLex", static_cast<void (*)(approximation_space_type&, const char*)>(&OrderLex<TDomain>), approxGrp);
	}
}

bool RegisterLibDisc_Domain(Registry& reg, string parentGroup)
{
//	GridLevel
	reg.add_class_<GridLevel>("GridLevel")
		.add_constructor()
		.add_constructor<void (*)(int)>("Level")
		.add_constructor<void (*)(int, std::string)>("Level, Type")
		.set_construct_as_smart_pointer(true);

//	IApproximationSpace
	{
	typedef IApproximationSpace T;
	reg.add_class_<T>("IApproximationSpace")
		.add_method("print_statistic", static_cast<void (T::*)(int) const>(&T::print_statistic))
		.add_method("print_statistic", static_cast<void (T::*)() const>(&T::print_statistic))
		.add_method("print_layout_statistic", static_cast<void (T::*)(int) const>(&T::print_layout_statistic))
		.add_method("print_layout_statistic", static_cast<void (T::*)() const>(&T::print_layout_statistic))
		.add_method("print_local_dof_statistic", static_cast<void (T::*)(int) const>(&T::print_local_dof_statistic))
		.add_method("print_local_dof_statistic", static_cast<void (T::*)() const>(&T::print_local_dof_statistic))
		.add_method("num_levels", &T::num_levels)
		.add_method("init_levels", &T::init_levels)
		.add_method("init_surfaces", &T::init_surfaces)
		.add_method("init_top_surface", &T::init_top_surface)
		.add_method("defragment", &T::defragment)
		.add_method("clear", &T::clear)
		.add_method("add_fct", static_cast<void (T::*)(const char*, const char*, int, const char*)>(&T::add_fct),
					"", "Name # Type|selection|value=[\"Lagrange\",\"DG\"] # Order # Subsets", "Adds a function to the Function Pattern",
					"currently no help available")
		.add_method("add_fct", static_cast<void (T::*)(const char*, const char*, int)>(&T::add_fct),
					"", "Name # Type|selection|value=[\"Lagrange\",\"DG\"] # Order", "Adds a function to the Function Pattern",
					"currently no help available");

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
		UG_LOG("### ERROR in RegisterLibDisc_Domain: "
				"Registration failed (using name " << ex.name << ").\n");
		UG_THROW_FATAL("Registration failed.");
	}

	return true;
}

}//	end of namespace ug
}//	end of namespace interface
