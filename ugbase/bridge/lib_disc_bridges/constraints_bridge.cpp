/*
 * constraints_bridge.cpp
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

#include "lib_disc/spatial_disc/domain_disc.h"
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"
#include "lib_disc/spatial_disc/constraints/constraint_interface.h"
#include "lib_disc/spatial_disc/constraints/dirichlet_boundary/lagrange_dirichlet_boundary.h"
#include "lib_disc/spatial_disc/constraints/continuity_constraints/p1_continuity_constraints.h"

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

//	suffix and tag
	string dimAlgSuffix = GetDomainSuffix<TDomain>();
	dimAlgSuffix.append(GetAlgebraSuffix<TAlgebra>());

	string dimAlgTag = GetDomainTag<TDomain>();
	dimAlgTag.append(GetAlgebraTag<TAlgebra>());

//	group string
	string approxGrp = parentGroup; approxGrp.append("/ApproximationSpace");
	string domDiscGrp = parentGroup; domDiscGrp.append("/SpatialDisc");

//	IDomainConstraint
	{
		std::string grp = parentGroup; grp.append("/Discretization/SpatialDisc");
		typedef IConstraint<TAlgebra> TBase;
		typedef IDomainConstraint<TDomain, TAlgebra> T;
		string name = string("IDomainConstraint").append(dimAlgSuffix);
		reg.add_class_<T, TBase>(name, grp);
		reg.add_class_to_group(name, "IDomainConstraint", dimAlgTag);
	}

//	OneSideP1Constraints
	{
		std::string grp = parentGroup; grp.append("/Discretization/SpatialDisc");
		typedef OneSideP1Constraints<TDomain, TAlgebra> T;
		typedef IDomainConstraint<TDomain, TAlgebra> baseT;
		string name = string("OneSideP1Constraints").append(dimAlgSuffix);
		reg.add_class_<T, baseT>(name, grp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "OneSideP1Constraints", dimAlgTag);
	}

//	SymP1Constraints
	{
		std::string grp = parentGroup; grp.append("/Discretization/SpatialDisc");
		typedef SymP1Constraints<TDomain, TAlgebra> T;
		typedef IDomainConstraint<TDomain, TAlgebra> baseT;
		string name = string("SymP1Constraints").append(dimAlgSuffix);
		reg.add_class_<T, baseT>(name, grp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "SymP1Constraints", dimAlgTag);
	}

//	DirichletBoundary
	{
		typedef DirichletBoundary<TDomain, TAlgebra> T;
		typedef IDomainConstraint<TDomain, TAlgebra> TBase;
		string name = string("DirichletBoundary").append(dimAlgSuffix);
		reg.add_class_<T, TBase>(name, domDiscGrp)
			.add_constructor()
			.add_method("add", static_cast<void (T::*)(SmartPtr<IPData<number, dim, bool> >, const char*, const char*)>(&T::add),
						"", "Value#Function#Subsets")
			.add_method("add", static_cast<void (T::*)(SmartPtr<IPData<number, dim> >, const char*, const char*)>(&T::add),
						"", "Value#Function#Subsets")
			.add_method("add", static_cast<void (T::*)(SmartPtr<IPData<MathVector<dim>, dim> >, const char*, const char*)>(&T::add),
						"", "Vector#Functions#Subsets")
			.add_method("add",static_cast<void (T::*)(number, const char*, const char*)>(&T::add),
						"", "ConstantValue#Function#Subsets")
#ifdef UG_FOR_LUA
			.add_method("add",static_cast<void (T::*)(const char*, const char*, const char*)>(&T::add),
						"", "LuaCallback#Function#Subsets")
#endif
			.add_method("clear", &T::clear)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "DirichletBoundary", dimAlgTag);
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
		UG_LOG("### ERROR in RegisterConstraints: "
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

bool RegisterConstraints(Registry& reg, string parentGroup)
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
		UG_LOG("### ERROR in RegisterConstraints: "
				"Registration failed (using name " << ex.name << ").\n");
		UG_THROW_FATAL("Registration failed.");
	}

	return true;
}

}//	end of namespace ug
}//	end of namespace interface
