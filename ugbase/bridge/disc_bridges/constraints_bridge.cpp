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
#include "bridge/bridge.h"
#include "bridge/util.h"
#include "bridge/util_domain_algebra_dependent.h"

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
namespace Constraints{


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
	typedef typename TAlgebra::vector_type vector_type;
	typedef typename TAlgebra::matrix_type matrix_type;
	typedef ApproximationSpace<TDomain> approximation_space_type;

//	IDomainConstraint
	{
		typedef IConstraint<TAlgebra> TBase;
		typedef IDomainConstraint<TDomain, TAlgebra> T;
		string name = string("IDomainConstraint").append(suffix);
		reg.add_class_<T, TBase>(name, grp);
		reg.add_class_to_group(name, "IDomainConstraint", tag);
	}

//	OneSideP1Constraints
	{
		typedef OneSideP1Constraints<TDomain, TAlgebra> T;
		typedef IDomainConstraint<TDomain, TAlgebra> baseT;
		string name = string("OneSideP1Constraints").append(suffix);
		reg.add_class_<T, baseT>(name, grp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "OneSideP1Constraints", tag);
	}

//	SymP1Constraints
	{
		typedef SymP1Constraints<TDomain, TAlgebra> T;
		typedef IDomainConstraint<TDomain, TAlgebra> baseT;
		string name = string("SymP1Constraints").append(suffix);
		reg.add_class_<T, baseT>(name, grp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "SymP1Constraints", tag);
	}

//	DirichletBoundary
	{
		typedef DirichletBoundary<TDomain, TAlgebra> T;
		typedef IDomainConstraint<TDomain, TAlgebra> TBase;
		string name = string("DirichletBoundary").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.add_method("add", static_cast<void (T::*)(SmartPtr<UserData<number, dim, bool> >, const char*, const char*)>(&T::add),
						"", "Value#Function#Subsets")
			.add_method("add", static_cast<void (T::*)(SmartPtr<UserData<number, dim> >, const char*, const char*)>(&T::add),
						"", "Value#Function#Subsets")
			.add_method("add", static_cast<void (T::*)(SmartPtr<UserData<MathVector<dim>, dim> >, const char*, const char*)>(&T::add),
						"", "Vector#Functions#Subsets")
			.add_method("add",static_cast<void (T::*)(number, const char*, const char*)>(&T::add),
						"", "ConstantValue#Function#Subsets")
			.add_method("set_approximation_space",static_cast<void (T::*)(SmartPtr<ApproximationSpace<TDomain> >)>(&T::set_approximation_space),
						"", "ApproximationSpace")
#ifdef UG_FOR_LUA
			.add_method("add",static_cast<void (T::*)(const char*, const char*, const char*)>(&T::add),
						"", "LuaCallback#Function#Subsets")
#endif
			.add_method("clear", &T::clear)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "DirichletBoundary", tag);
	}
}

};
}// namespace Constraints

void RegisterBridge_Constraints(Registry& reg, string grp)
{
	grp.append("/Discretization/SpatialDisc");
	typedef Constraints::Functionality Functionality;

	try{
		RegisterDomainAlgebraDependent<Functionality>(reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

}//	end of namespace bridge
}//	end of namespace ug
