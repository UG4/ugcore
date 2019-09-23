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
 * \defgroup constraints_bridge Constraints Bridge
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

//	IDomainConstraint
	{
		typedef IConstraint<TAlgebra> TBase;
		typedef IDomainConstraint<TDomain, TAlgebra> T;
		string name = string("IDomainConstraint").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_method("set_error_estimator", &T::set_error_estimator, "", "error estimator data object");
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
			.template add_constructor<void (*)()>()
			.template add_constructor<void (*)(bool)>()
#ifdef LAGRANGE_DIRICHLET_ADJ_TRANSFER_FIX
			.template add_constructor<void (*)(bool,bool)>()
#endif
			.add_method("add", static_cast<void (T::*)(SmartPtr<UserData<number, dim, bool> >, const char*, const char*)>(&T::add),
						"", "Value#Function#Subsets")
			.add_method("add", static_cast<void (T::*)(SmartPtr<UserData<number, dim, bool> >, const std::vector<std::string>&, const std::vector<std::string>&)>(&T::add),
						"", "Value#Function#Subsets")
			.add_method("add", static_cast<void (T::*)(SmartPtr<UserData<number, dim> >, const char*, const char*)>(&T::add),
						"", "Value#Function#Subsets")
			.add_method("add", static_cast<void (T::*)(SmartPtr<UserData<number, dim> >, const std::vector<std::string>&, const std::vector<std::string>&)>(&T::add),
						"", "Value#Function#Subsets")
			.add_method("add", static_cast<void (T::*)(SmartPtr<UserData<MathVector<dim>, dim> >, const char*, const char*)>(&T::add),
						"", "Vector#Functions#Subsets")
			.add_method("add", static_cast<void (T::*)(SmartPtr<UserData<MathVector<dim>, dim> >, const std::vector<std::string>&, const std::vector<std::string>&)>(&T::add),
						"", "Vector#Functions#Subsets")
			.add_method("add",static_cast<void (T::*)(number, const char*, const char*)>(&T::add),
						"", "ConstantValue#Function#Subsets")
			.add_method("add",static_cast<void (T::*)(number, const std::vector<std::string>&, const std::vector<std::string>&)>(&T::add),
						"", "ConstantValue#Function#Subsets")
			.add_method("add",static_cast<void (T::*)(const char*, const char*)>(&T::add),
						"", "Function#Subsets")
			.add_method("add",static_cast<void (T::*)(const std::vector<std::string>&, const std::vector<std::string>&)>(&T::add),
						"", "Function#Subsets")
			.add_method("invert_subset_selection",static_cast<void (T::*)()>(&T::invert_subset_selection),
						"", "")
			.add_method("set_approximation_space",static_cast<void (T::*)(SmartPtr<ApproximationSpace<TDomain> >)>(&T::set_approximation_space),
						"", "ApproximationSpace")
#ifdef UG_FOR_LUA
			.add_method("add",static_cast<void (T::*)(const char*, const char*, const char*)>(&T::add),
						"", "LuaCallback#Function#Subsets")
			.add_method("add",static_cast<void (T::*)(const char*, const std::vector<std::string>&, const std::vector<std::string>&)>(&T::add),
						"", "LuaCallback#Function#Subsets")
			.add_method("add",static_cast<void (T::*)(LuaFunctionHandle, const char*, const char*)>(&T::add),
						"", "LuaCallback#Function#Subsets")
			.add_method("add",static_cast<void (T::*)(LuaFunctionHandle, const std::vector<std::string>&, const std::vector<std::string>&)>(&T::add),
						"", "LuaCallback#Function#Subsets")
#endif
			.add_method("clear", &T::clear)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "DirichletBoundary", tag);
	}
}

};

// end group constraints_bridge
/// \}

}// namespace Constraints

/// \addtogroup constraints_bridge
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
