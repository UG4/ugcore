/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Raphael Prohl
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

// preconditioner
#include "lib_algebra/lib_algebra.h"
#include "lib_algebra/operator/preconditioner/projected_gauss_seidel/obstacles/obstacles.h"
#include "lib_disc/operator/non_linear_operator/truncated_monotone_mg/truncated_monotone_transfer.h"

using namespace std;

namespace ug{
namespace bridge{
namespace Obstacle{

/**
 * \defgroup obstacle_bridge Obstacle Bridge
 * \ingroup algebra_bridge
 * \{
 */

/**
 * Class exporting the functionality. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */


struct Functionality
{

template <typename TDomain, typename TAlgebra, typename TRegistry>
static void DomainAlgebra(TRegistry& reg, string grp)
{
	string suffix = GetDomainAlgebraSuffix<TDomain,TAlgebra>();
	string tag = GetDomainAlgebraTag<TDomain,TAlgebra>();

	typedef GridFunction<TDomain, TAlgebra> function_type;
	static const int dim = TDomain::dim;

//	Obstacle Classes

	//	IObstacleConstraint
	{
		typedef IObstacleConstraint<TDomain,TAlgebra> T;
		typedef IDomainConstraint<TDomain, TAlgebra> TBase;
		string name = string("IObstacleConstraint").append(suffix);
		reg.template add_class_<T,TBase>(name, grp)
			.add_method("add", static_cast<void (T::*)(SmartPtr<UserData<number, dim, bool> >, const char*)>(&T::add),
						"", "Value#Function")
			.add_method("add", static_cast<void (T::*)(SmartPtr<UserData<number, dim, bool> >, const char*, const char*)>(&T::add),
						"", "Value#Function#Subsets")
			.add_method("add", static_cast<void (T::*)(SmartPtr<UserData<number, dim> >, const char*)>(&T::add),
							"", "Value#Function")
			.add_method("add", static_cast<void (T::*)(SmartPtr<UserData<number, dim> >, const char*, const char*)>(&T::add),
						"", "Value#Function#Subsets")
			.add_method("add", static_cast<void (T::*)(SmartPtr<UserData<MathVector<dim>, dim> >, const char*)>(&T::add),
						"", "Vector#Functions")
			.add_method("add", static_cast<void (T::*)(SmartPtr<UserData<MathVector<dim>, dim> >, const char*, const char*)>(&T::add),
						"", "Vector#Functions#Subsets")
			.add_method("add",static_cast<void (T::*)(number, const char*)>(&T::add),
						"", "ConstantValue#Function")
			.add_method("add",static_cast<void (T::*)(number, const char*, const char*)>(&T::add),
						"", "ConstantValue#Function#Subsets")
#ifdef UG_FOR_LUA
			.add_method("add",static_cast<void (T::*)(const char*, const char*)>(&T::add),
						"", "LuaCallback#Function")
			.add_method("add",static_cast<void (T::*)(const char*, const char*, const char*)>(&T::add),
						"", "LuaCallback#Function#Subsets")
#endif
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "IObstacleConstraint", tag);
	}

	//	ScalarLowerObstacle
	{
		typedef ScalarLowerObstacle<TDomain,TAlgebra> T;
		typedef IObstacleConstraint<TDomain,TAlgebra> TBase;
		string name = string("ScalarLowerObstacle").append(suffix);
		reg.template add_class_<T,TBase>(name, grp)
			.template add_constructor<void (*)(const function_type&)>()
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ScalarLowerObstacle", tag);
	}

	//	ScalarUpperObstacle
	{
		typedef ScalarUpperObstacle<TDomain,TAlgebra> T;
		typedef IObstacleConstraint<TDomain,TAlgebra> TBase;
		string name = string("ScalarUpperObstacle").append(suffix);
		reg.template add_class_<T,TBase>(name, grp)
			.template add_constructor<void (*)(const function_type&)>()
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ScalarUpperObstacle", tag);
	}

	//	ObstacleInNormalDir
	{
		typedef ObstacleInNormalDir<TDomain,TAlgebra> T;
		typedef IObstacleConstraint<TDomain,TAlgebra> TBase;
		string name = string("ObstacleInNormalDir").append(suffix);
		reg.template add_class_<T,TBase>(name, grp)
			.template add_constructor<void (*)(const function_type&)>()
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ObstacleInNormalDir", tag);
	}


//	Projected Preconditioners

	//	IProjGaussSeidel
	{
		typedef IProjGaussSeidel<TDomain, TAlgebra> T;
		typedef GaussSeidelBase<TAlgebra> TBase;
		string name = string("IProjGaussSeidel").append(suffix);
		reg.template add_class_<T,TBase>(name, grp)
			.add_method("add_obstacle_constraint", &T::add_obstacle_constraint,
				"", "obstacle constraint", "adds an obstacle constraint")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "IProjGaussSeidel", tag);
	}

	//	ProjGaussSeidel
	{
		typedef ProjGaussSeidel<TDomain,TAlgebra> T;
		typedef IProjGaussSeidel<TDomain,TAlgebra> TBase;
		string name = string("ProjGaussSeidel").append(suffix);
		reg.template add_class_<T,TBase>(name, grp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ProjGaussSeidel", tag);
	}

	//	ProjBackwardGaussSeidel
	{
		typedef ProjBackwardGaussSeidel<TDomain,TAlgebra> T;
		typedef IProjGaussSeidel<TDomain,TAlgebra> TBase;
		string name = string("ProjBackwardGaussSeidel").append(suffix);
		reg.template add_class_<T,TBase>(name, grp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ProjBackwardGaussSeidel", tag);
	}

	//	ProjSymmetricGaussSeidel
	{
		typedef ProjSymmetricGaussSeidel<TDomain,TAlgebra> T;
		typedef IProjGaussSeidel<TDomain,TAlgebra> TBase;
		string name = string("ProjSymmetricGaussSeidel").append(suffix);
		reg.template add_class_<T,TBase>(name, grp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ProjSymmetricGaussSeidel", tag);
	}


//	transfer operators for truncated monotone multigrid
	{
		typedef TruncatedMonotoneTransfer<TDomain,TAlgebra> T;
		typedef StdTransfer<TDomain,TAlgebra> TBase;
		string name = string("TruncatedMonotoneTransfer").append(suffix);
		reg.template add_class_<T,TBase>(name, grp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "TruncatedMonotoneTransfer", tag);
	}
};

}; // end Functionality

// end group obstacle_bridge
/// \}

}// end Obstacle

/// \addtogroup obstacle_bridge
template <typename TRegistry=Registry>
void RegisterBridge_Obstacle_(TRegistry& reg, string grp)
{
	grp.append("/Algebra/Obstacle");
	typedef Obstacle::Functionality Functionality;

	try{
		RegisterDomainAlgebraDependent<Functionality, TRegistry>(reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

} // namespace bridge

UG_REGISTRY_DEFINE(RegisterBridge_Obstacle);
} // namespace ug
