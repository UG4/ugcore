/*
 * obstacle_bridge.cpp
 *
 *  Created on: 26.11.2013
 *      Author: raphaelprohl
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

template <typename TDomain, typename TAlgebra>
static void DomainAlgebra(Registry& reg, string grp)
{
	string suffix = GetDomainAlgebraSuffix<TDomain,TAlgebra>();
	string tag = GetDomainAlgebraTag<TDomain,TAlgebra>();

	typedef GridFunction<TDomain, TAlgebra> function_type;
	static const int dim = TDomain::dim;

//	Obstacle Classes

	//	IObstacleConstraint
	{
		typedef IObstacleConstraint<TDomain,TAlgebra> T;
		string name = string("IObstacleConstraint").append(suffix);
		reg.add_class_<T>(name, grp)
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
		reg.add_class_<T,TBase>(name, grp)
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
		reg.add_class_<T,TBase>(name, grp)
			.template add_constructor<void (*)(const function_type&)>()
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ScalarUpperObstacle", tag);
	}

	//	ObstacleInNormalDir
	/*{
		typedef ObstacleInNormalDir<TDomain,TAlgebra> T;
		typedef IObstacleConstraint<TAlgebra> TBase;
		string name = string("ObstacleInNormalDir").append(suffix);
		reg.add_class_<T,TBase>(name, grp)
			.template add_constructor<void (*)(const function_type&, const char*)>()
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ObstacleInNormalDir", tag);
	}*/


//	Projected Preconditioners

	//	IProjGaussSeidel
	{
		typedef IProjGaussSeidel<TDomain, TAlgebra> T;
		typedef GaussSeidelBase<TAlgebra> TBase;
		string name = string("IProjGaussSeidel").append(suffix);
		reg.add_class_<T,TBase>(name, grp)
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
		reg.add_class_<T,TBase>(name, grp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ProjGaussSeidel", tag);
	}

	//	ProjBackwardGaussSeidel
	{
		typedef ProjBackwardGaussSeidel<TDomain,TAlgebra> T;
		typedef IProjGaussSeidel<TDomain,TAlgebra> TBase;
		string name = string("ProjBackwardGaussSeidel").append(suffix);
		reg.add_class_<T,TBase>(name, grp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ProjBackwardGaussSeidel", tag);
	}

	//	ProjSymmetricGaussSeidel
	{
		typedef ProjSymmetricGaussSeidel<TDomain,TAlgebra> T;
		typedef IProjGaussSeidel<TDomain,TAlgebra> TBase;
		string name = string("ProjSymmetricGaussSeidel").append(suffix);
		reg.add_class_<T,TBase>(name, grp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ProjSymmetricGaussSeidel", tag);
	}
};

}; // end Functionality

// end group obstacle_bridge
/// \}

}// end Obstacle

/// \addtogroup obstacle_bridge
void RegisterBridge_Obstacle(Registry& reg, string grp)
{
	grp.append("/Algebra/Obstacle");
	typedef Obstacle::Functionality Functionality;

	try{
		RegisterDomainAlgebraDependent<Functionality>(reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

} // namespace bridge
} // namespace ug
