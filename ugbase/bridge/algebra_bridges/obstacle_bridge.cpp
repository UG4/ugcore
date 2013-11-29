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
#include "bridge/util_algebra_dependent.h"
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

template <typename TAlgebra>
static void Algebra(Registry& reg, string grp)
{
	string suffix = GetAlgebraSuffix<TAlgebra>();
	string tag = GetAlgebraTag<TAlgebra>();

//	typedefs for this algebra
	typedef typename TAlgebra::vector_type vector_type;
	typedef typename TAlgebra::matrix_type matrix_type;

	//	IObstacleConstraint
	{
		typedef IObstacleConstraint<TAlgebra> T;
		string name = string("IObstacleConstraint").append(suffix);
		reg.add_class_<T>(name, grp)
			.add_method("set_lower_obstacle", &T::set_lower_obstacle,
				"", "lower obstacle", "sets lower obstacle")
			.add_method("set_upper_obstacle", &T::set_upper_obstacle,
				"", "upper obstacle", "sets upper obstacle");
		reg.add_class_to_group(name, "IObstacleConstraint", tag);
	}

}

template <typename TDomain, typename TAlgebra>
static void DomainAlgebra(Registry& reg, string grp)
{
	string suffix = GetDomainAlgebraSuffix<TDomain,TAlgebra>();
	string tag = GetDomainAlgebraTag<TDomain,TAlgebra>();

	typedef GridFunction<TDomain, TAlgebra> function_type;

	//	Obstacle Classes for Projected Preconditioner

	//	ScalarObstacle
	{
		typedef ScalarObstacle<TDomain,TAlgebra> T;
		typedef IObstacleConstraint<TAlgebra> TBase;
		string name = string("ScalarObstacle").append(suffix);
		reg.add_class_<T,TBase>(name, grp)
			.template add_constructor<void (*)(const function_type&, const char*)>()
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ScalarObstacle", tag);
	}

	//	ObstacleInNormalDir
	{
		typedef ObstacleInNormalDir<TDomain,TAlgebra> T;
		typedef IObstacleConstraint<TAlgebra> TBase;
		string name = string("ObstacleInNormalDir").append(suffix);
		reg.add_class_<T,TBase>(name, grp)
			.template add_constructor<void (*)(const function_type&, const char*)>()
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ObstacleInNormalDir", tag);
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
		RegisterAlgebraDependent<Functionality>(reg,grp);
		RegisterDomainAlgebraDependent<Functionality>(reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

} // namespace bridge
} // namespace ug
