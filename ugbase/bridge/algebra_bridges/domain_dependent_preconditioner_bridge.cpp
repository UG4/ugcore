
/*
 * preconditioner_bridge.cpp
 *
 *  Created on: 03.05.2012
 *      Author: avogel
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
#include "lib_algebra/operator/preconditioner/line_smoothers.h"

using namespace std;

namespace ug{
namespace bridge{
namespace Preconditioner{

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
		
//	Line Gauss-Seidel
	{
		typedef LineGaussSeidel<TDomain,TAlgebra> T;
		typedef IPreconditioner<TAlgebra> TBase;
		string name = string("LineGaussSeidel").append(suffix);
		reg.add_class_<T,TBase>(name, grp, "Line Gauss-Seidel Preconditioner")
		.template add_constructor<void (*)(SmartPtr<ApproximationSpace<TDomain> >)>("Approximation Space")
		.add_method("update", &T::update, "", "update")
		.add_method("set_num_steps", static_cast<void (T::*)(size_t,size_t,size_t,size_t,size_t,size_t)>(&T::set_num_steps), "", "set_num_steps")
		.add_method("set_num_steps", static_cast<void (T::*)(size_t,size_t,size_t,size_t)>(&T::set_num_steps), "", "set_num_steps")
		.add_method("set_num_steps", static_cast<void (T::*)(size_t,size_t)>(&T::set_num_steps), "", "set_num_steps")
		.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "LineGaussSeidel", tag);
	}
	
//	Line Vanka
	{
		typedef LineVanka<TDomain,TAlgebra> T;
		typedef IPreconditioner<TAlgebra> TBase;
		string name = string("LineVanka").append(suffix);
		reg.add_class_<T,TBase>(name, grp, "LineVanka Preconditioner")
		.template add_constructor<void (*)(SmartPtr<ApproximationSpace<TDomain> >)>("Approximation Space")
		.add_method("update", &T::update, "", "update")
		.add_method("set_num_steps", static_cast<void (T::*)(size_t,size_t,size_t,size_t,size_t,size_t)>(&T::set_num_steps), "", "set_num_steps")
		.add_method("set_num_steps", static_cast<void (T::*)(size_t,size_t,size_t,size_t)>(&T::set_num_steps), "", "set_num_steps")
		.add_method("set_num_steps", static_cast<void (T::*)(size_t,size_t)>(&T::set_num_steps), "", "set_num_steps")
		.add_method("set_relax", &T::set_relax, "", "relax")
		.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "LineVanka",  tag);
	}
	
};

}; // end Functionality
}// end Preconditioner


void RegisterBridge_DomainDependentPreconditioner(Registry& reg, string grp)
{
	grp.append("/Algebra/Preconditioner");
	typedef Preconditioner::Functionality Functionality;

	try{		
		RegisterDomainAlgebraDependent<Functionality>(reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

} // namespace bridge
} // namespace ug
