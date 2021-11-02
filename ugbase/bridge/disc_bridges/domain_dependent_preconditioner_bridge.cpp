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

// preconditioner
#include "lib_algebra/lib_algebra.h"
#include "lib_disc/operator/preconditioner/line_smoothers.h"

using namespace std;

namespace ug{
namespace bridge{
namespace Preconditioner{

/**
 * \defgroup domdepprecond_bridge Domain Dependend Preconditioner Bridge
 * \ingroup precond_bridge
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

// end group domdepprecond_bridge
/// \}

}// end Preconditioner

/// \addtogroup domdepprecond_bridge
void RegisterBridge_DomainDependentPreconditioner(Registry& reg, string grp)
{
	grp.append("/Disc/Preconditioner");
	typedef Preconditioner::Functionality Functionality;

	try{		
		RegisterDomainAlgebraDependent<Functionality>(reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

} // namespace bridge
} // namespace ug
