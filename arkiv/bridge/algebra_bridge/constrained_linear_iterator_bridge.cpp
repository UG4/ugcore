/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
 * Author: Markus Breit
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

// linear iterators
//#include "lib_algebra/lib_algebra.h"
#include "lib_algebra/operator/interface/constrained_linear_iterator.h"
#include "lib_algebra/operator/preconditioner/jacobi.h"
#include "lib_algebra/operator/preconditioner/ilu.h"
#include "lib_algebra/operator/preconditioner/ilut.h"
#include "lib_algebra/operator/linear_solver/lu.h"
#include "lib_disc/operator/linear_operator/multi_grid_solver/mg_solver.h"

using namespace std;

namespace ug{
namespace bridge{

/**
 * \defgroup constr_lin_it_bridge Constrained linear iterator bridge
 * \ingroup algebra_bridge
 * \{
 */

/**
 * Class exporting the functionality. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */
struct CLI_Functionality
{

	
template <typename TDomain, typename TAlgebra>
static void DomainAlgebra(Registry& reg, string grp)
{
	string suffix = GetDomainAlgebraSuffix<TDomain,TAlgebra>();
	string tag = GetDomainAlgebraTag<TDomain,TAlgebra>();

	string precondGrp(grp);
	precondGrp.append("/Preconditioner");
	string solverGrp(grp);
	solverGrp.append("/Solver");
		
	// constrained Jacobi
	{
		using T = ConstrainedLinearIterator<TDomain, TAlgebra, Jacobi<TAlgebra> >;
		using TBase = Jacobi<TAlgebra>;
		string name = string("Jacobi_c").append(suffix);
		reg.add_class_<T,TBase>(name, precondGrp, "Jacobi preconditioner respecting constraints")
		.template add_constructor<void (*)(SmartPtr<IDomainDiscretization<TAlgebra> >)>("domain discretization")
			.add_method("set_time", &T::set_time, "", "time")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "Jacobi_c", tag);
	}

	// constrained ILU
	{
		using T = ConstrainedLinearIterator<TDomain, TAlgebra, ILU<TAlgebra> >;
		using TBase = ILU<TAlgebra>;
		string name = string("ILU_c").append(suffix);
		reg.add_class_<T,TBase>(name, precondGrp, "ILU preconditioner respecting constraints")
		.template add_constructor<void (*)(SmartPtr<IDomainDiscretization<TAlgebra> >)>("domain discretization")
			.add_method("set_time", &T::set_time, "", "time")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ILU_c", tag);
	}

	// constrained ILU-T
	{
		using T = ConstrainedLinearIterator<TDomain, TAlgebra, ILUTPreconditioner<TAlgebra> >;
		using TBase = ILUTPreconditioner<TAlgebra>;
		string name = string("ILUT_c").append(suffix);
		reg.add_class_<T,TBase>(name, precondGrp, "ILUT preconditioner respecting constraints")
		.template add_constructor<void (*)(SmartPtr<IDomainDiscretization<TAlgebra> >)>("domain discretization")
			.add_method("set_time", &T::set_time, "", "time")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ILUT_c", tag);
	}

	// constrained GMG
	{
		using T = ConstrainedLinearIterator<TDomain, TAlgebra, AssembledMultiGridCycle<TDomain, TAlgebra> >;
		using TBase = AssembledMultiGridCycle<TDomain, TAlgebra>;
		string name = string("GMG_c").append(suffix);
		reg.add_class_<T,TBase>(name, precondGrp, "GMG preconditioner respecting constraints")
		.template add_constructor<void (*)(SmartPtr<IDomainDiscretization<TAlgebra> >)>("domain discretization")
			.add_method("set_time", &T::set_time, "", "time")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "GMG_c", tag);
	}

	// constrained LU
	{
		using T = ConstrainedLinearIterator<TDomain, TAlgebra, LU<TAlgebra> >;
		using TBase = LU<TAlgebra>;
		string name = string("LU_c").append(suffix);
		reg.add_class_<T,TBase>(name, solverGrp, "LU solver respecting constraints")
		.template add_constructor<void (*)(SmartPtr<IDomainDiscretization<TAlgebra> >)>("domain discretization")
			.add_method("set_time", &T::set_time, "", "time")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "LU_c", tag);
	}
	
};

}; // end CLI_Functionality

// end group constr_lin_it_bridge
/// \}

/// \addtogroup constr_lin_it_bridge
void RegisterBridge_ConstrainedLinearIterator(Registry& reg, string grp)
{
	grp.append("/Algebra");

	try {RegisterDomainAlgebraDependent<CLI_Functionality>(reg,grp);}
	UG_REGISTRY_CATCH_THROW(grp);
}

} // namespace bridge
} // namespace ug
