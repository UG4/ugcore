/*
 * external_solvers_bridge.cpp
 *
 *  Created on: 08.01.2014
 *      Author: mrupp
 */

// extern headers
#include <iostream>
#include <sstream>
#include <string>

// include bridge
#include "bridge/bridge.h"
#include "bridge/util.h"
#include "bridge/util_algebra_dependent.h"

// preconditioner
#include "lib_algebra/lib_algebra.h"
#ifdef UG_PARALLEL
#include "lib_algebra/operator/linear_solver/external_solvers/external_solvers.h"
#include "lib_algebra/operator/linear_solver/external_solvers/super_lu.h"
#endif

#include "../util_overloaded.h"
using namespace std;

namespace ug{

namespace bridge{
namespace ExternalSolvers{


struct Functionality
{


/**
 * Function called for the registration of Algebra dependent parts.
 * All Functions and Classes depending on Algebra
 * are to be placed here when registering. The method is called for all
 * available Algebra types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TAlgebra>
static void Algebra(Registry& reg, string grp)
{
	string suffix = GetAlgebraSuffix<TAlgebra>();
	string tag = GetAlgebraTag<TAlgebra>();

//	typedefs for this algebra
	typedef typename TAlgebra::vector_type vector_type;
	typedef typename TAlgebra::matrix_type matrix_type;

#ifdef UG_SUPERLU
	// 	SuperLU Solver
	{
		typedef SuperLUSolver<TAlgebra> T;
		typedef ILinearOperatorInverse<vector_type> TBase;
		string name = string("SuperLU").append(suffix);
		reg.add_class_<T,TBase>(name, grp, "SuperLU")
			.add_constructor()
			.add_method("print_stat", &T::print_stat)
			.add_method("equil", &T::equil)
			.add_method("col_perm_natural", &T::col_perm_natural)
			.add_method("col_perm_mdo_ATA", &T::col_perm_mdo_ATA)
			.add_method("col_perm_mdo_AT_plus_A", &T::col_perm_mdo_AT_plus_A)
			.add_method("col_perm_approximate", &T::col_perm_approximate)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "SuperLU", tag);
	}
#else

	// reg.add_notexisting_class("SuperLU", "you need to uginstall SuperLU")
#endif
}


}; // end Functionality

// end group precond_bridge
/// \}

}// end Preconditioner

/// \addtogroup precond_bridge
void RegisterBridge_ExternalSolvers(Registry& reg, string grp)
{
	grp.append("/Algebra/ExternalSolvers/");
	typedef ExternalSolvers::Functionality Functionality;

	try{
		RegisterAlgebraDependent<Functionality>(reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

} // namespace bridge
} // namespace ug
