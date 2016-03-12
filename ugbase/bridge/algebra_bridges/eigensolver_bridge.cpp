/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
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
#include "bridge/util_algebra_dependent.h"

#include "lib_algebra/lib_algebra.h"

using namespace std;


#include "lib_algebra/operator/eigensolver/pinvit.h"
#include "lib_algebra/operator/preconditioner/operator_inverse_iterator.h"

namespace ug{
namespace bridge{
namespace Eigensolver{

/**
 * \defgroup eigensolver_bridge Eigensolver Bridge
 * \ingroup algebra_bridge
 * \{
 */

/**
 * Class exporting the functionality. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */
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

#ifdef LAPACK_AVAILABLE
	{
		string name = string("EigenSolver").append(suffix);
		typedef PINVIT<TAlgebra> T;
		typedef DebugWritingObject<TAlgebra> TBase;
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.add_method("add_vector", &T::add_vector, "", "vector")
			.add_method("set_preconditioner", &T::set_preconditioner, "", "Preconditioner")
			.add_method("set_linear_operator_A", &T::set_linear_operator_A,	"", "LinearOperatorA")
			.add_method("set_linear_operator_B", &T::set_linear_operator_B,	"", "LinearOperatorB")
			.add_method("set_max_iterations", &T::set_max_iterations, "", "precision")
			.add_method("num_eigenvalues", &T::num_eigenvalues, "number of eigenvalues")
			.add_method("get_eigenvalue", &T::get_eigenvalue, "eigenvalue i", "i")
			.add_method("set_precision", &T::set_precision, "", "precision")
			.add_method("set_relative_precision", &T::set_relative_precision, "", "precision")
			.add_method("get_eigenvector", &T::get_eigenvector, "eigenvector i", "i")
			.add_method("set_print_projected_eigenvectors", &T::set_print_projected_eigenvectors, "", "b")
			.add_method("set_print_projected_eigenvalues", &T::set_print_projected_eigenvalues, "", "b")
			.add_method("set_print_projected_eigenproblem", &T::set_print_projected_eigenproblem, "", "b")
			.add_method("set_additional_eigenvectors_to_keep", &T::set_additional_eigenvectors_to_keep, "", "i", "additional i smallest ritz-eigenvectors are added to the test space")
			.add_method("set_use_additional_corrections", &T::set_use_additional_corrections, "", "b")
			.add_method("set_debug_calc_projected_eigenvalues", &T::set_debug_calc_projected_eigenvalues, "", "b")
			.add_method("set_print_used_testvectors", &T::set_print_used_testvectors, "", "b")
			.add_method("set_pinvit", &T::set_pinvit, "", "iPINVIT", "1 = preconditioned inverse block iteration, 2 = preconditioned block gradient descent, 3 = LOBPCG")
			.add_method("set_linear_dependent_eps", &T::set_linear_dependent_eps)
			.add_method("set_laplacian", &T::set_laplacian)
			.add_method("get_max_deflection_of_a_mode", &T::get_max_deflection_of_a_mode, "maximal deflection of a mode", "composed mode and eigenvectors,")
			.add_method("set_store_defects", &T::set_store_defects)
			.add_method("set_store_lambdas", &T::set_store_lambdas)
			.add_method("get_defect", &T::get_defect)
			.add_method("get_lambda", &T::get_lambda)
			.add_method("get_iterations", &T::get_iterations)
			.add_method("apply", &T::apply);




		reg.add_class_to_group(name, "EigenSolver", tag);
	}
#endif
	{
		string name = string("OperatorInverseIterator").append(suffix);
		typedef typename TAlgebra::vector_type vector_type;
		typedef ILinearIterator<vector_type> TBase;
		typedef OperatorInverseIterator<TAlgebra> T;
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(SmartPtr<ILinearOperatorInverse<vector_type>  >)>( "opInv")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "OperatorInverseIterator", tag);
	}
}
}; // end Functionality

// end group eigensolver_bridge
/// \}

}// end Eigensolver


/// \addtogroup eigensolver_bridge
void RegisterBridge_Eigensolver(Registry& reg, string grp)
{
	grp.append("/Algebra/Solver");
	typedef Eigensolver::Functionality Functionality;

	try{
		RegisterAlgebraDependent<Functionality>(reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

} // namespace bridge
} // namespace ug
