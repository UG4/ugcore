/*
 * eigensolver_bridge.cpp
 *
 *  Created on: 07.07.2011
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

#include "lib_algebra/lib_algebra.h"

using namespace std;


#include "lib_algebra/operator/eigensolver/pinvit.h"
#include "lib_algebra/operator/preconditioner/operator_inverse_iterator.h"

namespace ug{
namespace bridge{
namespace Eigensolver{

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

//	typedefs for this algebra
	typedef typename TAlgebra::vector_type vector_type;
	typedef typename TAlgebra::matrix_type matrix_type;

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
			.add_method("get_eigenvector", &T::get_eigenvector, "eigenvector i", "i")
			.add_method("set_print_projected_eigenvectors", &T::set_print_projected_eigenvectors, "", "b")
			.add_method("set_print_projected_eigenvalues", &T::set_print_projected_eigenvalues, "", "b")
			.add_method("set_print_projected_eigenproblem", &T::set_print_projected_eigenproblem, "", "b")
			.add_method("set_additional_eigenvectors_to_keep", &T::set_additional_eigenvectors_to_keep, "", "i", "additional i smallest ritz-eigenvectors are added to the test space")
			.add_method("set_use_additional_corrections", &T::set_use_additional_corrections, "", "b")
			.add_method("set_debug_calc_projected_eigenvalues", &T::set_debug_calc_projected_eigenvalues, "", "b")
			.add_method("set_print_used_testvectors", &T::set_print_used_testvectors, "", "b")
			.add_method("set_pinvit", &T::set_pinvit, "", "iPINVIT", "1 = preconditioned inverse block iteration, 2 = preconditioned block gradient descent, 3 = LOBPCG")
			.add_method("apply", &T::apply);

		reg.add_class_to_group(name, "EigenSolver", tag);
	}
#endif
	{
		string name = string("OperatorInverseIterator").append(suffix);
		typedef ILinearIterator<typename TAlgebra::vector_type> TBase;
		typedef OperatorInverseIterator<TAlgebra> T;
		typedef typename TAlgebra::vector_type vector_type;
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(SmartPtr<ILinearOperatorInverse<vector_type>  >)>( "opInv")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "OperatorInverseIterator", tag);
	}
}
}; // end Functionality
}// end Eigensolver


void RegisterBridge_Eigensolver(Registry& reg, string grp)
{
	grp.append("/Algebra/Solver");
	typedef Eigensolver::Functionality Functionality;
//#ifdef UG_CPU_1
	typedef boost::mpl::list<CPUAlgebra> AlgList;
/*#else
	typedef boost::mpl::list<> AlgList;
#endif*/

	try{
		RegisterAlgebraDependent<Functionality>(reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

} // namespace bridge
} // namespace ug
