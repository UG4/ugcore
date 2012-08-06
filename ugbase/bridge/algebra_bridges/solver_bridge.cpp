/*
 * solver_bridge.cpp
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

// solver
#include "lib_algebra/lib_algebra.h"
#include "lib_algebra/operator/damping.h"
#include "lib_algebra/operator/linear_solver/linear_solver.h"
#include "lib_algebra/operator/linear_solver/cg.h"
#include "lib_algebra/operator/linear_solver/bicgstab.h"
#include "lib_algebra/operator/linear_solver/lu.h"
#ifdef UG_PARALLEL
#include "lib_algebra/operator/linear_solver/feti.h"
	#ifdef UG_HLIBPRO
	#include "lib_algebra/operator/linear_solver/hlibpro.h"
	#endif
#endif

using namespace std;

namespace ug{
namespace bridge{
namespace Solver{

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


// 	IDamping
	{
		typedef IDamping<vector_type> T;
		string name = string("IDamping").append(suffix);
		reg.add_class_<T>(name, grp);
		reg.add_class_to_group(name, "IDamping", tag);
	}

// 	MinimalResiduumDamping
	{
		typedef MinimalResiduumDamping<vector_type> T;
		typedef IDamping<vector_type> TBase;
		string name = string("MinimalResiduumDamping").append(suffix);
		reg.add_class_<T,TBase>(name, grp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "MinimalResiduumDamping", tag);
	}
	
// 	MinimalEngergyDamping
	{
		typedef MinimalEnergyDamping<vector_type> T;
		typedef IDamping<vector_type> TBase;
		string name = string("MinimalEnergyDamping").append(suffix);
		reg.add_class_<T,TBase>(name, grp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "MinimalEnergyDamping", tag);
	}

// 	LinearSolver
	{
		typedef LinearSolver<vector_type> T;
		typedef IPreconditionedLinearOperatorInverse<vector_type> TBase;
		string name = string("LinearSolver").append(suffix);
		reg.add_class_<T,TBase>(name, grp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "LinearSolver", tag);
	}

// 	CG Solver
	{
		typedef CG<vector_type> T;
		typedef IPreconditionedLinearOperatorInverse<vector_type> TBase;
		string name = string("CG").append(suffix);
		reg.add_class_<T,TBase>(name, grp, "Conjugate Gradient")
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "CG", tag);
	}

// 	BiCGStab Solver
	{
		typedef BiCGStab<vector_type> T;
		typedef IPreconditionedLinearOperatorInverse<vector_type> TBase;
		string name = string("BiCGStab").append(suffix);
		reg.add_class_<T,TBase>(name, grp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "BiCGStab", tag);
	}

// 	LU Solver
	{
		typedef LU<TAlgebra> T;
		typedef ILinearOperatorInverse<vector_type> TBase;
		string name = string("LU").append(suffix);
		reg.add_class_<T,TBase>(name, grp, "LU-Decomposition exact solver")
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "LU", tag);
	}

#ifdef UG_PARALLEL
// 	LocalSchurComplement
	{
		typedef LocalSchurComplement<TAlgebra> T;
		typedef ILinearOperator<vector_type> TBase;
		typedef DebugWritingObject<TAlgebra> TBase2;
		string name = string("LocalSchurComplement").append(suffix);
		reg.add_class_<	T, TBase, TBase2>(name, grp)
		.add_constructor()
		.add_method("set_matrix", &T::set_matrix,
					"", "Matrix")
		.add_method("set_dirichlet_solver", &T::set_dirichlet_solver,
					"", "Dirichlet Solver")
		// the following functions would normally not be executed from script
		.add_method("init", static_cast<void (T::*)()>(&T::init))
		.add_method("apply", &T::apply,
					"Success", "local SC times Vector#Vector")
		.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "LocalSchurComplement", tag);
	}

// 	FETISolver
	{
		typedef FETISolver<TAlgebra> T;
		typedef IMatrixOperatorInverse<matrix_type, vector_type> BaseT;
		typedef DebugWritingObject<TAlgebra> TBase2;
		string name = string("FETI").append(suffix);
		reg.add_class_<	T, BaseT,TBase2>(name, grp, "FETI Domain Decomposition Solver")
		.add_constructor()
		.add_method("set_neumann_solver", &T::set_neumann_solver,
					"", "Neumann Solver")
		.add_method("set_dirichlet_solver", &T::set_dirichlet_solver,
					"", "Dirichlet Solver")
		.add_method("set_coarse_problem_solver", &T::set_coarse_problem_solver,
					"", "Coarse Problem Solver")
		.add_method("set_domain_decomp_info", &T::set_domain_decomp_info)
		.add_method("print_statistic_of_inner_solver", &T::print_statistic_of_inner_solver)
		.add_method("set_debug", &T::set_debug)
		.add_method("test_layouts", &T::test_layouts)
		.add_method("set_test_one_to_many_layouts", &T::set_test_one_to_many_layouts)
		.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "FETI", tag);
	}
#endif

	// 	HLIBSolver
#ifdef UG_HLIBPRO
	{
		typedef HLIBSolver<TAlgebra> T;
		typedef ILinearOperatorInverse<vector_type> TBase;
		typedef DebugWritingObject<TAlgebra> TBase2;
		string name = string("HLIBSolver").append(suffix);
		reg.add_class_<	T, TBase, TBase2>(name, grp)
		.add_constructor()
		.add_method("set_hlib_nmin",         &T::set_hlib_nmin,
					"", "HLIB nmin")
		.add_method("set_hlib_accuracy_H",   &T::set_hlib_accuracy_H,
					"", "HLIB accuracy_H")
		.add_method("set_hlib_accuracy_LU",  &T::set_hlib_accuracy_LU,
					"", "HLIB accuracy_LU")
		.add_method("set_hlib_verbosity",    &T::set_hlib_verbosity,
					"", "HLIB verbosity")
		.add_method("set_clustering_method", &T::set_clustering_method,
					"", "Clustering")
		.add_method("set_ps_basename",       &T::set_ps_basename,
					"", "PostScript basename")
		.add_method("check_crs_matrix",      &T::check_crs_matrix,
					"", "Check CRS matrix")
		.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "HLIBSolver", tag);
	}
#endif

}
}; // end Functionality
}// end Solver

void RegisterBridge_Solver(Registry& reg, string grp)
{
	grp.append("/Algebra/Solver");
	typedef Solver::Functionality Functionality;

	try{
		RegisterAlgebraDependent<Functionality>(reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

} // namespace bridge
} // namespace ug
