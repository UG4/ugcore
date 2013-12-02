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
#include "bridge/util_algebra_dependent.h"

// solver
#include "lib_algebra/lib_algebra.h"
#include "lib_algebra/operator/damping.h"
#include "lib_algebra/operator/linear_solver/linear_solver.h"
#include "lib_algebra/operator/linear_solver/auto_linear_solver.h"
#include "lib_algebra/operator/linear_solver/analyzing_solver.h"
#include "lib_algebra/operator/linear_solver/cg.h"
#include "lib_algebra/operator/linear_solver/bicgstab.h"
#include "lib_algebra/operator/linear_solver/gmres.h"
#include "lib_algebra/operator/linear_solver/lu.h"
#include "lib_algebra/operator/linear_solver/agglomerating_solver.h"
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
 * \defgroup solver_bridge Solver Bridge
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
		reg.add_class_<T,TBase>(name, grp, "Minimal Residdum Damping (damping computed based on the minimal residuum)")
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "MinimalResiduumDamping", tag);
	}
	
// 	MinimalEngergyDamping
	{
		typedef MinimalEnergyDamping<vector_type> T;
		typedef IDamping<vector_type> TBase;
		string name = string("MinimalEnergyDamping").append(suffix);
		reg.add_class_<T,TBase>(name, grp, "Minimal Energy Damping (damping computed based on the minimal energy)")
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "MinimalEnergyDamping", tag);
	}

// 	LinearSolver
	{
		typedef LinearSolver<vector_type> T;
		typedef IPreconditionedLinearOperatorInverse<vector_type> TBase;
		string name = string("LinearSolver").append(suffix);
		reg.add_class_<T,TBase>(name, grp, "Linear Solver")
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "LinearSolver", tag);
	}


// 	AutoLinearSolver
	{
		typedef AutoLinearSolver<vector_type> T;
		typedef IPreconditionedLinearOperatorInverse<vector_type> TBase;
		string name = string("AutoLinearSolver").append(suffix);
		reg.add_class_<T,TBase>(name, grp, "Auto Linear Solver")
			.add_constructor()
			.template add_constructor<void (*)(double, double)>("reductionAlwaysAccept#worseThenAverage")
			.add_method("set_reduction_always_accept", &T::set_reduction_always_accept)
			.add_method("set_reinit_when_worse_then_average", &T::set_reinit_when_worse_then_average)
			.set_construct_as_smart_pointer(true);

		reg.add_class_to_group(name, "AutoLinearSolver", tag);
	}

// 	AnalyzingSolver
	{
		typedef AnalyzingSolver<matrix_type, vector_type> T;
		typedef  ILinearOperatorInverse<vector_type> TBase;
		string name = string("AnalyzingSolver").append(suffix);
		reg.add_class_<T,TBase>(name, grp, "AnalyzingSolver")
			.template add_constructor<void (*)(SmartPtr<ILinearOperatorInverse<vector_type> >)>("solver")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "AnalyzingSolver", tag);
	}

	// 	CG Solver
	{
		typedef CG<vector_type> T;
		typedef IPreconditionedLinearOperatorInverse<vector_type> TBase;
		string name = string("CG").append(suffix);
		reg.add_class_<T,TBase>(name, grp, "Conjugate Gradient Solver")
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "CG", tag);
	}

// 	BiCGStab Solver
	{
		typedef BiCGStab<vector_type> T;
		typedef IPreconditionedLinearOperatorInverse<vector_type> TBase;
		string name = string("BiCGStab").append(suffix);
		reg.add_class_<T,TBase>(name, grp, "BiCGStab Solver")
			.add_constructor()
			.add_method("set_restart", &T::set_restart)
			.add_method("set_min_orthogonality", &T::set_min_orthogonality)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "BiCGStab", tag);
	}

// 	GMRES Solver
	{
		typedef GMRES<vector_type> T;
		typedef IPreconditionedLinearOperatorInverse<vector_type> TBase;
		string name = string("GMRES").append(suffix);
		reg.add_class_<T,TBase>(name, grp, "GMRES Solver")
			.template add_constructor<void (*)(size_t)>("restart")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "GMRES", tag);
	}

// 	LU Solver
	{
		typedef LU<TAlgebra> T;
		typedef ILinearOperatorInverse<vector_type> TBase;
		string name = string("LU").append(suffix);
		reg.add_class_<T,TBase>(name, grp, "LU-Decomposition exact solver")
			.add_constructor()
			.add_method("set_minimum_for_sparse", &T::set_minimum_for_sparse, "", "N")
			.add_method("set_sort_sparse", &T::set_sort_sparse, "", "bSort", "if bSort=true, use a cuthill-mckey sorting to reduce fill-in in sparse LU. default true")
			.add_method("set_info", &T::set_info, "", "bInfo", "if true, sparse LU prints some fill-in info")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "LU", tag);
	}

// 	AgglomeratingSolver
	{
		typedef AgglomeratingSolver<TAlgebra> T;
		typedef ILinearOperatorInverse<vector_type> TBase;
		string name = string("AgglomeratingSolver").append(suffix);
		reg.add_class_<T,TBase>(name, grp, "AgglomeratingSolver")
			.template add_constructor<void (*)(SmartPtr<ILinearOperatorInverse<vector_type, vector_type> > )>("pLinOp")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "AgglomeratingSolver", tag);
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
		.add_method("set_neumann_solver", &T::set_neumann_solver, "", 
					"", "Neumann Solver")
		.add_method("set_dirichlet_solver", &T::set_dirichlet_solver, "",
					"", "Dirichlet Solver")
		.add_method("set_coarse_problem_solver", &T::set_coarse_problem_solver, "",
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

// end group solver_bridge
/// \}

}// end Solver

/// \addtogroup solver_bridge
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
