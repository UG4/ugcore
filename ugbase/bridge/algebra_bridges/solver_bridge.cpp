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
#include "bridge/util_algebra_dependent.h"
#include "bridge/util_domain_algebra_dependent.h" //todo: move to lib_disc

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
#include "lib_algebra/operator/linear_solver/distance_to_boundary_bruteforce.h"
#include "lib_algebra/operator/linear_solver/agglomerating_solver.h"
#include "lib_algebra/operator/linear_solver/debug_iterator.h"
#include "lib_algebra/operator/linear_solver/external_solvers/external_solvers.h"
#ifdef UG_PARALLEL
#include "lib_algebra/operator/linear_solver/feti.h"
#endif
#include "../util_overloaded.h"


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
	
//	Pre- and postprocess operations
	{
		typedef IPProcessVector<vector_type> T;
		string name = string("IPProcessVector").append(suffix);
		reg.add_class_<T>(name, grp)
			.add_method("apply", &T::apply, "applies the operation", "vector");
		reg.add_class_to_group(name, "IPProcessVector", tag);
	}

// 	LinearSolver
	{
		typedef LinearSolver<vector_type> T;
		typedef IPreconditionedLinearOperatorInverse<vector_type> TBase;
		string name = string("LinearSolver").append(suffix);
		reg.add_class_<T,TBase>(name, grp, "Linear Solver")
			.add_constructor()
			. ADD_CONSTRUCTOR( (SmartPtr<ILinearIterator<vector_type,vector_type> > ) )("precond")
			. ADD_CONSTRUCTOR( (SmartPtr<ILinearIterator<vector_type,vector_type> >, SmartPtr<IConvergenceCheck<vector_type> >) )("precond#convCheck")
			.add_method("add_postprocess_corr", &T::add_postprocess_corr, "adds a postprocess of the corrections", "op")
			.add_method("remove_postprocess_corr", &T::remove_postprocess_corr, "removes a postprocess of the corrections", "op")
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
			.ADD_CONSTRUCTOR( (double, double) )("reductionAlwaysAccept#worseThenAverage")
			.add_method("set_reduction_always_accept", &T::set_reduction_always_accept)
			.add_method("set_reinit_when_worse_then_average", &T::set_reinit_when_worse_then_average)
			.add_method("print_information", &T::print_information)
			.set_construct_as_smart_pointer(true);

		reg.add_class_to_group(name, "AutoLinearSolver", tag);
	}

// 	AnalyzingSolver
	{
		typedef AnalyzingSolver<matrix_type, vector_type> T;
		typedef  ILinearOperatorInverse<vector_type> TBase;
		string name = string("AnalyzingSolver").append(suffix);
		reg.add_class_<T,TBase>(name, grp, "AnalyzingSolver")
			.ADD_CONSTRUCTOR( (SmartPtr<ILinearOperatorInverse<vector_type> >) )("solver")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "AnalyzingSolver", tag);
	}

	//	Debug iteration
	{
		typedef DebugIterator<TAlgebra> T;
		typedef ILinearIterator<typename TAlgebra::vector_type> TBase;

		typedef ILinearIterator<typename TAlgebra::vector_type> base_type;
		typedef IPreconditionedLinearOperatorInverse<vector_type> solver_type;

		string name = string("DebugIterator").append(suffix);
		reg.add_class_<T,TBase>(name, grp, "DebugIterator")
		.add_constructor()
		. ADD_CONSTRUCTOR( (SmartPtr<base_type> pprecond, SmartPtr<solver_type> psolver) ) ("precond#solver")
		.add_method("set_preconditioner", &T::set_preconditioner)
		.add_method("set_solver", &T::set_solver)
		.add_method("set_debug", &T::set_debug)
		.add_method("set_solution", &T::set_solution)
		.add_method("set_random_bounds", &T::set_random_bounds)
		.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "DebugIterator", tag);
	}


	// 	CG Solver
	{
		typedef CG<vector_type> T;
		typedef IPreconditionedLinearOperatorInverse<vector_type> TBase;
		string name = string("CG").append(suffix);
		reg.add_class_<T,TBase>(name, grp, "Conjugate Gradient Solver")
			.add_constructor()
			. ADD_CONSTRUCTOR( (SmartPtr<ILinearIterator<vector_type,vector_type> > ) )("precond")
			. ADD_CONSTRUCTOR( (SmartPtr<ILinearIterator<vector_type,vector_type> >, SmartPtr<IConvergenceCheck<vector_type> >) )("precond#convCheck")
			.add_method("add_postprocess_corr", &T::add_postprocess_corr, "adds a postprocess of the corrections", "op")
			.add_method("remove_postprocess_corr", &T::remove_postprocess_corr, "removes a postprocess of the corrections", "op")
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
			. ADD_CONSTRUCTOR( (SmartPtr<ILinearIterator<vector_type,vector_type> > ) )("precond")
			. ADD_CONSTRUCTOR( (SmartPtr<ILinearIterator<vector_type,vector_type> >, SmartPtr<IConvergenceCheck<vector_type> >) )("precond#convCheck")
			.add_method("set_restart", &T::set_restart)
			.add_method("set_min_orthogonality", &T::set_min_orthogonality)
			.add_method("add_postprocess_corr", &T::add_postprocess_corr, "adds a postprocess of the corrections", "op")
			.add_method("remove_postprocess_corr", &T::remove_postprocess_corr, "removes a postprocess of the corrections", "op")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "BiCGStab", tag);
	}

// 	GMRES Solver
	{
		typedef GMRES<vector_type> T;
		typedef IPreconditionedLinearOperatorInverse<vector_type> TBase;
		string name = string("GMRES").append(suffix);
		reg.add_class_<T,TBase>(name, grp, "GMRES Solver")
			.ADD_CONSTRUCTOR( (size_t restar) )("restart")
			.add_method("add_postprocess_corr", &T::add_postprocess_corr, "adds a postprocess of the corrections", "op")
			.add_method("remove_postprocess_corr", &T::remove_postprocess_corr, "removes a postprocess of the corrections", "op")
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
			.add_method("set_show_progress", &T::set_show_progress, "", "onoff", "switches the progress indicator on/off")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "LU", tag);
	}

// 	AgglomeratingSolver
	{
		typedef AgglomeratingSolver<TAlgebra> T;
		typedef ILinearOperatorInverse<vector_type> TBase;
		string name = string("AgglomeratingSolver").append(suffix);
		reg.add_class_<T,TBase>(name, grp, "AgglomeratingSolver")
			.ADD_CONSTRUCTOR( (SmartPtr<ILinearOperatorInverse<vector_type, vector_type> > ) )("pLinOp")
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

	// ExternalSolver
	{
		typedef IExternalSolver<TAlgebra> T;
		string name = string("ExternalSolver").append(suffix);
		reg.add_class_<T>(name, grp)
			.add_method("set_disable_preprocessing", &T::set_disable_preprocessing, "", "", "")
			.add_method("enable_consistent_interfaces", &T::enable_consistent_interfaces, "", "", "");

	}


}


/* TODO: move to lib_disc */

/**
 * Function called for the registration of Domain and Algebra dependent parts.
 * All Functions and Classes depending on both Domain and Algebra
 * are to be placed here when registering. The method is called for all
 * available Domain and Algebra types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TDomain, typename TAlgebra>
static void DomainAlgebra(Registry& reg, string grp)
{
	string suffix = GetDomainAlgebraSuffix<TDomain,TAlgebra>();
	string tag = GetDomainAlgebraTag<TDomain,TAlgebra>();

//	typedefs for this algebra
	typedef typename TAlgebra::vector_type vector_type;

// 	DistanceToBoundaryBruteforce
	{
		typedef DistanceToBoundaryBruteforce<TAlgebra, TDomain> T;
		typedef ILinearOperatorInverse<vector_type> TBase;
		string name = string("DistanceToBoundaryBruteforce").append(suffix);
		reg.add_class_<T,TBase>(name, grp, "DistanceToBoundaryBruteforce")
			.add_constructor()
			.add_method("set_approximation_space", &T::set_approximation_space, "", "Approximation space")
			.add_method("select_inner", static_cast<void (T::*)(int)>(&T::select_inner))
			.add_method("select_boundary", static_cast<void (T::*)(int)>(&T::select_boundary))
			.add_method("set_level", static_cast<void (T::*)(int)>(&T::set_level))
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "DistanceToBoundaryBruteforce", tag);
	}
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
		RegisterDomainAlgebraDependent<Functionality>(reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

} // namespace bridge
} // namespace ug
