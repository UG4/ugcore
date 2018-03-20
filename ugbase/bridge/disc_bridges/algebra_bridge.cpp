/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
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

// extern libraries
#include <iostream>
#include <sstream>
#include <string>

#include "../../lib_disc/time_disc/composite_time_disc.h"
// include bridge
#include "bridge/bridge.h"
#include "bridge/util.h"
#include "bridge/util_domain_algebra_dependent.h"
#include "bridge/util_overloaded.h"

// discretization interfaces
#include "lib_algebra/active_set/active_set.h"
#include "lib_algebra/active_set/lagrange_multiplier_disc_interface.h"
#include "lib_algebra/operator/convergence_check.h"
#include "lib_algebra/operator/matrix_operator_functions.h"
#include "lib_disc/spatial_disc/domain_disc_interface.h"
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"
#include "lib_disc/spatial_disc/constraints/constraint_interface.h"
#include "lib_disc/time_disc/time_disc_interface.h"
#include "lib_disc/time_disc/theta_time_step.h"
#include "lib_disc/operator/linear_operator/assembled_linear_operator.h"
#include "lib_disc/operator/non_linear_operator/assembled_non_linear_operator.h"
#include "lib_disc/operator/non_linear_operator/line_search.h"
#include "lib_disc/operator/linear_operator/nested_iteration/nested_iteration.h"
#include "lib_disc/operator/non_linear_operator/newton_solver/newton.h"
#include "lib_disc/operator/non_linear_operator/nl_gauss_seidel/nl_gauss_seidel.h"
#include "lib_disc/operator/non_linear_operator/nl_jacobi/nl_jacobi.h"
#include "lib_disc/operator/composite_conv_check.h"
#include "lib_disc/spatial_disc/local_to_global/local_to_global_mapper.h"

using namespace std;


namespace ug{
namespace bridge{
namespace DiscAlgebra{

/**
 * \defgroup discalgebra_bridge Discretization Algebra Bridge
 * \ingroup disc_bridge
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
static void Algebra(Registry& reg, string parentGroup)
{
//	typedefs for Vector and Matrix
	typedef typename TAlgebra::vector_type vector_type;
	typedef typename TAlgebra::matrix_type matrix_type;

//	suffix and tag
	string suffix = GetAlgebraSuffix<TAlgebra>();
	string tag = GetAlgebraTag<TAlgebra>();

//	IConstraint
	{
		std::string grp = parentGroup; grp.append("/Discretization/SpatialDisc");
		typedef IConstraint<TAlgebra> T;
		string name = string("IConstraint").append(suffix);
		reg.add_class_<T>(name, grp);
		reg.add_class_to_group(name, "IConstraint", tag);
	}

//	IAssemble
	{
		std::string grp = parentGroup; grp.append("/Discretization/SpatialDisc");
		typedef IAssemble<TAlgebra> T;
		string name = string("IAssemble").append(suffix);
		reg.add_class_<T>(name, grp)
			.add_method("assemble_jacobian", static_cast<void (T::*)(matrix_type&, const vector_type&)>(&T::assemble_jacobian),"", "J(u)#u", "assembles jacobian on surface grid")
			.add_method("assemble_jacobian", static_cast<void (T::*)(matrix_type&, const vector_type&, const GridLevel&)>(&T::assemble_jacobian),"", "J(u)#u#GridLevel", "assembles jacobian on grid level")
			.add_method("assemble_defect", static_cast<void (T::*)(vector_type&, const vector_type&)>(&T::assemble_defect),"", "d(u)#u", "Assembles Defect on surface grid.")
			.add_method("assemble_defect", static_cast<void (T::*)(vector_type&, const vector_type&, const GridLevel&)>(&T::assemble_defect),"", "d(u)#u#GridLevel", "Assembles Defect on grid level")
			.add_method("assemble_linear", static_cast<void (T::*)(matrix_type&, vector_type&)>(&T::assemble_linear),"", "A#b", "Assembles Matrix and rhs on surface grid.")
			.add_method("assemble_linear", static_cast<void (T::*)(matrix_type&, vector_type&, const GridLevel&)>(&T::assemble_linear),"", "A#b#GridLevel", "Assembles Matrix and rhs on grid level.")
			.add_method("assemble_rhs", static_cast<void (T::*)(vector_type&, const vector_type&)>(&T::assemble_rhs),"", "rhs#u", "assembles right-hand side on surface grid")
			.add_method("assemble_rhs", static_cast<void (T::*)(vector_type&, const vector_type&, const GridLevel&)>(&T::assemble_rhs),"", "rhs#u#GridLevel", "assembles right-hand side on grid level")
			.add_method("assemble_rhs", static_cast<void (T::*)(vector_type&)>(&T::assemble_rhs),"", "rhs#GridLevel", "assembles right-hand side on surface grid for linear case")
			.add_method("assemble_rhs", static_cast<void (T::*)(vector_type&, const GridLevel&)>(&T::assemble_rhs),"", "rhs", "assembles right-hand side on grid level for linear case")
			.add_method("assemble_stiffness_matrix", static_cast<void (T::*)(matrix_type&, const vector_type&)>(&T::assemble_stiffness_matrix),"", "A#u", "assembles stiffness matrix on surface grid")
			.add_method("assemble_mass_matrix", static_cast<void (T::*)(matrix_type&, const vector_type&)>(&T::assemble_mass_matrix),"", "M#u", "assembles mass matrix on surface grid")
			.add_method("adjust_solution", static_cast<void (T::*)(vector_type&)>(&T::adjust_solution))
			.add_method("adjust_solution", static_cast<void (T::*)(vector_type&, const GridLevel&)>(&T::adjust_solution));
		reg.add_class_to_group(name, "IAssemble", tag);
	}

//	IDomainDiscretization
	{
		std::string grp = parentGroup; grp.append("/Discretization/SpatialDisc");
		typedef IAssemble<TAlgebra> TBase;
		typedef IDomainDiscretization<TAlgebra> T;
		string name = string("IDomainDiscretization").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_method("assemble_jacobian", static_cast<void (T::*)
			            (matrix_type&, ConstSmartPtr<VectorTimeSeries<vector_type> >,
			             number)>(&T::assemble_jacobian));
		reg.add_class_to_group(name, "IDomainDiscretization", tag);
	}

//	ITimeDiscretization
	{
		std::string grp = parentGroup; grp.append("/Discretization/TimeDisc");
		typedef IAssemble<TAlgebra>  TBase;
		typedef ITimeDiscretization<TAlgebra> T;
		string name = string("ITimeDiscretization").append(suffix);
		reg.add_class_<T,TBase>(name, grp)
			.add_method("prepare_step", &T::prepare_step, "", "", "prepares the assembling of defect/Jacobian for a time step")
			.add_method("prepare_step_elem", static_cast<void (T::*)(SmartPtr<VectorTimeSeries<vector_type> >, number)>(&T::prepare_step_elem))
			.add_method("finish_step", &T::finish_step, "", "", "finishes the assembling of defect/Jacobian for a time step")
			.add_method("finish_step_elem", static_cast<void (T::*)(SmartPtr<VectorTimeSeries<vector_type> >)>(&T::finish_step_elem))
			.add_method("finish_step_elem", static_cast<void (T::*)(SmartPtr<VectorTimeSeries<vector_type> >, const GridLevel&)>(&T::finish_step_elem))
			.add_method("num_stages", &T::num_stages, "the number of stages")
			.add_method("set_stage", &T::set_stage, "", "stage")
			.add_method("future_time", &T::future_time, "the future time point (i.e. the one that will be computed)")
			.add_method("num_prev_steps", &T::num_prev_steps);
		reg.add_class_to_group(name, "ITimeDiscretization", tag);
	}

//	MultiStepTimeDiscretization
	{
		std::string grp = parentGroup; grp.append("/Discretization/TimeDisc");
		typedef ITimeDiscretization<TAlgebra> TBase;
		typedef MultiStepTimeDiscretization<TAlgebra> T;
		string name = string("MultiStepTimeDiscretization").append(suffix);
		reg.add_class_<T,TBase>(name, grp)
			.add_method("calc_error", static_cast<void (T::*)(const vector_type&)>(&T::calc_error), "", "",
				"calculate error indicators for elements from error estimators of the elemDiscs")
			.add_method("calc_error", static_cast<void (T::*)(const vector_type&, vector_type&)>(&T::calc_error), "", "",
				"calculate error indicators for elements from error estimators of the elemDiscs")
			.add_method("invalidate_error", &T::invalidate_error, "", "Marks error indicators as invalid, "
				"which will prohibit refining and coarsening before a new call to calc_error.")
			.add_method("is_error_valid", &T::is_error_valid, "", "Returns whether error indicators are valid");
		reg.add_class_to_group(name, "MultiStepTimeDiscretization", tag);
	}

//	ThetaTimeStep
	{
		std::string grp = parentGroup; grp.append("/Discretization/TimeDisc");
		typedef MultiStepTimeDiscretization<TAlgebra> TBase;
		typedef ThetaTimeStep<TAlgebra> T;
		string name = string("ThetaTimeStep").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
				.template add_constructor<void (*)(SmartPtr<IDomainDiscretization<TAlgebra> >)>("Domain Discretization")
				.template add_constructor<void (*)(SmartPtr<IDomainDiscretization<TAlgebra> >,number)>("Domain Discretization#Theta")
				.template add_constructor<void (*)(SmartPtr<IDomainDiscretization<TAlgebra> >,const char*)>("Domain Discretization#Scheme")
				.add_method("set_theta", &T::set_theta, "", "Theta", "Theta value 1 = Impl; 0 = Expl")
				.add_method("set_scheme", &T::set_scheme, "", "Scheme|selection|value=[\"Theta\",\"Alexander\",\"FracStep\"]")
				.add_method("set_stage", &T::set_stage, "", "Stage")
				.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ThetaTimeStep", tag);
	}

//	BDF
	{
		std::string grp = parentGroup; grp.append("/Discretization/TimeDisc");
		typedef MultiStepTimeDiscretization<TAlgebra> TBase;
		typedef BDF<TAlgebra> T;
		string name = string("BDF").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
				.template add_constructor<void (*)(SmartPtr<IDomainDiscretization<TAlgebra> >)>("Domain Discretization")
				.template add_constructor<void (*)(SmartPtr<IDomainDiscretization<TAlgebra> >,int)>("Domain Discretization#Order")
				.add_method("set_order", &T::set_order, "", "Order")
				.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "BDF", tag);
	}

//	SDIRK
	{
		std::string grp = parentGroup; grp.append("/Discretization/TimeDisc");
		typedef MultiStepTimeDiscretization<TAlgebra> TBase;
		typedef SDIRK<TAlgebra> T;
		string name = string("SDIRK").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
				.template add_constructor<void (*)(SmartPtr<IDomainDiscretization<TAlgebra> >)>("Domain Discretization")
				.template add_constructor<void (*)(SmartPtr<IDomainDiscretization<TAlgebra> >,int)>("Domain Discretization#Order")
				.add_method("set_order", &T::set_order, "", "Order")
				.add_method("set_stage", &T::set_stage, "", "Stage")
				.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "SDIRK", tag);
	}

//	CompositeTimeDiscretization
	{
		std::string grp = parentGroup; grp.append("/Discretization/TimeDisc");
		typedef IAssemble<TAlgebra> TBase;
		typedef CompositeTimeDiscretization<TAlgebra> T;
		string name = string("CompositeTimeDiscretization").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.add_method("add_time_disc", &T::add_time_disc, "", "time discretization")
			.add_method("prepare_step", &T::prepare_step, "", "", "prepares the assembling of defect/Jacobian for a time step")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "CompositeTimeDiscretization", tag);
	}

//	AssembledLinearOperator
	{
		std::string grp = parentGroup; grp.append("/Discretization");
		typedef AssembledLinearOperator<TAlgebra> T;
		typedef MatrixOperator<matrix_type, vector_type> TBase;
		string name = string("AssembledLinearOperator").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.template add_constructor<void (*)(SmartPtr<IAssemble<TAlgebra> >)>("Assembling Routine")
			.template add_constructor<void (*)(SmartPtr<IAssemble<TAlgebra> >, const GridLevel&)>("AssemblingRoutine#GridLevel")
			.add_method("set_discretization", &T::set_discretization)
			.add_method("set_level", &T::set_level)
			.add_method("set_dirichlet_values", &T::set_dirichlet_values)
			.add_method("init_op_and_rhs", &T::init_op_and_rhs)
			.add_method("level", &T::level)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "AssembledLinearOperator", tag);
	}
	

//	NewtonSolver
	{
		std::string grp = parentGroup; grp.append("/Discretization/Nonlinear");
		typedef NewtonSolver<TAlgebra> T;
		typedef IOperatorInverse<vector_type> TBase;
		typedef DebugWritingObject<TAlgebra> TBase2;
		string name = string("NewtonSolver").append(suffix);
		reg.add_class_<T, TBase, TBase2>(name, grp)
			.add_constructor()
			.template add_constructor<void (*)(SmartPtr<IOperator<vector_type> >)>("Operator")
			.template add_constructor<void (*)(SmartPtr<IAssemble<TAlgebra> >)>("AssemblingRoutine")
			.add_method("set_linear_solver", &T::set_linear_solver, "", "linSolver")
			.add_method("set_convergence_check", &T::set_convergence_check, "", "convCheck")
			.add_method("set_line_search", &T::set_line_search, "", "lineSeach")
			.add_method("disable_line_search", &T::disable_line_search)
			.add_method("line_search", &T::line_search, "lineSeach", "")
			.add_method("init", &T::init, "success", "op")
			.add_method("prepare", &T::prepare, "success", "u")
			.add_method("apply", &T::apply, "success", "u")
			.add_method("print_average_convergence", &T::print_average_convergence)
			.add_method("clear_average_convergence", &T::clear_average_convergence)
			.add_method("num_newton_steps", &T::num_newton_steps, "number of newton steps in history")
			.add_method("num_linsolver_calls", &T::num_linsolver_calls, "number of linsolver calls in iNewtonStep", "iNewtonStep")
			.add_method("num_linsolver_steps", &T::num_linsolver_steps, "number of linsolver steps in newton step iNewtonStep", "iNewtonStep")
			.add_method("average_linear_steps", &T::average_linear_steps, "average number of linsolver steps per linsolver call in newton step iNewtonStep", "iNewtonStep")
			.add_method("total_linsolver_calls", &T::total_linsolver_calls, "total number of linsolver calls", "")
			.add_method("total_linsolver_steps", &T::total_linsolver_steps, "total number of linsolver steps", "")
			.add_method("total_average_linear_steps", &T::total_average_linear_steps, "total average number of linsolver steps per linsolver call", "")
			.add_method("last_num_newton_steps", &T::last_num_newton_steps, "Number of newton steps performed in last iteration")
			.add_method("add_inner_step_update", &T::add_inner_step_update, "data update called before every linsolver step", "")
			.add_method("clear_inner_step_update", &T::clear_inner_step_update, "clear inner step update", "")
			.add_method("add_step_update", &T::add_step_update, "data update called before every Newton step", "")
			.add_method("clear_step_update", &T::clear_step_update, "clear step update", "")
			.add_method("config_string", &T::config_string)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NewtonSolver", tag);
	}

	//	NonlinearJacobiSolver
	{
		std::string grp = parentGroup; grp.append("/Discretization/Nonlinear");
		typedef NLJacobiSolver<TAlgebra> T;
		typedef IOperatorInverse<vector_type> TBase;
		typedef DebugWritingObject<TAlgebra> TBase2;
		string name = string("NLJacobiSolver").append(suffix);
		reg.add_class_<T, TBase, TBase2>(name, grp)
			.add_constructor()
			.template add_constructor<void (*)(SmartPtr<IConvergenceCheck<vector_type> >)>("ConvCheck")
			.add_method("set_convergence_check", &T::set_convergence_check, "", "convCheck")
			.add_method("set_damp", &T::set_damp, "", "setDampingFactor")
			.add_method("init", &T::init, "success", "op")
			.add_method("prepare", &T::prepare, "success", "u")
			.add_method("apply", &T::apply, "success", "u")
			.add_method("config_string", &T::config_string)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NLJacobiSolver", tag);
	}

//	AssembledOperator
	{
		std::string grp = parentGroup; grp.append("/Discretization/Nonlinear");
		typedef AssembledOperator<TAlgebra> T;
		typedef IOperator<vector_type> TBase;
		string name = string("AssembledOperator").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.template add_constructor<void (*)(SmartPtr<IAssemble<TAlgebra> >)>("AssemblingRoutine")
			.template add_constructor<void (*)(SmartPtr<IAssemble<TAlgebra> >, const GridLevel&)>("AssemblingRoutine#GridLevel")
			.add_method("set_discretization", &T::set_discretization, "", "ass")
			.add_method("set_level", &T::set_level, "", "gridLevel")
			.add_method("level", &T::level, "gridLevel", "")
			.add_method("init", &T::init)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "AssembledOperator", tag);
	}


//	ILocalToGlobalMapper
	{
		std::string grp = parentGroup; grp.append("/Discretization/SpatialDisc");
		typedef ILocalToGlobalMapper<TAlgebra> T;
		string name = string("ILocalToGlobalMapper").append(suffix);
		reg.add_class_<T>(name, grp);
		reg.add_class_to_group(name, "ILocalToGlobalMapper", tag);
	}

//	some functions
	{
		std::string grp = parentGroup; grp.append("/Discretization");
		reg.add_function("AssembleLinearOperatorRhsAndSolution",
						 &ug::AssembleLinearOperatorRhsAndSolution<TAlgebra>, grp);
	}

//	some functions
	{
		std::string grp = parentGroup; grp.append("/Algebra/Operation");
		reg.add_function("MatIdentity", &MatIdentity<vector_type, vector_type, matrix_type>, grp,
				"", "opOut", "sets matrix to identity");
		reg.add_function("MatAdd", &MatAdd<vector_type, vector_type, matrix_type>, grp,
				"", "res#alpha1#A1#alpha2#A2", "calculates res = alpha1*A1 + alpha2*A2");
		reg.add_function("MatScale", &MatScale<vector_type, vector_type, matrix_type>, grp,
				"", "mat#alpha", "calculates mat = mat*alpha");
		reg.add_function("MatTranspose", &MatTranspose<vector_type, vector_type, matrix_type>, grp,
				"", "matTranspose#mat", "calculates matTranspose = mat^T");
	}

//	ILineSearch
	{
		std::string grp = parentGroup; grp.append("/Discretization/Nonlinear");
		typedef ILineSearch<vector_type> T;
		string name = string("ILineSearch").append(suffix);
		reg.add_class_<T>(name, grp);
		reg.add_class_to_group(name, "ILineSearch", tag);
	}

//	StandardLineSearch
	{
		std::string grp = parentGroup; grp.append("/Discretization/Nonlinear");
		typedef StandardLineSearch<vector_type> T;
		typedef ILineSearch<vector_type> TBase;
		string name = string("StandardLineSearch").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.template add_constructor<void (*)(int, number, number, bool)>()
			.template add_constructor<void (*)(int, number, number, bool, bool)>()
			.add_method("set_maximum_steps", &T::set_maximum_steps, "Max. steps of the line search", "steps")
			.add_method("set_lambda_start", &T::set_lambda_start, "Initial step scaling factor", "start")
			.add_method("set_reduce_factor", &T::set_reduce_factor, "Step reduction factor", "factor")
			.add_method("set_suff_descent_factor", &T::set_suff_descent_factor, "Sufficient descent", "factor")
			.add_method("set_accept_best", &T::set_accept_best, "If to accept the best", "bAcceptBest")
			.add_method("set_check_all", &T::set_check_all, "If to test all the steps", "bCheckAll")
			.add_method("set_maximum_defect", &T::set_maximum_defect, "Max. norm of the defect", "maxDef")
			.add_method("set_verbose", &T::set_verbose, "Verbosity level", "verboseLevel")
			.add_method("set_offset", &T::set_offset, "Indentation", "strOffset")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "StandardLineSearch", tag);
	}

// PreviousSolutions
	{
		std::string grp = parentGroup; grp.append("/Discretization/TimeDisc");
		string name = string("SolutionTimeSeries").append(suffix);
		typedef VectorTimeSeries<vector_type> T;
		reg.add_class_<T>(name, grp)
			.add_constructor()
			.add_method("size", &T::size, "number of time steps handled")
			.add_method("push_discard_oldest", &T::push_discard_oldest, "oldest solution", "vec#time", "adds new time point, oldest solution is discarded and returned")
			.add_method("push", &T::push, "", "vec#time", "adds new time point, not discarding the oldest")
			.add_method("solution", static_cast<ConstSmartPtr<vector_type> (T::*)(size_t) const>(&T::solution),
					"the local vector for the i'th time point", "i")
			.add_method("oldest", static_cast<SmartPtr<vector_type> (T::*)()>(&T::oldest),
					"oldest solution")
			.add_method("latest", static_cast<SmartPtr<vector_type> (T::*)()>(&T::latest),
					"latest solution")
			.add_method("time", &T::time, "point in time for solution", "i")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "SolutionTimeSeries", tag);
	}
}

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
static void DomainAlgebra(Registry& reg, string parentGroup)
{
	//	typedefs for Vector and Matrix
	typedef typename TAlgebra::vector_type vector_type;

	string suffix = GetDomainAlgebraSuffix<TDomain,TAlgebra>();
	string tag = GetDomainAlgebraTag<TDomain,TAlgebra>();

	// 	CompositeConvCheck
	{
		std::string &grp = parentGroup;
		typedef CompositeConvCheck<vector_type, TDomain> T;
		typedef IConvergenceCheck<vector_type> TBase;
		string name = string("CompositeConvCheck").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(SmartPtr<ApproximationSpace<TDomain> >)>("ApproximationSpace")
			.template add_constructor<void (*)(SmartPtr<ApproximationSpace<TDomain> >, int, number, number)>("ApproximationSpace#maxSteps#minDefect#relReduction")
			.add_method("set_level", (void (T::*)(int)) &T::set_level,
			            "", "grid_level", "sets grid level where defect vectors come from")
			.add_method("set_time_measurement", &T::set_time_measurement,
						"", "", "whether to perform a time measurement or not")
			.add_method("set_adaptive", &T::set_adaptive,
						"", "", "whether adaptive meshing is used or not")
			.add_method("set_maximum_steps", &T::set_maximum_steps, 
			            "", "maximum steps|default|min=0;value=100")
			.add_method("set_component_check", (void (T::*)(const std::string&,const std::vector<number>&,const std::vector<number>&)) &T::set_component_check,
						"", "Components#minDefect#relReduction")
			.add_method("set_component_check", (void (T::*)(const std::vector<std::string>&,const std::vector<number>&,const std::vector<number>&)) &T::set_component_check,
						"", "Components#minDefect#relReduction")
			.add_method("set_component_check", (void (T::*)(const std::vector<std::string>&,const number,const number)) &T::set_component_check,
						"", "Components#minDefect#relReduction")
			.add_method("set_component_check", (void (T::*)(const std::string&,const number,const number)) &T::set_component_check,
						"", "Components#minDefect#relReduction")
			.add_method("set_all_component_check", (void (T::*)(const number,const number)) &T::set_all_component_check,
						"", "minDefect#relReduction")
			.add_method("set_group_check", (void (T::*)(const std::vector<std::string>&,const number,const number)) &T::set_group_check,
						"", "ComponentGroup#minDefect#relReduction")
			.add_method("set_group_check", (void (T::*)(const std::string&,const number,const number)) &T::set_group_check,
						"", "ComponentGroup#minDefect#relReduction")
			.add_method("disable_rest_check", &T::disable_rest_check,
						"", "")
			.add_method("set_rest_check", &T::set_rest_check,
						"", "minDefect#relReduction")
			.add_method("set_verbose", &T::set_verbose, "", "verbosity")
			.add_method("set_supress_unsuccessful", &T::set_supress_unsuccessful,
						"", "", "whether to always report success when iteration is finished")
			.add_method("defect", (number (T::*)(size_t) const) &T::defect, 
			            "defect", "function index", "returns the current defect")
			.add_method("step", &T::step,  
			            "step", "", "returns the current number of steps")
			.add_method("reduction", (number (T::*)(size_t) const) &T::reduction, 
			            "reduction", "function index", "returns the current relative reduction for a function")
			.add_method("rate", (number (T::*)(size_t) const) &T::rate, 
			            "rate", "function index", "returns the current convergence rate for a function")
			.add_method("avg_rate", (number (T::*)(size_t) const) &T::avg_rate, 
			            "avg_rate", "function index", "returns the averaged convergence rate for a function")
			.add_method("iteration_ended", &T::iteration_ended)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "CompositeConvCheck", tag);
	}

	//	NonlinearGaussSeidelSolver
	{
		std::string grp = parentGroup; grp.append("/Discretization/Nonlinear");
		typedef NLGaussSeidelSolver<TDomain, TAlgebra> T;
		typedef IOperatorInverse<vector_type> TBase;
		typedef DebugWritingObject<TAlgebra> TBase2;
		string name = string("NLGaussSeidelSolver").append(suffix);
		reg.add_class_<T, TBase, TBase2>(name, grp)
			.add_constructor()
			.template add_constructor<void (*)(SmartPtr<ApproximationSpace<TDomain> >,
				SmartPtr<IConvergenceCheck<vector_type> >)> ("ApproxSpaceConvCheck")
			.add_method("set_approximation_space", &T::set_approximation_space, "", "approxSpace")
			.add_method("set_convergence_check", &T::set_convergence_check, "", "convCheck")
			.add_method("set_damp", &T::set_damp, "", "setDampingFactor")
			.add_method("set_constraint", &T::set_constraint, "", "setConstraint")
			.add_method("init", &T::init, "success", "op")
			.add_method("prepare", &T::prepare, "success", "u")
			.add_method("apply", &T::apply, "success", "u")
			.add_method("config_string", &T::config_string)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NLGaussSeidelSolver", tag);
	}

	//	NestedIterationSolver
	{
		std::string grp = parentGroup; grp.append("/Discretization");
		typedef NestedIterationSolver<TDomain, TAlgebra> T;
		typedef IOperatorInverse<vector_type> TBase;
		typedef DebugWritingObject<TAlgebra> TBase2;
		string name = string("NestedIterationSolver").append(suffix);
		reg.add_class_<T, TBase, TBase2>(name, grp)
							.add_constructor()
							.template add_constructor<void (*)(SmartPtr<IOperator<vector_type> >)>("Operator")
							.template add_constructor<void (*)(SmartPtr<IAssemble<TAlgebra> >)>("AssemblingRoutine")
							.template add_constructor<void (*)(SmartPtr<IAssemble<TAlgebra> >, SmartPtr<IAssemble<TAlgebra> >)>("AssemblingRoutine#ErrorRoutine")
							.add_method("set_linear_solver", &T::set_linear_solver, "", "linSolver")
						// 	.add_method("set_convergence_check", &T::set_convergence_check, "", "convCheck")
							.add_method("init", &T::init, "success", "op")
							.add_method("prepare", &T::prepare, "success", "u")
							.add_method("apply", &T::apply, "success", "u")
							.add_method("top_level", &T::top_level)
							.add_method("set_top_level", &T::set_top_level)
							.add_method("base_level", &T::base_level)
							.add_method("set_base_level", &T::set_base_level)
							.add_method("set_refiner", &T::set_refiner)
							.add_method("set_refinement_marking", &T::set_refinement_marking)
							.add_method("set_tolerance", &T::set_tolerance)
							.add_method("set_max_steps", &T::set_max_steps)
							.add_method("last_error", &T::last_error)
							.add_method("set_debug_elem_error", &T::set_debug_elem_error)

							.add_method("use_adaptive_refinement", &T::use_adaptive_refinement)
							.add_method("enable_adaptive_refinement", &T::enable_adaptive_refinement)
							.add_method("disable_adaptive_refinement", &T::disable_adaptive_refinement)
							.add_method("config_string", &T::config_string)
							.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NestedIterationSolver", tag);
	}

	//	LagrangeMultiplierDisc base class
	typedef GridFunction<TDomain, TAlgebra> function_type;
	{
		std::string &grp = parentGroup;
		typedef ILagrangeMultiplierDisc<TDomain, function_type> T;
		string name = string("ILagrangeMultiplierDisc").append(suffix);
		reg.add_class_<T>(name, grp)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ILagrangeMultiplierDisc", tag);
	}

	//	ActiveSet
	{
		std::string &grp = parentGroup;
		typedef ActiveSet<TDomain, TAlgebra> T;
		string name = string("ActiveSet").append(suffix);
		reg.add_class_<T>(name, grp)
			.add_constructor()
			.add_method("set_obstacle", &T::set_obstacle, "", "sets limiting "
					"obstacle constraint")
			.add_method("set_lagrange_multiplier_disc", &T::set_lagrange_multiplier_disc, "",
					"discretization to compute the lagrange multiplier")
			.add_method("prepare", &T::prepare, "", "prepare")
			.add_method("active_index", &T::active_index, "", "",
					"is index active or not, stores activeSetList")
			.add_method("set_dirichlet_rows", &T::set_dirichlet_rows, "", "",
					"sets dirichlet rows for active DoFs")
			.add_method("lagrange_multiplier", &T::lagrange_multiplier, "", "",
					"computes lagrange multiplier")
			.add_method("residual_lagrange_mult", &T::residual_lagrange_mult, "", "",
					"computes lagrange multiplier")
			.add_method("check_conv", &T::check_conv, "", "",
					"activeIndexSet changed or not")
			.add_method("check_ineq", &T::check_inequ, "", "",
						"check if all inequalities are fulfilled")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ActiveSet", tag);
	}

}

static void Common(Registry& reg, string parentGroup)
{
	std::string grp = parentGroup; grp.append("/Discretization/Nonlinear");
	//  INewtonUpdate
	{
		typedef INewtonUpdate T;
		string name = string("INewtonUpdate");
		reg.add_class_<T>(name, grp);
	}
}

}; // end Functionality

// end group discalgebra_bridge
/// \}

}// end DiscAlgebra

/// \addtogroup discalgebra_bridge
void RegisterBridge_DiscAlgebra(Registry& reg, string grp)
{
	typedef DiscAlgebra::Functionality Functionality;

	try{
		RegisterAlgebraDependent<Functionality>(reg,grp);
		RegisterDomainAlgebraDependent<Functionality>(reg,grp);
		RegisterCommon<Functionality>(reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

} // end namespace ug
} // end namespace ug
