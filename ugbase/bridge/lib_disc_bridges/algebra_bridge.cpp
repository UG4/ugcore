/*
 * algebra_bridge.cpp
 *
 *  Created on: 22.09.2010
 *      Author: andreasvogel
 */

// extern libraries
#include <iostream>
#include <sstream>
#include <string>

// registry
#include "registry/registry.h"
#include "bridge/bridge.h"

// algebra chooser
#include "lib_algebra/cpu_algebra_types.h"
#include "lib_algebra/operator/matrix_operator_functions.h"

// discretization interfaces
#include "lib_disc/spatial_disc/domain_disc_interface.h"
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"
#include "lib_disc/spatial_disc/constraints/constraint_interface.h"
#include "lib_disc/time_disc/time_disc_interface.h"

// time discretization implementation
#include "lib_disc/time_disc/theta_time_step.h"

// operator interfaces
#include "lib_disc/operator/linear_operator/assembled_linear_operator.h"
#include "lib_disc/operator/non_linear_operator/assembled_non_linear_operator.h"

// newton solver
#include "lib_disc/operator/non_linear_operator/line_search.h"
#include "lib_disc/operator/non_linear_operator/newton_solver/newton.h"

using namespace std;

namespace ug
{
namespace bridge
{

template <typename TAlgebra>
static void Register__Algebra(Registry& reg, string parentGroup)
{
//	typedefs for Vector and Matrix
	typedef typename TAlgebra::vector_type vector_type;
	typedef typename TAlgebra::matrix_type matrix_type;

//	get group string
	string grp = parentGroup; grp.append("/Discretization");

//	suffix and tag
	string algSuffix = GetAlgebraSuffix<TAlgebra>();
	string algTag = GetAlgebraTag<TAlgebra>();

	try{
//	IConstraint
	{
		std::string grp = parentGroup; grp.append("/Discretization/SpatialDisc");
		typedef IConstraint<TAlgebra> T;
		string name = string("IConstraint").append(algSuffix);
		reg.add_class_<T>(name, grp);
		reg.add_class_to_group(name, "IConstraint", algTag);
	}

//	IAssemble
	{
		std::string grp = parentGroup; grp.append("/Discretization/SpatialDisc");
		typedef IAssemble<TAlgebra> T;
		string name = string("IAssemble").append(algSuffix);
		reg.add_class_<T>(name, grp)
			.add_method("assemble_jacobian", static_cast<void (T::*)(matrix_type&, const vector_type&)>(&T::assemble_jacobian),
					"", "J(u)#u", "assembles jacobian on surface grid")
			.add_method("assemble_defect", static_cast<void (T::*)(vector_type&, const vector_type&)>(&T::assemble_defect),
					"", "d(u)#u", "Assembles Defect at a given Solution u.")
			.add_method("assemble_linear", static_cast<void (T::*)(matrix_type&, vector_type&)>(&T::assemble_linear),
					"", "A#b", "Assembles Matrix and rhs on surface grid.")
			.add_method("assemble_stiffness_matrix", static_cast<void (T::*)(matrix_type&, const vector_type&)>(&T::assemble_stiffness_matrix),
					"", "A#u", "assembles stiffness matrix on surface grid")
			.add_method("assemble_mass_matrix", static_cast<void (T::*)(matrix_type&, const vector_type&)>(&T::assemble_mass_matrix),
					"", "M#u", "assembles mass matrix on surface grid")
			.add_method("assemble_rhs", static_cast<void (T::*)(vector_type&, const vector_type&)>(&T::assemble_rhs),
					"", "rhs#u", "assembles right-hand side on surface grid")
			.add_method("adjust_solution", static_cast<void (T::*)(vector_type&)>(&T::adjust_solution));
		reg.add_class_to_group(name, "IAssemble", algTag);
	}

//	IDomainDiscretization
	{
		std::string grp = parentGroup; grp.append("/Discretization/SpatialDisc");
		typedef IAssemble<TAlgebra> TBase;
		typedef IDomainDiscretization<TAlgebra> T;
		string name = string("IDomainDiscretization").append(algSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_method("assemble_jacobian", static_cast<void (T::*)
			            (matrix_type&, ConstSmartPtr<VectorTimeSeries<vector_type> >,
			             number)>(&T::assemble_jacobian));
		reg.add_class_to_group(name, "IDomainDiscretization", algTag);
	}

//	ITimeDiscretization
	{
		std::string grp = parentGroup; grp.append("/Discretization/TimeDisc");
		typedef IAssemble<TAlgebra>  TBase;
		typedef ITimeDiscretization<TAlgebra> T;
		string name = string("ITimeDiscretization").append(algSuffix);
		reg.add_class_<T,TBase>(name, grp)
			.add_method("prepare_step", &T::prepare_step, "", "prepares the assembling of Defect/Jacobian for a time step")
			.add_method("prepare_step_elem", static_cast<void (T::*)(SmartPtr<VectorTimeSeries<vector_type> >, number)>(&T::prepare_step_elem))
			.add_method("finish_step_elem", static_cast<void (T::*)(SmartPtr<VectorTimeSeries<vector_type> >, number)>(&T::finish_step_elem))
			.add_method("num_stages", &T::num_stages, "the number of stages")
			.add_method("set_stage", &T::set_stage, "", "stage")
			.add_method("future_time", &T::future_time, "the future time point (i.e. the one that will be computed)")
			.add_method("num_prev_steps", &T::num_prev_steps);
		reg.add_class_to_group(name, "ITimeDiscretization", algTag);
	}

//	ThetaTimeStep
	{
		std::string grp = parentGroup; grp.append("/Discretization/TimeDisc");
		typedef ITimeDiscretization<TAlgebra> TBase;
		typedef ThetaTimeStep<TAlgebra> T;
		string name = string("ThetaTimeStep").append(algSuffix);
		reg.add_class_<T, TBase>(name, grp)
				.template add_constructor<void (*)(SmartPtr<IDomainDiscretization<TAlgebra> >)>("Domain Discretization")
				.template add_constructor<void (*)(SmartPtr<IDomainDiscretization<TAlgebra> >,number)>("Domain Discretization#Theta")
				.template add_constructor<void (*)(SmartPtr<IDomainDiscretization<TAlgebra> >,const char*)>("Domain Discretization#Scheme")
				.add_method("set_theta", &T::set_theta, "", "Theta (1 = Impl; 0 = Expl)")
				.add_method("set_scheme", &T::set_scheme, "", "Scheme|selection|value=[\"Theta\",\"Alexander\",\"FracStep\"]")
				.add_method("set_stage", &T::set_stage, "", "Stage")
				.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ThetaTimeStep", algTag);
	}

//	BDF
	{
		std::string grp = parentGroup; grp.append("/Discretization/TimeDisc");
		typedef ITimeDiscretization<TAlgebra> TBase;
		typedef BDF<TAlgebra> T;
		string name = string("BDF").append(algSuffix);
		reg.add_class_<T, TBase>(name, grp)
				.template add_constructor<void (*)(SmartPtr<IDomainDiscretization<TAlgebra> >)>("Domain Discretization")
				.template add_constructor<void (*)(SmartPtr<IDomainDiscretization<TAlgebra> >,int)>("Domain Discretization#Order")
				.add_method("set_order", &T::set_order, "", "Order")
				.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "BDF", algTag);
	}

//	AssembledLinearOperator
	{
		typedef AssembledLinearOperator<TAlgebra> T;
		typedef MatrixOperator<matrix_type, vector_type> TBase;
		string name = string("AssembledLinearOperator").append(algSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.template add_constructor<void (*)(IAssemble<TAlgebra>&)>("Assembling Routine")
			.add_method("set_discretization", &T::set_discretization)
			.add_method("set_level", &T::set_level)
			.add_method("set_dirichlet_values", &T::set_dirichlet_values)
			.add_method("init_op_and_rhs", &T::init_op_and_rhs)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "AssembledLinearOperator", algTag);
	}

//	NewtonSolver
	{
		std::string grp = parentGroup; grp.append("/Discretization/Nonlinear");
		typedef NewtonSolver<TAlgebra> T;
		typedef IOperatorInverse<vector_type> TBase;
		typedef DebugWritingObject<TAlgebra> TBase2;
		string name = string("NewtonSolver").append(algSuffix);
		reg.add_class_<T, TBase, TBase2>(name, grp)
			.add_constructor()
			.add_method("set_linear_solver", &T::set_linear_solver, "", "linSolver")
			.add_method("set_convergence_check", &T::set_convergence_check, "", "convCheck")
			.add_method("set_line_search", &T::set_line_search, "", "lineSeach")
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
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NewtonSolver", algTag);
	}

//	AssembledOperator
	{
		std::string grp = parentGroup; grp.append("/Discretization/Nonlinear");
		typedef AssembledOperator<TAlgebra> T;
		typedef IOperator<vector_type> TBase;
		string name = string("AssembledOperator").append(algSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.template add_constructor<void (*)(IAssemble<TAlgebra>&)>("Assembling Routine")
			.add_method("set_discretization", &T::set_discretization, "", "ass")
			.add_method("set_level", &T::set_level, "", "gridLevel")
			.add_method("init", &T::init)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "AssembledOperator", algTag);
	}

//	some functions
	{
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
	}

//	ILineSearch
	{
		std::string grp = parentGroup; grp.append("/Discretization/Nonlinear");
		typedef ILineSearch<vector_type> T;
		string name = string("ILineSearch").append(algSuffix);
		reg.add_class_<T>(name, grp);
		reg.add_class_to_group(name, "ILineSearch", algTag);
	}

//	StandardLineSearch
	{
		std::string grp = parentGroup; grp.append("/Discretization/Nonlinear");
		typedef StandardLineSearch<vector_type> T;
		typedef ILineSearch<vector_type> TBase;
		string name = string("StandardLineSearch").append(algSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.add_method("set_maximum_steps", &T::set_maximum_steps, "", "steps")
			.add_method("set_lambda_start", &T::set_lambda_start, "", "start")
			.add_method("set_reduce_factor", &T::set_reduce_factor, "", "factor")
			.add_method("set_accept_best", &T::set_accept_best, "", "bAcceptBest")
			.add_method("set_maximum_defect", &T::set_maximum_defect, "", "maxDef")
			.add_method("set_verbose", &T::set_verbose, "", "verboseLevel")
			.add_method("set_offset", &T::set_offset, "", "strOffset")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "StandardLineSearch", algTag);
	}

// PreviousSolutions
	{
		std::string grp = parentGroup; grp.append("/Discretization/TimeDisc");
		string name = string("SolutionTimeSeries").append(algSuffix);
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
		reg.add_class_to_group(name, "SolutionTimeSeries", algTag);
	}

	} catch(UG_REGISTRY_ERROR_RegistrationFailed ex)
	{
		UG_LOG("### ERROR in RegisterLibDisc_Algebra: "
				"Registration failed (using name " << ex.name << ").\n");
		UG_THROW_FATAL("Registration failed.");
	}
}

bool RegisterLibDisc_Algebra(Registry& reg, string parentGroup)
{
#ifdef UG_CPU_1	
	Register__Algebra<CPUAlgebra>(reg, parentGroup);
#endif
#ifdef UG_CPU_2
	Register__Algebra<CPUBlockAlgebra<2> >(reg, parentGroup);
#endif
#ifdef UG_CPU_3
	Register__Algebra<CPUBlockAlgebra<3> >(reg, parentGroup);
#endif
#ifdef UG_CPU_4
	Register__Algebra<CPUBlockAlgebra<4> >(reg, parentGroup);
#endif
#ifdef UG_CPU_VAR
	Register__Algebra<CPUVariableBlockAlgebra>(reg, parentGroup);
#endif

	return true;
}

} // end namespace ug
} // end namespace ug
