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
			.add_method("assemble_jacobian", static_cast<void (T::*)(matrix_type&, const vector_type&)>(&T::assemble_jacobian))
			.add_method("assemble_defect", static_cast<void (T::*)(vector_type&, const vector_type&)>(&T::assemble_defect))
			.add_method("assemble_linear", static_cast<void (T::*)(matrix_type&, vector_type&)>(&T::assemble_linear))
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
			            (matrix_type&, const VectorTimeSeries<vector_type>&,
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
			.add_method("prepare_step", &T::prepare_step)
			.add_method("prepare_step_elem", static_cast<void (T::*)(VectorTimeSeries<typename TAlgebra::vector_type>&, number)>(&T::prepare_step_elem))
			.add_method("finish_step_elem", static_cast<void (T::*)(VectorTimeSeries<typename TAlgebra::vector_type>&, number)>(&T::finish_step_elem))
			.add_method("num_stages", &T::num_stages)
			.add_method("set_stage", &T::set_stage)
			.add_method("future_time", &T::future_time)
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
			.add_method("set_linear_solver", &T::set_linear_solver)
			.add_method("set_convergence_check", &T::set_convergence_check)
			.add_method("set_line_search", &T::set_line_search)
			.add_method("init", &T::init)
			.add_method("prepare", &T::prepare)
			.add_method("apply", &T::apply)
			.add_method("print_average_convergence", &T::print_average_convergence)
			.add_method("clear_average_convergence", &T::clear_average_convergence)
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
			.add_method("set_discretization", &T::set_discretization)
			.add_method("set_level", &T::set_level)
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
		reg.add_function("MatIdentity", &MatIdentity<vector_type, vector_type, matrix_type>, grp);
		reg.add_function("MatAdd", &MatAdd<vector_type, vector_type, matrix_type>, grp);
		reg.add_function("MatScale", &MatScale<vector_type, vector_type, matrix_type>, grp);
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
			.add_method("set_maximum_steps", &T::set_maximum_steps)
			.add_method("set_lambda_start", &T::set_lambda_start)
			.add_method("set_reduce_factor", &T::set_reduce_factor)
			.add_method("set_accept_best", &T::set_accept_best)
			.add_method("set_maximum_defect", &T::set_maximum_defect)
			.add_method("set_verbose", &T::set_verbose)
			.add_method("set_offset", &T::set_offset)
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
			.add_method("size", &T::size)
			.add_method("push_discard_oldest", &T::push_discard_oldest)
			.add_method("push", &T::push)
			.add_method("solution", static_cast<const vector_type&(T::*)(size_t) const>(&T::solution))
			.add_method("oldest", static_cast<vector_type& (T::*)()>(&T::oldest))
			.add_method("latest", static_cast<vector_type& (T::*)()>(&T::latest))
			.add_method("time", &T::time);
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
