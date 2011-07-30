/*
 * lib_disc_bridge_domain_independent.cpp
 *
 *  Created on: 22.09.2010
 *      Author: andreasvogel
 */

// extern libraries
#include <iostream>
#include <sstream>

// registry
#include "registry/registry.h"
#include "ug_bridge/ug_bridge.h"

// algebra chooser
#include "lib_algebra/algebra_selector.h"
#include "lib_algebra/algebra_types.h"
#include "lib_algebra/operator/matrix_operator_functions.h"

// lib_discretization part
#include "lib_discretization/dof_manager/dof_distribution.h"
#include "lib_discretization/dof_manager/p1conform/p1conform.h"
#include "lib_discretization/dof_manager/conform/conform.h"

// discretization interfaces
#include "lib_discretization/spatial_discretization/domain_discretization_interface.h"
#include "lib_discretization/spatial_discretization/elem_disc/elem_disc_interface.h"
#include "lib_discretization/spatial_discretization/constraints/constraint_interface.h"
#include "lib_discretization/time_discretization/time_discretization_interface.h"

// time discretization implementation
#include "lib_discretization/time_discretization/theta_time_step.h"

// domain discretization implementation
#include "lib_discretization/spatial_discretization/domain_discretization.h"

// post processes
#include "lib_discretization/spatial_discretization/constraints/continuity_constraints/p1_continuity_constraints.h"

// operator interfaces
#include "lib_discretization/operator/linear_operator/assembled_linear_operator.h"
#include "lib_discretization/operator/non_linear_operator/assembled_non_linear_operator.h"

// newton solver
#include "lib_discretization/operator/non_linear_operator/line_search.h"
#include "lib_discretization/operator/non_linear_operator/newton_solver/newton.h"

namespace ug
{
extern enum_AlgebraType g_AlgebraType;

namespace bridge
{

template <typename TAlgebra, typename TDoFDistribution>
static bool RegisterLibDiscAlgebra__Algebra_DoFDistribution(Registry& reg, string parentGroup)
{
//	typedefs for Vector and Matrix
	typedef typename TAlgebra::vector_type vector_type;
	typedef typename TAlgebra::matrix_type matrix_type;

//	get group string
	string grp = parentGroup; grp.append("/Discretization");

//	suffix and tag
	string algDDSuffix = GetAlgebraSuffix<TAlgebra>();
	algDDSuffix.append(GetDoFDistributionSuffix<TDoFDistribution>());
	string algDDTag = GetAlgebraTag<TAlgebra>();
	algDDTag.append(GetDoFDistributionTag<TDoFDistribution>());

	try{

//	Base class
	{
		typedef IConstraint<TDoFDistribution, TAlgebra> T;
		string name = string("IConstraint").append(algDDSuffix);
		reg.add_class_<T>(name, grp);
		reg.add_class_to_group(name, "IConstraint", algDDTag);
	}

//	OneSideP1ConstraintsPostProcess
	{
		typedef OneSideP1ConstraintsPostProcess<TDoFDistribution, TAlgebra> T;
		typedef IConstraint<TDoFDistribution, TAlgebra> baseT;
		string name = string("OneSideP1Constraints").append(algDDSuffix);
		reg.add_class_<T, baseT>(name, grp)
			.add_constructor();
		reg.add_class_to_group(name, "OneSideP1Constraints", algDDTag);
	}

//	SymP1ConstraintsPostProcess
	{
		typedef SymP1ConstraintsPostProcess<TDoFDistribution, TAlgebra> T;
		typedef IConstraint<TDoFDistribution, TAlgebra> baseT;
		string name = string("SymP1Constraints").append(algDDSuffix);
		reg.add_class_<T, baseT>(name, grp)
			.add_constructor();
		reg.add_class_to_group(name, "SymP1Constraints", algDDTag);
	}

//	IAssemble
	{
		typedef IAssemble<TDoFDistribution, TAlgebra> T;
		string name = string("IAssemble").append(algDDSuffix);
		reg.add_class_<T>(name, grp)
			.add_method("assemble_jacobian", static_cast<bool (T::*)(matrix_type&, const vector_type&,
								const IDoFDistribution<TDoFDistribution>&)>(&T::assemble_jacobian));
		reg.add_class_to_group(name, "IAssemble", algDDTag);
	}

//	IDomainDiscretization
	{
		typedef IAssemble<TDoFDistribution, TAlgebra> TBase;
		typedef IDomainDiscretization<TDoFDistribution, TAlgebra> T;
		string name = string("IDomainDiscretization").append(algDDSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_method("assemble_jacobian", static_cast<bool (T::*)(matrix_type&, const vector_type&,
							number, const SolutionTimeSeries<vector_type>&, const IDoFDistribution<TDoFDistribution>&, number, number)>(&T::assemble_jacobian));
		reg.add_class_to_group(name, "IDomainDiscretization", algDDTag);
	}

//	DomainDiscretization
	{
		typedef IDomainDiscretization<TDoFDistribution, TAlgebra> TBase;
		typedef DomainDiscretization<TDoFDistribution, TAlgebra> T;
		string name = string("DomainDiscretization").append(algDDSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.add_method("add|interactive=false", static_cast<bool (T::*)(IConstraint<TDoFDistribution, TAlgebra>&)>(&T::add),
						"", "Post Process")
			.add_method("add|interactive=false", static_cast<bool (T::*)(IElemDisc&)>(&T::add),
						"", "Discretization")
			.add_method("assemble_mass_matrix", &T::assemble_mass_matrix)
			.add_method("assemble_stiffness_matrix", &T::assemble_stiffness_matrix)
			.add_method("assemble_rhs", &T::assemble_rhs);
		reg.add_class_to_group(name, "DomainDiscretization", algDDTag);
	}

//	ITimeDiscretization
	{
		typedef IAssemble<TDoFDistribution, TAlgebra>  TBase;
		typedef ITimeDiscretization<TDoFDistribution, TAlgebra> T;
		string name = string("ITimeDiscretization").append(algDDSuffix);
		reg.add_class_<T,TBase>(name, grp)
			.add_method("prepare_step", &T::prepare_step)
			.add_method("num_prev_steps", &T::num_prev_steps);
		reg.add_class_to_group(name, "ITimeDiscretization", algDDTag);
	}

//	ThetaTimeDiscretization
	{
		typedef ITimeDiscretization<TDoFDistribution, TAlgebra> TBase;
		typedef ThetaTimeDiscretization<TDoFDistribution, TAlgebra> T;
		string name = string("ThetaTimeDiscretization").append(algDDSuffix);
		reg.add_class_<T, TBase>(name, grp)
				.add_constructor()
				.add_method("set_domain_discretization|interactive=false", &T::set_domain_discretization,
							"", "Domain Discretization")
				.add_method("set_theta|interactive=false", &T::set_theta,
							"", "Theta (0.0 = Impl; 1.0 = Expl)")
				.add_method("prepare_step", &T::prepare_step)
				.add_method("num_prev_steps", &T::num_prev_steps);
		reg.add_class_to_group(name, "ThetaTimeDiscretization", algDDTag);
	}

//	AssembledLinearOperator
	{
		typedef AssembledLinearOperator<TDoFDistribution, TAlgebra> T;
		typedef MatrixOperator<vector_type, vector_type, matrix_type> TBase;
		string name = string("AssembledLinearOperator").append(algDDSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.add_method("set_discretization", &T::set_discretization)
			.add_method("set_dof_distribution", &T::set_dof_distribution)
			.add_method("set_dirichlet_values", &T::set_dirichlet_values)
			.add_method("init_op_and_rhs", &T::init_op_and_rhs);
		reg.add_class_to_group(name, "AssembledLinearOperator", algDDTag);
	}

//	NewtonSolver
	{
		typedef NewtonSolver<TDoFDistribution, TAlgebra> T;
		typedef IOperatorInverse<vector_type, vector_type> TBase;
		string name = string("NewtonSolver").append(algDDSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.add_method("set_linear_solver", &T::set_linear_solver)
			.add_method("set_convergence_check", &T::set_convergence_check)
			.add_method("set_line_search", &T::set_line_search)
			.add_method("init", &T::init)
			.add_method("prepare", &T::prepare)
			.add_method("apply", &T::apply)
			.add_method("set_debug", &T::set_debug);
		reg.add_class_to_group(name, "NewtonSolver", algDDTag);
	}

//	AssembledOperator
	{
		typedef AssembledOperator<TDoFDistribution, TAlgebra> T;
		typedef IOperator<vector_type, vector_type> TBase;
		string name = string("AssembledOperator").append(algDDSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.add_method("set_discretization", &T::set_discretization)
			.add_method("set_dof_distribution", &T::set_dof_distribution)
			.add_method("init", &T::init);
		reg.add_class_to_group(name, "AssembledOperator", algDDTag);
	}

//	some functions
	{
		reg.add_function("AssembleLinearOperatorRhsAndSolution",
						 &AssembleLinearOperatorRhsAndSolution<TDoFDistribution, TAlgebra>);
	}

	} catch(UG_REGISTRY_ERROR_RegistrationFailed ex)
	{
		UG_LOG("### ERROR in RegisterLibDiscAlgebra__Algebra_DoFDistribution: "
				"Registration failed (using name " << ex.name << ").\n");
		return false;
	}

	return true;
}


template <typename TDoFDistribution>
static bool RegisterLibDiscAlgebra__DoFDistribution(Registry& reg, string parentGroup)
{
//	get group string
	std::string grp = parentGroup; grp.append("/Discretization");

//	TDoFDistribution
	{
		typedef IDoFDistribution<TDoFDistribution> T;
		reg.add_class_<T>(GetDoFDistributionSuffix<TDoFDistribution>(), grp);
	}

	return true;
}

template <typename TAlgebra>
static bool RegisterLibDiscAlgebra__Algebra(Registry& reg, string parentGroup)
{
//	typedefs for Vector and Matrix
	typedef typename TAlgebra::vector_type vector_type;
	typedef typename TAlgebra::matrix_type matrix_type;

//	get group string
	std::string grp = parentGroup; grp.append("/Discretization");

//	suffix and tag
	string algSuffix = GetAlgebraSuffix<TAlgebra>();
	string algTag = GetAlgebraTag<TAlgebra>();

	try{

//	some functions
	{
		//reg.add_function("MatAdd", &MatAdd<vector_type, vector_type, matrix_type>);
		reg.add_function("MatIdentity", &MatIdentity<vector_type, vector_type, matrix_type>);
		reg.add_function("MatAdd", &MatAdd<vector_type, vector_type, matrix_type>);
		reg.add_function("MatScale", &MatScale<vector_type, vector_type, matrix_type>);
	}

//	ILineSearch
	{
		typedef ILineSearch<vector_type> T;
		string name = string("ILineSearch").append(algSuffix);
		reg.add_class_<T>(name, grp);
		reg.add_class_to_group(name, "ILineSearch", algTag);
	}

//	StandardLineSearch
	{
		typedef StandardLineSearch<vector_type> T;
		typedef ILineSearch<vector_type> TBase;
		string name = string("StandardLineSearch").append(algSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.add_method("set_maximum_steps", &T::set_maximum_steps)
			.add_method("set_lambda_start", &T::set_lambda_start)
			.add_method("set_reduce_factor", &T::set_reduce_factor)
			.add_method("set_accept_best", &T::set_accept_best)
			.add_method("set_verbose_level", &T::set_verbose_level)
			.add_method("set_offset", &T::set_offset);
		reg.add_class_to_group(name, "StandardLineSearch", algTag);
	}

// PreviousSolutions
	{
		string name = string("SolutionTimeSeries").append(algSuffix);
		typedef SolutionTimeSeries<vector_type> T;
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
		UG_LOG("### ERROR in RegisterLibDiscAlgebra__Algebra: "
				"Registration failed (using name " << ex.name << ").\n");
		return false;
	}


//	add DoFDistributionType
	bool bReturn = true;
	bReturn &= RegisterLibDiscAlgebra__Algebra_DoFDistribution<TAlgebra, P1DoFDistribution<false> >(reg, parentGroup);
//	bReturn &= RegisterLibDiscAlgebra__Algebra_DoFDistribution<TAlgebra, P1DoFDistribution<true> >(reg, parentGroup);
//	bReturn &= RegisterLibDiscAlgebra__Algebra_DoFDistribution<TAlgebra, DoFDistribution >(reg, parentGroup);

	return bReturn;
}

bool RegisterLibDisc_Algebra(Registry& reg, string parentGroup)
{
	bool bReturn = true;
	bReturn &= RegisterLibDiscAlgebra__Algebra<CPUAlgebra>(reg, parentGroup);
//	bReturn &= RegisterLibDiscAlgebra__Algebra<CPUBlockAlgebra<2> >(reg, parentGroup);
	bReturn &= RegisterLibDiscAlgebra__Algebra<CPUBlockAlgebra<3> >(reg, parentGroup);
//	bReturn &= RegisterLibDiscAlgebra__Algebra<CPUBlockAlgebra<4> >(reg, parentGroup);
//	bReturn &= RegisterLibDiscAlgebra__Algebra<CPUVariableBlockAlgebra>(reg, parentGroup);

	bReturn &= RegisterLibDiscAlgebra__DoFDistribution<P1DoFDistribution<false> >(reg, parentGroup);
//	bReturn &= RegisterLibDiscAlgebra__DoFDistribution<P1DoFDistribution<true> >(reg, parentGroup);
//	bReturn &= RegisterLibDiscAlgebra__DoFDistribution<DoFDistribution >(reg, parentGroup);

	return bReturn;
}

} // end namespace ug
} // end namespace ug
