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

// algebra chooser
#include "lib_algebra/algebra_selector.h"
#include "lib_algebra/algebra_types.h"
#include "lib_algebra/operator/matrix_operator_functions.h"

// lib_discretization part
#include "lib_discretization/dof_manager/dof_distribution.h"

// discretization interfaces
#include "lib_discretization/spatial_discretization/domain_discretization_interface.h"
#include "lib_discretization/spatial_discretization/elem_disc/elem_disc_interface.h"
#include "lib_discretization/spatial_discretization/post_process/post_process_interface.h"
#include "lib_discretization/time_discretization/time_discretization_interface.h"

// time discretization implementation
#include "lib_discretization/time_discretization/theta_time_step.h"

// domain discretization implementation
#include "lib_discretization/spatial_discretization/domain_discretization.h"

// post processes
#include "lib_discretization/spatial_discretization/post_process/constraints/p1_constraints_post_process.h"

// operator interfaces
#include "lib_discretization/operator/linear_operator/assembled_linear_operator.h"
#include "lib_discretization/operator/non_linear_operator/assembled_non_linear_operator.h"

// newton solver
#include "lib_discretization/operator/non_linear_operator/line_search.h"
#include "lib_discretization/operator/non_linear_operator/newton_solver/newton.h"

#include "lib_discretization/dof_manager/p1conform/p1conform.h"

namespace ug
{
extern enum_AlgebraType g_AlgebraType;

namespace bridge
{

template <typename TAlgebra, typename TDoFDistribution>
bool RegisterLibDiscForAlgebra(Registry& reg, const char* parentGroup)
{
	typedef TDoFDistribution dof_distribution_type;
	typedef TAlgebra algebra_type;
	typedef typename algebra_type::vector_type vector_type;
	typedef typename algebra_type::matrix_type matrix_type;

	try
	{
	//	get group string
		std::string grp = parentGroup; grp.append("/Discretization");

	//	P1ConformDoFDistribution
		{
			typedef IDoFDistribution<TDoFDistribution> T;
			reg.add_class_<T>("P1ConformDoFDistribution", grp.c_str());
		}

	//	Base class
		reg.add_class_<IPostProcess<dof_distribution_type, algebra_type> >("IPostProcess", grp.c_str());

	//	OneSideP1ConstraintsPostProcess
		{
			typedef OneSideP1ConstraintsPostProcess<dof_distribution_type, algebra_type> T;
			typedef IPostProcess<dof_distribution_type, algebra_type> baseT;
			reg.add_class_<T, baseT>("OneSideP1Constraints", grp.c_str())
				.add_constructor();
		}

	//	SymP1ConstraintsPostProcess
		{
			typedef SymP1ConstraintsPostProcess<dof_distribution_type, algebra_type> T;
			typedef IPostProcess<dof_distribution_type, algebra_type> baseT;
			reg.add_class_<T, baseT>("SymP1Constraints", grp.c_str())
				.add_constructor();
		}

	//	DomainDiscretization
		{
			typedef IAssemble<dof_distribution_type, algebra_type> TBase;
			typedef IDomainDiscretization<dof_distribution_type, algebra_type>	TIDomDisc;
			typedef DomainDiscretization<dof_distribution_type, algebra_type> T;

			reg.add_class_<TBase>("IAssemble", grp.c_str())
				.add_method("assemble_jacobian", static_cast<bool (TBase::*)(matrix_type&, const vector_type&,
									const IDoFDistribution<TDoFDistribution>&)>(&TBase::assemble_jacobian));

			reg.add_class_<TIDomDisc, TBase>("IDomainDiscretization", grp.c_str())
//				.add_method("assemble_jacobian", static_cast<bool (TIDomDisc::*)(matrix_type&, const vector_type&,
//												const IDoFDistribution<TDoFDistribution>&)>(&TIDomDisc::assemble_jacobian))
				.add_method("assemble_jacobian", static_cast<bool (TIDomDisc::*)(matrix_type&, const vector_type&,
								number, const SolutionTimeSeries<vector_type>&, const IDoFDistribution<TDoFDistribution>&, number, number)>(&TIDomDisc::assemble_jacobian));


		//	\TODO: There seems to be an error in the Lua-Skript parsing of
		//		   overloaded functions, or an error in the registry process
		//		   when virtual functions are in the game and are implemented
		//		   on several levels of the class hierarchy.
			reg.add_class_<T, TIDomDisc>("DomainDiscretization", grp.c_str())
				.add_constructor()
				.add_method("add_post_process|interactive=false", &T::add_post_process,
							"", "Post Process")
				.add_method("add_elem_disc|interactive=false", (bool (T::*)(IElemDisc&)) &T::add_elem_disc,
							"", "Discretization")
				.add_method("assemble_mass_matrix", &T::assemble_mass_matrix)
				.add_method("assemble_stiffness_matrix", &T::assemble_stiffness_matrix)
				.add_method("assemble_rhs", &T::assemble_rhs);
//				.add_method("assemble_jacobian", (bool (T::*)(matrix_type&, const vector_type&,
//													const IDoFDistribution<TDoFDistribution>&))&T::assemble_jacobian)
//				.add_method("assemble_jacobian", (bool (T::*)(matrix_type&, const vector_type&,
//							number, const SolutionTimeSeries<vector_type>&, const IDoFDistribution<TDoFDistribution>&, number, number )) &T::assemble_jacobian);

		}

	//	ITimeDiscretization
		{
			typedef IAssemble<dof_distribution_type, algebra_type>  TBase;
			typedef ITimeDiscretization<dof_distribution_type, algebra_type> T;
			reg.add_class_<T,TBase>("ITimeDiscretization", grp.c_str())
				.add_method("prepare_step", &T::prepare_step)
				.add_method("num_prev_steps", &T::num_prev_steps);

		}

	//	ThetaTimeDiscretization
		{
			typedef ITimeDiscretization<dof_distribution_type, algebra_type> TBase;
			typedef ThetaTimeDiscretization<dof_distribution_type, algebra_type> T;
			reg.add_class_<T, TBase>("ThetaTimeDiscretization", grp.c_str())
					.add_constructor()
					.add_method("set_domain_discretization|interactive=false", &T::set_domain_discretization,
								"", "Domain Discretization")
					.add_method("set_theta|interactive=false", &T::set_theta,
								"", "Theta (0.0 = Impl; 1.0 = Expl)")
					.add_method("prepare_step", &T::prepare_step)
					.add_method("num_prev_steps", &T::num_prev_steps);
		}

	//	AssembledLinearOperator
		{
			typedef AssembledLinearOperator<dof_distribution_type, algebra_type> T;

			reg.add_class_<T, MatrixOperator<vector_type, vector_type, matrix_type> >
							("AssembledLinearOperator", grp.c_str())
				.add_constructor()
				.add_method("set_discretization", &T::set_discretization)
				.add_method("set_dof_distribution", &T::set_dof_distribution)
				.add_method("set_dirichlet_values", &T::set_dirichlet_values)
				.add_method("init_op_and_rhs", &T::init_op_and_rhs);


			//reg.add_function("MatAdd", &MatAdd<vector_type, vector_type, matrix_type>);
			reg.add_function("MatIdentity", &MatIdentity<vector_type, vector_type, matrix_type>);
			reg.add_function("MatAdd", &MatAdd<vector_type, vector_type, matrix_type>);
			reg.add_function("MatScale", &MatScale<vector_type, vector_type, matrix_type>);

			reg.add_function("AssembleLinearOperatorRhsAndSolution", &AssembleLinearOperatorRhsAndSolution<dof_distribution_type, algebra_type>);
		}


	//	AssembledOperator
		{
			typedef AssembledOperator<dof_distribution_type, algebra_type> T;

			reg.add_class_<T, IOperator<vector_type, vector_type> >
							("AssembledOperator", grp.c_str())
				.add_constructor()
				.add_method("set_discretization", &T::set_discretization)
				.add_method("set_dof_distribution", &T::set_dof_distribution)
				.add_method("init", &T::init);
		}

	//	StandardLineSearch
		{
			typedef StandardLineSearch<vector_type> T;

			reg.add_class_<ILineSearch<vector_type> >("ILineSearch", grp.c_str());

			reg.add_class_<StandardLineSearch<vector_type>,
							ILineSearch<vector_type> >("StandardLineSearch", grp.c_str())
				.add_constructor()
				.add_method("set_maximum_steps", &T::set_maximum_steps)
				.add_method("set_lambda_start", &T::set_lambda_start)
				.add_method("set_reduce_factor", &T::set_reduce_factor)
				.add_method("set_accept_best", &T::set_accept_best)
				.add_method("set_verbose_level", &T::set_verbose_level)
				.add_method("set_offset", &T::set_offset);
		}

	// PreviousSolutions
		{
			typedef SolutionTimeSeries<vector_type> T;
			reg.add_class_<T>("SolutionTimeSeries", grp.c_str())
				.add_constructor()
				.add_method("size", &T::size)
				.add_method("push_discard_oldest", &T::push_discard_oldest)
				.add_method("push", &T::push)
				.add_method("solution", (const vector_type&(T::*)(size_t) const)&T::solution)
				.add_method("oldest", (vector_type& (T::*)()) &T::oldest)
				.add_method("latest", (vector_type& (T::*)()) &T::latest)
				.add_method("time", &T::time);

		}

	//	NewtonSolver
		{
			typedef NewtonSolver<dof_distribution_type, algebra_type> T;

			reg.add_class_<T, IOperatorInverse<vector_type, vector_type> >("NewtonSolver", grp.c_str())
				.add_constructor()
				.add_method("set_linear_solver", &T::set_linear_solver)
				.add_method("set_convergence_check", &T::set_convergence_check)
				.add_method("set_line_search", &T::set_line_search)
				.add_method("init", &T::init)
				.add_method("prepare", &T::prepare)
				.add_method("apply", &T::apply)
				.add_method("set_debug", &T::set_debug);

		}
	}
	catch(UG_REGISTRY_ERROR_RegistrationFailed ex)
	{
		UG_LOG("### ERROR in RegisterLibDiscretizationInterfaceForAlgebra: "
				"Registration failed (using name " << ex.name << ").\n");
		return false;
	}

	return true;
}

bool RegisterDynamicLibDiscAlgebra(Registry& reg, int algebra_type, const char* parentGroup)
{
	bool bReturn = true;

	switch(algebra_type)
	{
	case eCPUAlgebra:		 		bReturn &= RegisterLibDiscForAlgebra<CPUAlgebra, P1DoFDistribution<false> >(reg, parentGroup); break;
//	case eCPUBlockAlgebra2x2: 		bReturn &= RegisterLibDiscForAlgebra<CPUBlockAlgebra<2>, P1DoFDistribution<true> >(reg, parentGroup); break;
	case eCPUBlockAlgebra3x3: 		bReturn &= RegisterLibDiscForAlgebra<CPUBlockAlgebra<3>, P1DoFDistribution<true> >(reg, parentGroup); break;
//	case eCPUBlockAlgebra4x4: 		bReturn &= RegisterLibDiscForAlgebra<CPUBlockAlgebra<4>, P1DoFDistribution<true> >(reg, parentGroup); break;
//	case eCPUVariableBlockAlgebra: 	bReturn &= RegisterLibDiscForAlgebra<CPUVariableBlockAlgebra, P1DoFDistribution<true> >(reg, parentGroup); break;
	default: UG_ASSERT(0, "Unsupported Algebra Type");
				UG_LOG("Unsupported Algebra Type requested.\n");
				return false;
	}

	return bReturn;
}


} // end namespace ug
} // end namespace ug
