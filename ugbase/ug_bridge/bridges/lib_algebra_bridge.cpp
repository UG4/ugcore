/*
 * lib_algebra_bridge.cpp
 *
 *  Created on: 11.10.2010
 *      Author: andreasvogel
 */

#include "../ug_bridge.h"
#include "lib_algebra/lib_algebra.h"
#include "lib_discretization/lib_discretization.h"


namespace ug
{
namespace bridge
{

template <typename TAlgebra>
void RegisterAlgebraType(Registry& reg)
{
//	typedefs for this algebra
	typedef TAlgebra algebra_type;
	typedef typename algebra_type::vector_type vector_type;
	typedef typename algebra_type::matrix_type matrix_type;


	//	Vector
	{
		reg.add_class_<vector_type>("Vector")
			.add_constructor();
	}

	//	Matrix
	{
		reg.add_class_<matrix_type>("Matrix")
			.add_constructor();
	}

	// Base Classes
	{
		reg.add_class_<ILinearOperator<vector_type, vector_type> >("ILinearOperator");
		reg.add_class_<	IMatrixOperator<vector_type, vector_type, matrix_type>,
						ILinearOperator<vector_type, vector_type> >("IMatrixOperator");

		reg.add_class_<ILinearIterator<vector_type, vector_type> >("ILinearIterator");
		reg.add_class_<	IPreconditioner<algebra_type>,
						ILinearIterator<vector_type, vector_type> >("IPreconditioner");

		reg.add_class_<ILinearOperatorInverse<vector_type, vector_type> >("ILinearOperatorInverse");

		reg.add_class_<IOperator<vector_type, vector_type> >("IOperator");
		reg.add_class_<IOperatorInverse<vector_type, vector_type> >("IOperatorInverse");
	}

	// Preconditioner
	{
	//	Jacobi
		reg.add_class_<	JacobiPreconditioner<algebra_type>,
						IPreconditioner<algebra_type> >("JacobiPreconditioner")
			.add_constructor()
			.add_method("set_damp", &JacobiPreconditioner<algebra_type>::set_damp);

	//	GaussSeidel
		reg.add_class_<	GSPreconditioner<algebra_type>,
						IPreconditioner<algebra_type> >("GSPreconditioner")
			.add_constructor();

	//	Symmetric GaussSeidel
		reg.add_class_<	SGSPreconditioner<algebra_type>,
						IPreconditioner<algebra_type> >("SGSPreconditioner")
			.add_constructor();

	//	Backward GaussSeidel
		reg.add_class_<	BGSPreconditioner<algebra_type>,
						IPreconditioner<algebra_type> >("BGSPreconditioner")
			.add_constructor();

	//	ILU
		reg.add_class_<	ILUPreconditioner<algebra_type>,
						IPreconditioner<algebra_type> >("ILUPreconditioner")
			.add_constructor();

	//	ILU Threshold
		reg.add_class_<	ILUTPreconditioner<algebra_type>,
						IPreconditioner<algebra_type> >("ILUTPreconditioner")
			.add_constructor();

	//	AMG

		typedef Domain<2, MultiGrid, MGSubsetHandler> domain_type;
		typedef P1ConformDoFDistribution dof_distribution_type;
		typedef MartinAlgebra algebra_type;

	#ifdef UG_PARALLEL
		typedef ParallelGridFunction<GridFunction<domain_type, dof_distribution_type, TAlgebra> > function_type;
	#else
		typedef GridFunction<domain_type, dof_distribution_type, TAlgebra> function_type;
	#endif

		reg.add_class_<	amg<algebra_type>, IPreconditioner<algebra_type> > ("AMGPreconditioner")
			.add_constructor()

			.add_method("set_nu1", &amg<algebra_type>::set_nu1)
			.add_method("set_nu2", &amg<algebra_type>::set_nu2)
			.add_method("set_gamma", &amg<algebra_type>::set_gamma)
			.add_method("set_theta", &amg<algebra_type>::set_theta)
			.add_method("set_max_levels", &amg<algebra_type>::set_max_levels)

			.add_method("tostring", &amg<algebra_type>::tostring)

			.add_method("set_aggressive_coarsening_A_2", &amg<algebra_type>::set_aggressive_coarsening_A_2)
			.add_method("set_aggressive_coarsening_A_1", &amg<algebra_type>::set_aggressive_coarsening_A_1)

			.add_method("set_presmoother", &amg<algebra_type>::set_presmoother)
			.add_method("set_postsmoother", &amg<algebra_type>::set_postsmoother)
			.add_method("set_base_solver", &amg<algebra_type>::set_base_solver)

			.add_method("set_debug", &amg<algebra_type>::set_debug<function_type>);
	}


	// todo: Solvers should be independent of type and placed in general section
	{
	// 	LinearSolver
		reg.add_class_<	LinearSolver<algebra_type>,
						ILinearOperatorInverse<vector_type, vector_type> >("LinearSolver")
			.add_constructor()
			.add_method("set_preconditioner", &LinearSolver<algebra_type>::set_preconditioner)
			.add_method("set_convergence_check", &LinearSolver<algebra_type>::set_convergence_check);

	// 	CG Solver
		reg.add_class_<	CGSolver<algebra_type>,
						ILinearOperatorInverse<vector_type, vector_type> >("CGSolver")
			.add_constructor()
			.add_method("set_preconditioner", &CGSolver<algebra_type>::set_preconditioner)
			.add_method("set_convergence_check", &CGSolver<algebra_type>::set_convergence_check);

	// 	BiCGStab Solver
		reg.add_class_<	BiCGStabSolver<algebra_type>,
						ILinearOperatorInverse<vector_type, vector_type> >("BiCGStabSolver")
			.add_constructor()
			.add_method("set_preconditioner", &BiCGStabSolver<algebra_type>::set_preconditioner)
			.add_method("set_convergence_check", &BiCGStabSolver<algebra_type>::set_convergence_check);

	// 	BiCGStab Solver
		reg.add_class_<	LapackLUSolver<algebra_type>,
						ILinearOperatorInverse<vector_type, vector_type> >("LapackLUSolver")
			.add_constructor();
	}
}


void RegisterLibAlgebraInterface(Registry& reg)
{

	// StandardConvCheck
	{
		reg.add_class_<IConvergenceCheck>("IConvergenceCheck");

		reg.add_class_<StandardConvCheck, IConvergenceCheck>("StandardConvergenceCheck")
			.add_constructor()
			.add_method("set_maximum_steps", &StandardConvCheck::set_maxiumum_steps)
			.add_method("set_minimum_defect", &StandardConvCheck::set_minimum_defect)
			.add_method("set_reduction", &StandardConvCheck::set_reduction)
			.add_method("set_verbose_level", &StandardConvCheck::set_verbose_level);
	}

	// register martin algebra
	RegisterAlgebraType<MartinAlgebra>(reg);

}

} // end namespace bridge
} // end namespace ug
