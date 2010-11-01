/*
 * lib_algebra_bridge.cpp
 *
 *  Created on: 11.10.2010
 *      Author: andreasvogel
 */

#include "../ug_bridge.h"
#include "lib_algebra/lib_algebra.h"
#include "lib_discretization/lib_discretization.h"
#include "lib_algebra/operator/preconditioner/amg/amg.h"
#include <iostream>
#include <sstream>

namespace ug
{
namespace bridge
{

template <typename TAlgebra>
void RegisterAlgebraType(Registry& reg, const char* parentGroup)
{
//	get group string (use same as parent)
	std::string grp = std::string(parentGroup);

//	typedefs for this algebra
	typedef TAlgebra algebra_type;
	typedef typename algebra_type::vector_type vector_type;
	typedef typename algebra_type::matrix_type matrix_type;


	//	Vector
	{
		reg.add_class_<vector_type>("Vector", grp.c_str())
			.add_constructor()
			.add_method("print", &vector_type::p);

		//reg.add_function("VecScaleAssign", , "",
			//	"dest, alpha1, vec1", "dest = alpha1*vec1");
		reg.add_function("VecScaleAdd2", &VecScaleAdd2<vector_type>, "",
				"dest, alpha1, vec1, alpha2, vec2", "dest = alpha1*vec1 + alpha2*vec2");
		reg.add_function("VecScaleAdd3", &VecScaleAdd3<vector_type>, "",
				"dest, alpha1, vec1, alpha2, vec2, alpha3, vec3", "dest = alpha1*vec1 + alpha2*vec2 + alpha3*vec3");
	}

	//	Matrix
	{
		reg.add_class_<matrix_type>("Matrix", grp.c_str())
			.add_constructor()
			.add_method("print", &matrix_type::p);
	}

	// Base Classes
	{
		reg.add_class_<ILinearOperator<vector_type, vector_type> >("ILinearOperator", grp.c_str());
		reg.add_class_<	IMatrixOperator<vector_type, vector_type, matrix_type>,
						ILinearOperator<vector_type, vector_type> >("IMatrixOperator", grp.c_str());

		reg.add_class_<ILinearIterator<vector_type, vector_type> >("ILinearIterator", grp.c_str());
		reg.add_class_<	IPreconditioner<algebra_type>,
						ILinearIterator<vector_type, vector_type> >("IPreconditioner", grp.c_str());

		reg.add_class_<ILinearOperatorInverse<vector_type, vector_type> >("ILinearOperatorInverse", grp.c_str());

		reg.add_class_<IOperator<vector_type, vector_type> >("IOperator", grp.c_str());
		reg.add_class_<IOperatorInverse<vector_type, vector_type> >("IOperatorInverse", grp.c_str());

		reg.add_class_<IProlongationOperator<vector_type, vector_type> >("IProlongationOperator", grp.c_str());
		reg.add_class_<IProjectionOperator<vector_type, vector_type> >("IProjectionOperator", grp.c_str());
	}

	// Preconditioner
	{
	//	Jacobi
		reg.add_class_<	JacobiPreconditioner<algebra_type>,
						IPreconditioner<algebra_type> >("JacobiPreconditioner", grp.c_str())
			.add_constructor()
			.add_method("set_damp", &JacobiPreconditioner<algebra_type>::set_damp, "", "damp");

	//	GaussSeidel
		reg.add_class_<	GSPreconditioner<algebra_type>,
						IPreconditioner<algebra_type> >("GSPreconditioner", grp.c_str())
			.add_constructor();

	//	Symmetric GaussSeidel
		reg.add_class_<	SGSPreconditioner<algebra_type>,
						IPreconditioner<algebra_type> >("SGSPreconditioner", grp.c_str())
			.add_constructor();

	//	Backward GaussSeidel
		reg.add_class_<	BGSPreconditioner<algebra_type>,
						IPreconditioner<algebra_type> >("BGSPreconditioner", grp.c_str())
			.add_constructor();

	//	ILU
		reg.add_class_<	ILUPreconditioner<algebra_type>,
						IPreconditioner<algebra_type> >("ILUPreconditioner", grp.c_str())
			.add_constructor();

	//	ILU Threshold
/*		reg.add_class_<	ILUTPreconditioner<algebra_type>,
						IPreconditioner<algebra_type> >("ILUTPreconditioner", grp.c_str())
			.add_constructor();
*/
	//	AMG
		typedef Domain<2, MultiGrid, MGSubsetHandler> domain_type;
		typedef P1ConformDoFDistribution dof_distribution_type;

	#ifdef UG_PARALLEL
		typedef ParallelGridFunction<GridFunction<domain_type, dof_distribution_type, TAlgebra> > function_type;
	#else
		typedef GridFunction<domain_type, dof_distribution_type, TAlgebra> function_type;
	#endif

//todo: existance of AMGPreconditioner class should not depend on defines.
	#ifdef LAPACK_AVAILABLE
	#ifdef BLAS_AVAILABLE
		reg.add_class_<	amg<algebra_type>, IPreconditioner<algebra_type> > ("AMGPreconditioner", grp.c_str())
			.add_constructor()

			.add_method("set_nu1", &amg<algebra_type>::set_nu1, "", "nu1", "sets nr. of presmoothing steps")
			.add_method("set_nu2", &amg<algebra_type>::set_nu2, "", "nu2", "sets nr. of postsmoothing steps")
			.add_method("set_gamma", &amg<algebra_type>::set_gamma, "", "gamma", "sets gamma in multigrid cycle")
			.add_method("set_theta", &amg<algebra_type>::set_theta, "", "theta")
			.add_method("set_max_levels", &amg<algebra_type>::set_max_levels, "", "max_levels", "sets max nr of AMG levels")

			.add_method("tostring", &amg<algebra_type>::tostring)

			.add_method("set_aggressive_coarsening_A_2", &amg<algebra_type>::set_aggressive_coarsening_A_2)
			.add_method("set_aggressive_coarsening_A_1", &amg<algebra_type>::set_aggressive_coarsening_A_1)

			.add_method("set_presmoother", &amg<algebra_type>::set_presmoother, "", "presmoother")
			.add_method("set_postsmoother", &amg<algebra_type>::set_postsmoother, "", "postsmoother")
			.add_method("set_base_solver", &amg<algebra_type>::set_base_solver, "", "basesmoother");

//			// todo: I comment this out, since grid function das not yet been registered to registery at this point.
//			.add_method("set_debug", &amg<algebra_type>::set_debug<function_type>, "", "u",
//					"sets the internal positions of each node");
	#endif
	#endif

	}


	// todo: Solvers should be independent of type and placed in general section
	{
	// 	LinearSolver
		reg.add_class_<	LinearSolver<algebra_type>,
						ILinearOperatorInverse<vector_type, vector_type> >("LinearSolver", grp.c_str())
			.add_constructor()
			.add_method("set_preconditioner", &LinearSolver<algebra_type>::set_preconditioner, "", "preconditioner")
			.add_method("set_convergence_check", &LinearSolver<algebra_type>::set_convergence_check, "", "check");

	// 	CG Solver
		reg.add_class_<	CGSolver<algebra_type>,
						ILinearOperatorInverse<vector_type, vector_type> >("CGSolver", grp.c_str())
			.add_constructor()
			.add_method("set_preconditioner", &CGSolver<algebra_type>::set_preconditioner, "", "preconditioner")
			.add_method("set_convergence_check", &CGSolver<algebra_type>::set_convergence_check, "", "check");

	// 	BiCGStab Solver
		reg.add_class_<	BiCGStabSolver<algebra_type>,
						ILinearOperatorInverse<vector_type, vector_type> >("BiCGStabSolver", grp.c_str())
			.add_constructor()
			.add_method("set_preconditioner", &BiCGStabSolver<algebra_type>::set_preconditioner, "", "preconditioner")
			.add_method("set_convergence_check", &BiCGStabSolver<algebra_type>::set_convergence_check, "", "check");

//todo: existance of LapackLUSolver class should not depend on defines.
	#ifdef LAPACK_AVAILABLE
	#ifdef BLAS_AVAILABLE
	// 	BiCGStab Solver
		reg.add_class_<	LapackLUSolver<algebra_type>,
						ILinearOperatorInverse<vector_type, vector_type> >("LapackLUSolver", grp.c_str())
			.add_constructor()
			.add_method("set_convergence_check", &LapackLUSolver<algebra_type>::set_convergence_check, "", "check");
	#endif
	#endif
	}
}


bool RegisterLibAlgebraInterface(Registry& reg, const char* parentGroup)
{
	try
	{
	//	get group string
		std::stringstream groupString; groupString << parentGroup << "/Algebra";
		std::string grp = groupString.str();

		UG_LOG("RegisterLibAlgebraInterface\n");
		// StandardConvCheck
		{
			reg.add_class_<IConvergenceCheck>("IConvergenceCheck", grp.c_str());

			reg.add_class_<StandardConvCheck, IConvergenceCheck>("StandardConvergenceCheck", grp.c_str())
				.add_constructor()
				.add_method("set_maximum_steps", &StandardConvCheck::set_maxiumum_steps,
						"", "maximum_steps")
				.add_method("set_minimum_defect", &StandardConvCheck::set_minimum_defect,
						"", "minimum_defect")
				.add_method("set_reduction", &StandardConvCheck::set_reduction,
						"", "reduction")
				.add_method("set_verbose_level", &StandardConvCheck::set_verbose_level,
						"", "verbose_level");
		}

	// register algebra
		RegisterAlgebraType<CPUAlgebra>(reg, grp.c_str());
	//	RegisterAlgebraType<Block2x2Algebra>(reg, grp.c_str());

	}
	catch(UG_REGISTRY_ERROR_RegistrationFailed ex)
	{
		UG_LOG("### ERROR in RegisterLibAlgebraInterface: "
				"Registration failed (using name " << ex.name << ").\n");
		return false;
	}

	return true;
}

} // end namespace bridge
} // end namespace ug
