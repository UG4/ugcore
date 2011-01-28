/*
 * lib_algebra_bridge.cpp
 *
 *  Created on: 11.10.2010
 *      Author: andreasvogel
 */
// #define UG_USE_AMG // temporary switch until AMG for systems works again

#include "../ug_bridge.h"
#include "lib_algebra/lib_algebra.h"
#include "lib_discretization/lib_discretization.h"
#include <iostream>
#include <sstream>

namespace ug
{
extern enum_AlgebraType g_AlgebraType;
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
			.add_method("set|hide=true", (bool (vector_type::*)(number))&vector_type::set,
									"Success", "Number")
			.add_method("size|hide=true", (size_t (vector_type::*)())&vector_type::size,
									"Size", "")
			.add_method("set_random|hide=true", (bool (vector_type::*)(number))&vector_type::set_random,
									"Success", "Number")
			.add_method("print|hide=true", &vector_type::p);

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
			.add_method("print|hide=true", &matrix_type::p);
	}

	// Debug Writer (abstract base class)
	{
		typedef IDebugWriter<algebra_type> T;
		reg.add_class_<T>("DebugWriter", grp.c_str());
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
		reg.add_class_<IMatrixOperatorInverse<vector_type, vector_type,matrix_type>,
						ILinearOperatorInverse<vector_type, vector_type> >("IMatrixOperatorInverse", grp.c_str());

		reg.add_class_<IOperator<vector_type, vector_type> >("IOperator", grp.c_str());
		reg.add_class_<IOperatorInverse<vector_type, vector_type> >("IOperatorInverse", grp.c_str());

		reg.add_class_<IProlongationOperator<vector_type, vector_type> >("IProlongationOperator", grp.c_str());
		reg.add_class_<IProjectionOperator<vector_type, vector_type> >("IProjectionOperator", grp.c_str());
	}

	// Preconditioner
	{
	//	get group string
		std::stringstream grpSS2; grpSS2 << grp << "/Preconditioner";
		std::string grp2 = grpSS2.str();

	//	Jacobi
		reg.add_class_<	JacobiPreconditioner<algebra_type>,
						IPreconditioner<algebra_type> >("Jacobi", grp2.c_str())
			.add_constructor()
			.add_method("set_damp", &JacobiPreconditioner<algebra_type>::set_damp, "", "damp");

	//	GaussSeidel
		reg.add_class_<	GSPreconditioner<algebra_type>,
						IPreconditioner<algebra_type> >("GaussSeidel", grp2.c_str())
			.add_constructor();

	//	Symmetric GaussSeidel
		reg.add_class_<	SGSPreconditioner<algebra_type>,
						IPreconditioner<algebra_type> >("SymmetricGaussSeidel", grp2.c_str())
			.add_constructor();

	//	Backward GaussSeidel
		reg.add_class_<	BGSPreconditioner<algebra_type>,
						IPreconditioner<algebra_type> >("BackwardGaussSeidel", grp2.c_str())
			.add_constructor();

	//	ILU
		reg.add_class_<	ILUPreconditioner<algebra_type>,
						IPreconditioner<algebra_type> >("ILU", grp2.c_str())
			.add_constructor();

	//	ILU Threshold
		reg.add_class_<	ILUTPreconditioner<algebra_type>,
						IPreconditioner<algebra_type> >("ILUT", grp2.c_str())
			.add_constructor()
			.add_method("set_threshold", &ILUTPreconditioner<algebra_type>::set_threshold,
						"", "threshold", "sets threshold of incomplete LU factorisation"); // added 01122010ih


#ifdef UG_USE_AMG
	//	AMG
		typedef Domain<2, MultiGrid, MGSubsetHandler> domain_type;
		typedef P1ConformDoFDistribution dof_distribution_type;

	#ifdef UG_PARALLEL
		typedef ParallelGridFunction<GridFunction<domain_type, dof_distribution_type, TAlgebra> > function_type;
	#else
		typedef GridFunction<domain_type, dof_distribution_type, TAlgebra> function_type;
	#endif

//todo: existance of AMGPreconditioner class should not depend on defines.
		reg.add_class_<	amg_base<algebra_type>, IPreconditioner<algebra_type> > ("AMGBase", grp.c_str())
			.add_method("set_num_presmooth", &amg_base<algebra_type>::set_num_presmooth, "", "nu1", "sets nr. of presmoothing steps")
			.add_method("set_num_postsmooth", &amg_base<algebra_type>::set_num_postsmooth, "", "nu2", "sets nr. of postsmoothing steps")
			.add_method("set_cycle_type", &amg_base<algebra_type>::set_cycle_type, "", "gamma", "sets cycle type in multigrid cycle")
			.add_method("set_presmoother", &amg_base<algebra_type>::set_presmoother, "", "presmoother")
			.add_method("set_postsmoother", &amg_base<algebra_type>::set_postsmoother, "", "postsmoother")
			.add_method("set_base_solver", &amg_base<algebra_type>::set_base_solver, "", "basesmoother")
			.add_method("set_max_levels", &amg_base<algebra_type>::set_max_levels, "", "max_levels", "sets max nr of AMG levels")
			.add_method("set_debug", (bool (amg_base<algebra_type>::*)(function_type&)) &amg_base<algebra_type>::set_debug,"", "u",
					"sets the internal positions of each node")
			.add_method("check", &amg_base<algebra_type>::check, "", "")
			.add_method("set_max_nodes_for_exact", &amg_base<algebra_type>::set_max_nodes_for_exact, "", "")
			.add_method("set_max_fill_before_exact", &amg_base<algebra_type>::set_max_fill_before_exact, "", "")
			.add_method("set_matrix_write_path", &amg_base<algebra_type>::set_matrix_write_path, "", "")
			.add_method("set_fsmoothing", &amg_base<algebra_type>::set_fsmoothing, "", "");

		reg.add_class_<	amg<algebra_type>, amg_base<algebra_type> > ("AMGPreconditioner", grp.c_str())
			.add_constructor()
			.add_method("set_theta", &amg<algebra_type>::set_theta, "", "theta")
			.add_method("tostring", &amg<algebra_type>::tostring)
			.add_method("enable_aggressive_coarsening_A_2", &amg<algebra_type>::enable_aggressive_coarsening_A_2)
			.add_method("enable_aggressive_coarsening_A_1", &amg<algebra_type>::enable_aggressive_coarsening_A_1)
			.add_method("disable_aggressive_coarsening", &amg<algebra_type>::disable_aggressive_coarsening);

		reg.add_class_<	famg<algebra_type>, amg_base<algebra_type> > ("FAMGPreconditioner", grp.c_str())
			.add_constructor()
			.add_method("tostring", &famg<algebra_type>::tostring)
			.add_method("set_aggressive_coarsening", &famg<algebra_type>::set_aggressive_coarsening)
			.add_method("set_delta", &famg<algebra_type>::set_delta)
			.add_method("set_theta", &famg<algebra_type>::set_theta)
			.add_method("set_damping_for_smoother_in_interpolation_calculation",
					&famg<algebra_type>::set_damping_for_smoother_in_interpolation_calculation)
			.add_method("set_testvector_zero_at_dirichlet", &famg<algebra_type>::set_testvector_zero_at_dirichlet)
			.add_method("set_testvector_damps", &famg<algebra_type>::set_testvector_damps)
			;
#endif
	}
/*
	{
#ifdef LAPACK_AVAILABLE
		string subgroup = grp; // + string("/Preconditioner");

		reg.add_class_<	PINVIT<algebra_type> >("EigenSolver", subgroup.c_str())
			.add_constructor()
			.add_method("add_vector", &PINVIT<algebra_type>::add_vector,
						"", "vector")
			.add_method("set_preconditioner|interactive=false", &PINVIT<algebra_type>::set_preconditioner,
						"", "Preconditioner||invokeOnChange=true")
			.add_method("set_linear_operator_A|interactive=false", &PINVIT<algebra_type>::set_linear_operator_A,
						"", "LinearOperatorA|invokeOnChange=true")
			.add_method("set_linear_operator_B|interactive=false", &PINVIT<algebra_type>::set_linear_operator_B,
						"", "LinearOperatorB|invokeOnChange=true")
			.add_method("set_max_iterations|interactive=false", &PINVIT<algebra_type>::set_max_iterations,
							"", "precision|invokeOnChange=true")
			.add_method("set_precision|interactive=false", &PINVIT<algebra_type>::set_precision,
							"", "precision|invokeOnChange=true")
			.add_method("apply", &PINVIT<algebra_type>::apply);
#endif
	}
*/
	// todo: Solvers should be independent of type and placed in general section
	{
	//	get group string
		std::stringstream grpSS3; grpSS3 << grp << "/Solver";
		std::string grp3 = grpSS3.str();

	// 	LinearSolver
		reg.add_class_<	LinearSolver<algebra_type>,
						ILinearOperatorInverse<vector_type, vector_type> >("LinearSolver", grp3.c_str())
			.add_constructor()
			.add_method("set_preconditioner|interactive=false", &LinearSolver<algebra_type>::set_preconditioner,
						"", "Preconditioner||invokeOnChange=true")
			.add_method("set_convergence_check|interactive=false", &LinearSolver<algebra_type>::set_convergence_check,
						"", "Check||invokeOnChange=true");

	// 	CG Solver
		reg.add_class_<	CGSolver<algebra_type>,
						ILinearOperatorInverse<vector_type, vector_type> >("CG", grp3.c_str())
			.add_constructor()
			.add_method("set_preconditioner|interactive=false", &CGSolver<algebra_type>::set_preconditioner,
						"", "Preconditioner||invokeOnChange=true")
			.add_method("set_convergence_check|interactive=false", &CGSolver<algebra_type>::set_convergence_check,
						"", "Check||invokeOnChange=true");

	// 	BiCGStab Solver
		reg.add_class_<	BiCGStabSolver<algebra_type>,
						ILinearOperatorInverse<vector_type, vector_type> >("BiCGStab", grp3.c_str())
			.add_constructor()
			.add_method("set_preconditioner|interactive=false", &BiCGStabSolver<algebra_type>::set_preconditioner,
						"", "Preconditioner||invokeOnChange=true")
			.add_method("set_convergence_check|interactive=false", &BiCGStabSolver<algebra_type>::set_convergence_check,
						"", "Check||invokeOnChange=true");

	// 	LUSolver
		reg.add_class_<	LUSolver<algebra_type>,
						ILinearOperatorInverse<vector_type, vector_type> >("LU", grp3.c_str())
			.add_constructor()
			.add_method("set_convergence_check|interactive=false", &LUSolver<algebra_type>::set_convergence_check,
						"", "Check||invokeOnChange=true");

	// 	DirichletDirichletSolver
#ifdef UG_PARALLEL
		{
			typedef DirichletDirichletSolver<algebra_type> T;
			typedef ILinearOperatorInverse<vector_type, vector_type> BaseT;
			reg.add_class_<	T, BaseT >("DirichletDirichlet", grp3.c_str())
			.add_constructor()
			.add_method("set_convergence_check|interactive=false", &T::set_convergence_check,
						"", "Check||invokeOnChange=true")
			.add_method("set_theta|interactive=false", &T::set_theta,
						"", "Theta||invokeOnChange=true")
			.add_method("set_neumann_solver|interactive=false", &T::set_neumann_solver,
						"", "Neumann Solver||invokeOnChange=true")
			.add_method("set_dirichlet_solver|interactive=false", &T::set_dirichlet_solver,
						"", "Dirichlet Solver||invokeOnChange=true")
			.add_method("set_debug", &T::set_debug);
		}
#endif
	// 	LocalSchurComplement
#ifdef UG_PARALLEL
		{
			typedef LocalSchurComplement<algebra_type> T;
			reg.add_class_<	T, ILinearOperator<vector_type, vector_type> >("LocalSchurComplement", grp3.c_str())
			.add_constructor()
			.add_method("set_matrix|interactive=false", &T::set_matrix,
						"", "Matrix||invokeOnChange=true")
			.add_method("set_dirichlet_solver|interactive=false", &T::set_dirichlet_solver,
						"", "Dirichlet Solver||invokeOnChange=true")
			.add_method("set_debug", &T::set_debug)
			// the following functions would normally not be executed from script
			.add_method("init", (bool (T::*)())&T::init)
			.add_method("apply", &T::apply,
						"Success", "local SC times Vector", "Vector");
		}
#endif
	// 	FETISolver
#ifdef UG_PARALLEL
		{
			typedef FETISolver<algebra_type> T;
			typedef IMatrixOperatorInverse<vector_type, vector_type, matrix_type> BaseT;
			reg.add_class_<	T, BaseT >("FETI", grp3.c_str())
			.add_constructor()
			.add_method("set_convergence_check|interactive=false", &T::set_convergence_check,
						"", "Check||invokeOnChange=true")
			.add_method("set_neumann_solver|interactive=false", &T::set_neumann_solver,
						"", "Neumann Solver||invokeOnChange=true")
			.add_method("set_dirichlet_solver|interactive=false", &T::set_dirichlet_solver,
						"", "Dirichlet Solver||invokeOnChange=true")
			.add_method("set_domain_decomp_info", &T::set_domain_decomp_info)
			.add_method("set_debug", &T::set_debug);
		}
#endif

	}

}



bool RegisterStaticLibAlgebraInterface(Registry& reg, const char* parentGroup)
{
	try
	{
	//	get group string
		std::stringstream groupString; groupString << parentGroup << "/Algebra";
		std::string grp = groupString.str();

		// Chooser Interface
		reg.add_class_<	AlgebraTypeChooserInterface >("AlgebraTypeChooserInterface", grp.c_str());
		reg.add_class_<	CPUAlgebraChooser,AlgebraTypeChooserInterface >("CPUAlgebraChooser", grp.c_str())
			.add_constructor()
			.add_method("set_fixed_blocksize", &CPUAlgebraChooser::set_fixed_blocksize, "", "blocksize")
			.add_method("set_variable_blocksize", &CPUAlgebraChooser::set_variable_blocksize);

		// StandardConvCheck
		reg.add_class_<IConvergenceCheck>("IConvergenceCheck", grp.c_str());

		reg.add_class_<StandardConvCheck, IConvergenceCheck>("StandardConvergenceCheck", grp.c_str())
			.add_constructor()
			.add_method("set_maximum_steps|interactive=false", &StandardConvCheck::set_maximum_steps,
					"", "Maximum Steps||invokeOnChange=true")
			.add_method("set_minimum_defect|interactive=false", &StandardConvCheck::set_minimum_defect,
					"", "Minimum Defect||invokeOnChange=true")
			.add_method("set_reduction|interactive=false", &StandardConvCheck::set_reduction,
					"", "Reduction||invokeOnChange=true")
			.add_method("set_verbose_level|interactive=false", &StandardConvCheck::set_verbose_level,
					"", "Verbose||invokeOnChange=true");

	}
	catch(UG_REGISTRY_ERROR_RegistrationFailed ex)
	{
		UG_LOG("### ERROR in RegisterLibAlgebraInterface: "
				"Registration failed (using name " << ex.name << ").\n");
		return false;
	}

	return true;
}

bool RegisterDynamicLibAlgebraInterface(Registry& reg, int algebra_type, const char* parentGroup)
{
	try
	{
		//	get group string
		std::stringstream groupString; groupString << parentGroup << "/Algebra";
		std::string grp = groupString.str();

		// register algebra
		switch(algebra_type)
		{
		case eCPUAlgebra:		 		RegisterAlgebraType<CPUAlgebra >(reg, grp.c_str()); break;
//		case eCPUBlockAlgebra2x2: 		RegisterAlgebraType<CPUBlockAlgebra<2> >(reg, grp.c_str()); break;
		case eCPUBlockAlgebra3x3: 		RegisterAlgebraType<CPUBlockAlgebra<3> >(reg, grp.c_str()); break;
//		case eCPUBlockAlgebra4x4: 		RegisterAlgebraType<CPUBlockAlgebra<4> >(reg, grp.c_str()); break;
//		case eCPUVariableBlockAlgebra: 	RegisterAlgebraType<CPUVariableBlockAlgebra>(reg, grp.c_str()); break;
		default: UG_ASSERT(0, "In RegisterDynamicLibAlgebraInterface: " << algebra_type << " is unsupported algebra type");
					UG_LOG("In RegisterDynamicLibAlgebraInterface: " << algebra_type << " is unsupported algebra type");
					return false;
		}
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
