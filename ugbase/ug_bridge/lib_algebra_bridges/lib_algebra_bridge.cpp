/*
 * lib_algebra_bridge.cpp
 *
 *  Created on: 11.10.2010
 *      Author: andreasvogel
 */

// extern headers
#include <iostream>
#include <sstream>
#include <string>

// bridge
#include "ug_bridge/ug_bridge.h"

// algebra includes
#include "lib_algebra_bridge.h"

// \todo: extract only really needed includes
// all parts of lib algebra
#include "lib_algebra/lib_algebra.h"
#include "lib_algebra/operator/operator_impl.h"

using namespace std;

namespace ug
{
namespace bridge
{

//! calculates dest = alpha1*v1 + alpha2*v2
template<typename vector_t>
inline void VecScaleAdd2(vector_t &dest, double alpha1, const vector_t &v1, double alpha2, const vector_t &v2)
{
	VecScaleAdd(dest, alpha1, v1, alpha2, v2);
}

//! calculates dest = alpha1*v1 + alpha2*v2 + alpha3*v3
template<typename vector_t>
inline void VecScaleAdd3(vector_t &dest, double alpha1, const vector_t &v1, double alpha2, const vector_t &v2, double alpha3, const vector_t &v3)
{
	VecScaleAdd(dest, alpha1, v1, alpha2, v2, alpha3, v3);
}

template <typename TAlgebra>
struct cRegisterAlgebraType
{
static bool reg(Registry& reg, string parentGroup)
{
//	get group string (use same as parent)
	string grp = string(parentGroup);

//	typedefs for this algebra
	typedef typename TAlgebra::vector_type vector_type;
	typedef typename TAlgebra::matrix_type matrix_type;

//	suffix and tag
	string algSuffix = GetAlgebraSuffix<TAlgebra>();
	string algTag = GetAlgebraTag<TAlgebra>();

//	Vector
	{
		string name = string("Vector").append(algSuffix);
		reg.add_class_<vector_type>(name, grp)
		.add_constructor()
		.add_method("set|hide=true", static_cast<bool (vector_type::*)(number)>(&vector_type::set),
								"Success", "Number")
		.add_method("size|hide=true", static_cast<size_t (vector_type::*)() const>(&vector_type::size),
								"Size", "")
		.add_method("set_random|hide=true", static_cast<bool (vector_type::*)(number, number)>(&vector_type::set_random),
								"Success", "Number")
		.add_method("print|hide=true", &vector_type::p);
		reg.add_class_to_group(name, "Vector", algTag);

		reg.add_function("VecScaleAssign", static_cast<void (*)(vector_type&, number, const vector_type &)>(&VecScaleAssign<vector_type>));
		reg.add_function("VecAssign", static_cast<void (*)(vector_type&,const vector_type &)>(&VecAssign<vector_type>));
		reg.add_function("VecScaleAdd2", /*(void (*)(vector_type&, number, const vector_type&, number, const vector_type &)) */
				&VecScaleAdd2<vector_type>, "", "alpha1*vec1 + alpha2*vec2",
				"dest#alpha1#vec1#alpha2#vec2");
		reg.add_function("VecScaleAdd3", /*(void (*)(vector_type&, number, const vector_type&, number, const vector_type &, number, const vector_type &))*/
				&VecScaleAdd3<vector_type>, "", "alpha1*vec1 + alpha2*vec2 + alpha3*vec3",
				"dest#alpha1#vec1#alpha2#vec2#alpha3#vec3");
	}

//	Matrix
	{
		string name = string("Matrix").append(algSuffix);
		reg.add_class_<matrix_type>(name, grp)
			.add_constructor()
			.add_method("print|hide=true", &matrix_type::p);
		reg.add_class_to_group(name, "Matrix", algTag);
	}

//	ApplyLinearSolver
	{
		reg.add_function( "ApplyLinearSolver",
						  &ApplyLinearSolver<vector_type>, grp);
	}

// Debug Writer (abstract base class)
	{
		typedef IDebugWriter<TAlgebra> T;
		string name = string("IDebugWriter").append(algSuffix);
		reg.add_class_<T>(name, grp);
		reg.add_class_to_group(name, "IDebugWriter", algTag);
	}

// IVectorWriter (abstract base class)
	{
		typedef IVectorWriter<vector_type> T;
		string name = string("IVectorWriter").append(algSuffix);
		reg.add_class_<T>(name, grp)
				.add_method("update", &T::update, "", "v", "updates the vector v");
		reg.add_class_to_group(name, "IVectorWriter", algTag);
	}

/////////////////////////
//	Base Classes
/////////////////////////

//	ILinearOperator
	{
		typedef ILinearOperator<vector_type, vector_type> T;
		string name = string("ILinearOperator").append(algSuffix);
		reg.add_class_<T>(name, grp)
			.add_method("init", static_cast<bool(T::*)()>(&T::init));
		reg.add_class_to_group(name, "ILinearOperator", algTag);
	}

// 	MatrixOperator
	{
		typedef ILinearOperator<vector_type, vector_type> TBase;
		typedef MatrixOperator<vector_type, vector_type, matrix_type> T;
		string name = string("MatrixOperator").append(algSuffix);
		reg.add_class_<T, TBase, matrix_type>(name, grp)
			.add_constructor();
		reg.add_class_to_group(name, "MatrixOperator", algTag);
	}

//	ILinearIterator
	{
		typedef ILinearIterator<vector_type, vector_type> T;
		string name = string("ILinearIterator").append(algSuffix);
		reg.add_class_<T>(name, grp);
		reg.add_class_to_group(name, "ILinearIterator", algTag);
	}

//	IPreconditioner
	{
		typedef IPreconditioner<TAlgebra> T;
		typedef ILinearIterator<vector_type, vector_type>  TBase;
		string name = string("IPreconditioner").append(algSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_method("set_debug", &T::set_debug, "sets a debug writer", "d");
		reg.add_class_to_group(name, "IPreconditioner", algTag);
	}

//	ILinearOperatorInverse
	{
		typedef ILinearOperatorInverse<vector_type, vector_type> T;
		string name = string("ILinearOperatorInverse").append(algSuffix);
		reg.add_class_<T>(name, grp)
			.add_method("init", static_cast<bool(T::*)(ILinearOperator<vector_type,vector_type>&)>(&T::init))
			.add_method("apply_return_defect", &T::apply_return_defect, "Success", "u#f",
					"Solve A*u = f, such that u = A^{-1} f by iterating u := u + B(f - A*u),  f := f - A*u becomes new defect")
			.add_method("apply", &T::apply, "Success", "u#f", "Solve A*u = f, such that u = A^{-1} f by iterating u := u + B(f - A*u), f remains constant");
		reg.add_class_to_group(name, "ILinearOperatorInverse", algTag);
	}

//	IMatrixOperatorInverse
	{
		typedef ILinearOperatorInverse<vector_type, vector_type>  TBase;
		typedef IMatrixOperatorInverse<vector_type, vector_type, matrix_type> T;
		string name = string("IMatrixOperatorInverse").append(algSuffix);
		reg.add_class_<T, TBase>(name, grp);
		reg.add_class_to_group(name, "IMatrixOperatorInverse", algTag);
	}

//	IOperator
	{
		typedef IOperator<vector_type, vector_type> T;
		string name = string("IOperator").append(algSuffix);
		reg.add_class_<T>(name, grp);
		reg.add_class_to_group(name, "IOperator", algTag);
	}

//	IOperatorInverse
	{
		typedef IOperatorInverse<vector_type, vector_type> T;
		string name = string("IOperatorInverse").append(algSuffix);
		reg.add_class_<T>(name, grp);
		reg.add_class_to_group(name, "IOperatorInverse", algTag);
	}

//	IProlongationOperator
	{
		typedef IProlongationOperator<vector_type, vector_type> T;
		string name = string("IProlongationOperator").append(algSuffix);
		reg.add_class_<T>(name, grp);
		reg.add_class_to_group(name, "IProlongationOperator", algTag);
	}

//	IProjectionOperator
	{
		typedef IProjectionOperator<vector_type, vector_type> T;
		string name = string("IProjectionOperator").append(algSuffix);
		reg.add_class_<T>(name, grp);
		reg.add_class_to_group(name, "IProjectionOperator", algTag);
	}

//////////////////////
// Preconditioner
//////////////////////
//	get group string
	stringstream grpSS2; grpSS2 << grp << "/Preconditioner";
	string grp2 = grpSS2.str();

//	Jacobi
	{
		typedef Jacobi<TAlgebra> T;
		typedef IPreconditioner<TAlgebra> TBase;
		string name = string("Jacobi").append(algSuffix);
		reg.add_class_<T,TBase>(name, grp2, "Jacobi Preconditioner")
			.add_constructor()
			.add_method("set_damp", &T::set_damp, "", "damp");
		reg.add_class_to_group(name, "Jacobi", algTag);
	}

//	GaussSeidel
	{
		typedef GaussSeidel<TAlgebra> T;
		typedef IPreconditioner<TAlgebra> TBase;
		string name = string("GaussSeidel").append(algSuffix);
		reg.add_class_<T,TBase>(name, grp2, "Gauss-Seidel Preconditioner")
		.add_constructor();
		reg.add_class_to_group(name, "GaussSeidel", algTag);
	}

//	Symmetric GaussSeidel
	{
		typedef SymmetricGaussSeidel<TAlgebra> T;
		typedef IPreconditioner<TAlgebra> TBase;
		string name = string("SymmetricGaussSeidel").append(algSuffix);
		reg.add_class_<T,TBase>(name, grp2)
			.add_constructor();
		reg.add_class_to_group(name, "SymmetricGaussSeidel", algTag);
	}

//	Backward GaussSeidel
	{
		typedef BackwardGaussSeidel<TAlgebra> T;
		typedef IPreconditioner<TAlgebra> TBase;
		string name = string("BackwardGaussSeidel").append(algSuffix);
		reg.add_class_<T,TBase>(name, grp2)
			.add_constructor();
		reg.add_class_to_group(name, "BackwardGaussSeidel", algTag);
	}

//	ILU
	{
		typedef ILU<TAlgebra> T;
		typedef IPreconditioner<TAlgebra> TBase;
		string name = string("ILU").append(algSuffix);
		reg.add_class_<T,TBase>(name, grp2, "Incomplete LU Decomposition")
			.add_constructor()
			.add_method("set_beta", &T::set_beta, "", "beta");
		reg.add_class_to_group(name, "ILU", algTag);
	}

//	ILU Threshold
	{
		typedef ILUTPreconditioner<TAlgebra> T;
		typedef IPreconditioner<TAlgebra> TBase;
		string name = string("ILUT").append(algSuffix);
		reg.add_class_<T,TBase>(name, grp2, "Incomplete LU Decomposition with threshold")
			.add_constructor()
			.add_method("set_threshold", &T::set_threshold,
						"", "threshold", "sets threshold of incomplete LU factorisation");
		reg.add_class_to_group(name, "ILUT", algTag);
	}

//	LinearIteratorProduct
	{
		typedef LinearIteratorProduct<vector_type, vector_type> T;
		typedef ILinearIterator<vector_type, vector_type> TBase;
		string name = string("LinearIteratorProduct").append(algSuffix);
		reg.add_class_<T,TBase>(name, grp,
						"Linear Iterator consisting of several LinearIterations")
				.add_constructor()
				.add_method("add_iteration", &T::add_iterator, "Add an iterator");
		reg.add_class_to_group(name, "LinearIteratorProduct", algTag);
	}

///////////////////////
//	Solver
///////////////////////

//	get group string
	stringstream grpSS3; grpSS3 << grp << "/Solver";
	string grp3 = grpSS3.str();

// 	LinearSolver
	{
		typedef LinearSolver<TAlgebra> T;
		typedef ILinearOperatorInverse<vector_type, vector_type> TBase;
		string name = string("LinearSolver").append(algSuffix);
		reg.add_class_<T,TBase>(name, grp3)
			.add_constructor()
			.add_method("set_preconditioner|interactive=false", &T::set_preconditioner,
						"", "Preconditioner")
			.add_method("set_convergence_check|interactive=false", &T::set_convergence_check,
						"", "Check")
			.add_method("set_compute_fresh_defect_when_finished", &T::set_compute_fresh_defect_when_finished)
			.add_method("set_debug", &T::set_debug);
		reg.add_class_to_group(name, "LinearSolver", algTag);
	}

// 	CG Solver
	{
		typedef CG<TAlgebra> T;
		typedef ILinearOperatorInverse<vector_type, vector_type> TBase;
		string name = string("CG").append(algSuffix);
		reg.add_class_<T,TBase>(name, grp3, "Conjugate Gradient")
			.add_constructor()
			.add_method("set_preconditioner|interactive=false", &T::set_preconditioner,
						"", "Preconditioner")
			.add_method("set_convergence_check|interactive=false", &T::set_convergence_check,
						"", "Check");
		reg.add_class_to_group(name, "CG", algTag);
	}

// 	BiCGStab Solver
	{
		typedef BiCGStab<TAlgebra> T;
		typedef ILinearOperatorInverse<vector_type, vector_type> TBase;
		string name = string("BiCGStab").append(algSuffix);
		reg.add_class_<T,TBase>(name, grp3)
			.add_constructor()
			.add_method("set_preconditioner|interactive=false", &T::set_preconditioner,
						"", "Preconditioner")
			.add_method("set_convergence_check|interactive=false", &T::set_convergence_check,
						"", "Check");
		reg.add_class_to_group(name, "BiCGStab", algTag);
	}

// 	LU Solver
	{
		typedef LU<TAlgebra> T;
		typedef ILinearOperatorInverse<vector_type, vector_type> TBase;
		string name = string("LU").append(algSuffix);
		reg.add_class_<T,TBase>(name, grp3, "LU-Decomposition exact solver")
			.add_constructor()
			.add_method("set_convergence_check|interactive=false", &T::set_convergence_check,
						"", "Check");
		reg.add_class_to_group(name, "LU", algTag);
	}

	// 	DirichletDirichletSolver
#ifdef UG_PARALLEL
	{
		typedef DirichletDirichletSolver<TAlgebra> T;
		typedef ILinearOperatorInverse<vector_type, vector_type> BaseT;
		string name = string("DirichletDirichlet").append(algSuffix);
		reg.add_class_<	T, BaseT >(name, grp3, "Dirichlet-Dirichlet Domain Decomposition Algorithm")
		.add_constructor()
		.add_method("set_convergence_check|interactive=false", &T::set_convergence_check,
					"", "Check")
		.add_method("set_theta|interactive=false", &T::set_theta,
					"", "Theta", "set damping factor theta")
		.add_method("set_neumann_solver|interactive=false", &T::set_neumann_solver,
					"", "Neumann Solver")
		.add_method("set_dirichlet_solver|interactive=false", &T::set_dirichlet_solver,
					"", "Dirichlet Solver")
		.add_method("set_debug", &T::set_debug);
		reg.add_class_to_group(name, "DirichletDirichlet", algTag);
	}

// 	LocalSchurComplement
	{
		typedef LocalSchurComplement<TAlgebra> T;
		typedef ILinearOperator<vector_type, vector_type> TBase;
		string name = string("LocalSchurComplement").append(algSuffix);
		reg.add_class_<	T, TBase>(name, grp3)
		.add_constructor()
		.add_method("set_matrix|interactive=false", &T::set_matrix,
					"", "Matrix")
		.add_method("set_dirichlet_solver|interactive=false", &T::set_dirichlet_solver,
					"", "Dirichlet Solver")
		.add_method("set_debug", &T::set_debug, "", "d")
		// the following functions would normally not be executed from script
		.add_method("init", (bool (T::*)())&T::init)
		.add_method("apply", &T::apply,
					"Success", "local SC times Vector#Vector");
		reg.add_class_to_group(name, "LocalSchurComplement", algTag);
	}

// 	FETISolver
	{
		typedef FETISolver<TAlgebra> T;
		typedef IMatrixOperatorInverse<vector_type, vector_type, matrix_type> BaseT;
		string name = string("FETI").append(algSuffix);
		reg.add_class_<	T, BaseT >(name, grp3, "FETI Domain Decomposition Solver")
		.add_constructor()
		.add_method("set_convergence_check|interactive=false", &T::set_convergence_check,
					"", "Check")
		.add_method("set_neumann_solver|interactive=false", &T::set_neumann_solver,
					"", "Neumann Solver")
		.add_method("set_dirichlet_solver|interactive=false", &T::set_dirichlet_solver,
					"", "Dirichlet Solver")
		.add_method("set_coarse_problem_solver|interactive=false", &T::set_coarse_problem_solver,
					"", "Coarse Problem Solver")
		.add_method("set_domain_decomp_info", &T::set_domain_decomp_info)
		.add_method("print_statistic_of_inner_solver", &T::print_statistic_of_inner_solver)
		.add_method("set_debug", &T::set_debug)
		.add_method("test_layouts", &T::test_layouts);
		reg.add_class_to_group(name, "FETI", algTag);
	}
#endif

	// 	HLIBSolver
#ifdef USE_HLIBPRO
	{
		typedef HLIBSolver<TAlgebra> T;
		typedef ILinearOperatorInverse<vector_type, vector_type> TBase;
		string name = string("HLIBSolver").append(algSuffix);
		reg.add_class_<	T, TBase>(name, grp3)
		.add_constructor()
		.add_method("set_convergence_check|interactive=false", &T::set_convergence_check,
					"", "Check")
		.add_method("set_hlib_nmin|interactive=false",         &T::set_hlib_nmin,
					"", "HLIB nmin")
		.add_method("set_hlib_accuracy_H|interactive=false",   &T::set_hlib_accuracy_H,
					"", "HLIB accuracy_H")
		.add_method("set_hlib_accuracy_LU|interactive=false",  &T::set_hlib_accuracy_LU,
					"", "HLIB accuracy_LU")
		.add_method("set_hlib_verbosity|interactive=false",    &T::set_hlib_verbosity,
					"", "HLIB verbosity")
		.add_method("set_clustering_method|interactive=false", &T::set_clustering_method,
					"", "Clustering")
		.add_method("set_ps_basename|interactive=false",       &T::set_ps_basename,
					"", "PostScript basename")
		.add_method("check_crs_matrix|interactive=false",      &T::check_crs_matrix,
					"", "Check CRS matrix")
		.add_method("set_debug", &T::set_debug);
		reg.add_class_to_group(name, "HLIBSolver", algTag);
	}
#endif

	return true;
}
};



static bool RegisterLibAlgebra__Common(Registry& reg, string parentGroup)
{
	try
	{
//	get group string
	stringstream groupString; groupString << parentGroup << "/Algebra";
	string grp = groupString.str();

// 	AlgebraSelector Interface
	reg.add_class_<	IAlgebraTypeSelector>("IAlgebraTypeSelector", grp);
	reg.add_class_<	CPUAlgebraSelector, IAlgebraTypeSelector>("CPUAlgebraSelector", grp)
		.add_constructor()
		.add_method("set_fixed_blocksize", &CPUAlgebraSelector::set_fixed_blocksize, "", "blocksize")
		.add_method("set_variable_blocksize", &CPUAlgebraSelector::set_variable_blocksize);

// 	IConvergenceCheck
	reg.add_class_<IConvergenceCheck>("IConvergenceCheck", grp);

// 	StandardConvCheck
	reg.add_class_<StandardConvCheck, IConvergenceCheck>("StandardConvergenceCheck", grp)
		.add_constructor()
		.add_method("set_maximum_steps|interactive=false", &StandardConvCheck::set_maximum_steps,
				"", "Maximum Steps")
		.add_method("set_minimum_defect|interactive=false", &StandardConvCheck::set_minimum_defect,
				"", "Minimum Defect")
		.add_method("set_reduction|interactive=false", &StandardConvCheck::set_reduction,
				"", "Reduction")
		.add_method("set_verbose_level|interactive=false", &StandardConvCheck::set_verbose_level,
				"", "Verbose")
		.add_method("defect|interactive=false", &StandardConvCheck::defect, "defect", "", "returns the current defect")
		.add_method("step|interactive=false", &StandardConvCheck::step, "step", "", "returns the current number of steps")
		.add_method("reduction|interactive=false", &StandardConvCheck::reduction, "reduction", "", "returns the current relative reduction")
		.add_method("iteration_ended|interactive=false", &StandardConvCheck::iteration_ended)
		.add_method("previous_defect|interactive=false", &StandardConvCheck::previous_defect);

// IPositionProvider (abstract base class)
	{
		reg.add_class_<IPositionProvider<1> >("IPositionProvider1d", grp);
		reg.add_class_<IPositionProvider<2> >("IPositionProvider2d", grp);
		reg.add_class_<IPositionProvider<3> >("IPositionProvider3d", grp);
		reg.add_class_to_group("IPositionProvider1d", "IPositionProvider", GetDomainTag<1>());
		reg.add_class_to_group("IPositionProvider2d", "IPositionProvider", GetDomainTag<2>());
		reg.add_class_to_group("IPositionProvider3d", "IPositionProvider", GetDomainTag<3>());
	}

	}catch(UG_REGISTRY_ERROR_RegistrationFailed ex)
	{
		UG_LOG("### ERROR in RegisterLibAlgebra: "
				"Registration failed (using name " << ex.name << ").\n");
		return false;
	}

	return true;
}


bool RegisterAMG(Registry& reg, string parentGroup);
bool RegisterEigensolver(Registry& reg, string parentGroup);

bool RegisterLibAlgebra(Registry& reg, string parentGroup)
{
// 	switch moved to lib_algebra_bridge.h
	RegisterAlgebraClass<cRegisterAlgebraType>(reg, parentGroup);
	RegisterLibAlgebra__Common(reg, parentGroup);

	RegisterAMG(reg, parentGroup);
	RegisterEigensolver(reg, parentGroup);

	return true;
}



} // end namespace bridge
} // end namespace ug
