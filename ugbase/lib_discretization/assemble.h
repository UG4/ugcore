/*
 * assemble.h
 *
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISCRETIZATION__ASSEMBLE__
#define __H__UG__LIB_DISCRETIZATION__ASSEMBLE__

#include "lib_discretization/dof_manager/dof_distribution.h"

namespace ug{

/// Return types of assembling
enum IAssembleReturn {
	IAssemble_OK = 0,
	IAssemble_NOT_IMPLEMENTED = 1,
	IAssemble_ERROR = 2,
	IAssemble_NON_LINEAR = 3,
	IAssemble_TIME_INDEPENDENT = 4

};

/**
 * \brief Assemblings.
 *
 *	Interfaces to assemble problems.
 *
 * \defgroup lib_disc_assemble Assembling
 * \ingroup lib_discretization
 */

/// \ingroup lib_disc_assemble
/// @{

//////////////////////////////////////////////////////////////////////////
//	IAssemble
///	Interface providing Jacobian and Defect of a discretization.
/**
 * Interface to generate Jacobi-Matrices and Defect-Vectors for a nonlinear
 * problem and to compute Matrix and Right-Hand Side for a linear problem.
 *
 * The Interface can be used directly to compute Jacobian and Defect
 * (resp. Mass-Stiffness-matrix_type and right-hand side) in case of nonlinear
 * (resp. linear) problem. Furthermore, the interface is used in the
 * NewtonSolver and MultiGridSolver. The NewtonSolver uses the functions
 * assemble_defect+assemble_jacobian or assemble_jacobian_defect.
 * The MultiGridSolver uses the function assemble_jacobian in order to
 * assemble coarse grid Matrices.
 *
 * To give an idea of the usage, an implementation should obey this rules:
 *
 * - Time dependent nonlinear problem \f$ \partial_t u + A(u) = f \f$.
 * 	 Using a l-step time stepping, the defect will be
 *   \f[
 *   	d(u^k) = \sum_{i=0}^{l-1} s_{m,i} M(u^{k-i})
 *   							+ s_{a,i} \{A(u^{k-i}) - f\}
 *   \f]
 *   and assemble_linear will return false.
 *
 * - Time dependent linear problem \f$ \partial_t u + A*u = f \f$.
 * 	 Using a l-step time stepping, the defect will be
 * 	 \f[
 * 		d(u^k) = \sum_{i=0}^{l-1} s_{m,i} M*u^{k-i} + s_{a,i} \{A*u^{k-i} - f \}
 * 		J(u^k) = s_{m,0}*M + s_{a,0} A
 *   \f]
 *   and assemble_linear will compute \f$ A = s_{m,0}*M + s_{a,0} A \f$ and
 *   \f$ b = - \{ \sum_{i=1}^{l-1} s_{m,i} M*u^{k-i}
 *   			+ s_{a,i} \{A*u^{k-i} - f\} \} + s_{a,0} * f \f$
 *
 * - Stationary Non-linear Problem \f$ A(u) = f \f$. Then, the defect will be
 *   \f[
 *   	d(u) = A(u) - f
 *   \f]
 *   and assemble_linear will return false.
 *
 * - Stationary linear Problem \f$ A*u = f \f$. Then, the defect will be
 *   \f[
 *   	d(u) = A*u - f
 *   	J(u) = A
 *   \f]
 *   and assemble_linear will compute \f$ A = A \f$ and \f$ b = f \f$.
 *
 *	\tparam		TDoFDistribution		DoF Distribution
 *	\tparam		TAlgebra				Algebra type
 */
template <	typename TDoFDistribution,
			typename TAlgebra>
class IAssemble {
	public:
	// 	Algebra type
		typedef TAlgebra algebra_type;

	// 	Type of algebra matrix
		typedef typename TAlgebra::matrix_type matrix_type;

	// 	Type of algebra vector
		typedef typename TAlgebra::vector_type vector_type;

	// 	Type of DoF Distribution
		typedef IDoFDistribution<TDoFDistribution> dof_distribution_type;

	public:
		/// assembles Jacobian (or Approximation of Jacobian)
		/**
		 * Assembles Jacobian at a given iterate u.
		 *
		 * \param[out] 	J 			Jacobian J(u) matrix to be filled
		 * \param[in]  	u 			Current iterate
		 * \param[in]	dofDistr	DoF Distribution
		 */
		virtual IAssembleReturn assemble_jacobian(matrix_type& J,
		                                          const vector_type& u,
		                                          const dof_distribution_type& dofDistr)
		{return IAssemble_NOT_IMPLEMENTED;}

		/// assembles Defect
		/**
		 * Assembles Defect at a given Solution u.
		 *
		 * \param[out] 	d 			Defect d(u) to be filled
		 * \param[in] 	u 			Current iterate
		 * \param[in]	dofDistr	DoF Distribution
		 */
		virtual IAssembleReturn assemble_defect(vector_type& d,
		                                        const vector_type& u,
		                                        const dof_distribution_type& dofDistr)
		{return IAssemble_NOT_IMPLEMENTED;}

		/// Assembles Matrix and Right-Hand-Side for a linear problem
		/**
		 * Assembles matrix_type and Right-Hand-Side for a linear problem
		 *
		 * \param[out] 	A 			Mass-/Stiffness- Matrix
		 * \param[out] 	b 			Right-Hand-Side
		 * \param[in] 	u 			Current iterate
		 * \param[in]	dofDistr	DoF Distribution
		 *
		 * \return 	IAssemble_OK 		if problem is linear and assembling successful
		 * 			IAssemble_ERROR 	if problem is linear and an error occurred during assembling
		 * 			IAssemble_NONLINEAR if problem is non-linear
		 */
		virtual IAssembleReturn assemble_linear(matrix_type& A,
		                                        vector_type& b,
		                                        const vector_type& u,
		                                        const dof_distribution_type& dofDistr)
		{return IAssemble_NOT_IMPLEMENTED;}

		/// sets dirichlet values in solution vector
		/**
		 * Sets dirichlet values of the NumericalSolution u when components
		 * are dirichlet
		 *
		 * \param[out] 	u			Numerical Solution
		 * \param[in]	dofDistr	DoF Distribution
		 *
		 * \return 	IAssemble_OK 				if function is implemented and assembling successful
		 * 			IAssemble_NOT_IMPLEMENTED 	if function has not been implemented
		 * 			IAssemble_ERROR 			if function is implemented and an error occurred during assembling
		 */
		virtual IAssembleReturn assemble_solution(vector_type& u,
		                                          const dof_distribution_type& dofDistr)
		{return IAssemble_NOT_IMPLEMENTED;}


	/// forces the assembling to consider the grid as regular
		virtual void force_regular_grid(bool bForce) = 0;

	//	todo: Remove, iff possible
		virtual size_t num_fct() const = 0;

	/// Virtual Destructor
		virtual ~IAssemble(){};
};

/// @}

}; // name space ug

#endif /* __H__UG__LIB_DISCRETIZATION__ASSEMBLE__ */
