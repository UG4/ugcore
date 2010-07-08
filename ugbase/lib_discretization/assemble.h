#ifndef __H__LIB_DISCRETIZATION__ASSEMBLE__
#define __H__LIB_DISCRETIZATION__ASSEMBLE__


// other ug4 modules
#include "lib_algebra/lib_algebra.h"

// library intern includes
#include "lib_discretization/domain.h"
#include "lib_discretization/function_spaces/grid_function_space.h"

namespace ug{

enum IAssembleReturn {

	IAssemble_OK = 0,
	IAssemble_NOT_IMPLEMENTED = 1,
	IAssemble_ERROR = 2,
	IAssemble_NON_LINEAR = 3,
	IAssemble_TIME_INDEPENDENT = 4

};


////////////////////////////////////////////////////////////////////////////////////////////////
//	IAssemble
///	Interface providing Jacobian and Defect (in case of linear problems also right-hand side) of a discretization.
/**
 * Interface to generate Jacobi - Matrices and Defect - vector_types for a nonlinear problem
 * and to compute matrix_type and right-hand side for a linear problem
 *
 * The Interface can be used directly to compute Jacobian and Defect (resp. Mass-Stiffness-matrix_type and right-hand side)
 * in case of nonlinear (resp. linear) problem.
 * Furthermore, the interface is used in the NewtonSolver and MultiGridSolver.
 * The NewtonSolver uses the functions assemble_defect+assemble_jacobian or assemble_jacobian_defect.
 * The MultiGridSolver uses the function assemble_jacobian in order to assemble coarse grid Matrices.
 *
 * To give an idea of the usage, an implementation should obey this rules:
 *
 * - time dependent nonlinear problem \$ \partial_t u + A(u) = f \$. Using a l-step time stepping, the defect will be
 *   \[
 *   	d(u^k) = \sum_{i=0}^{l-1} s_{m,i} M(u^{k-i}) + s_{a,i} \{A(u^{k-i}) - f\}
 *   \]
 *   and assemble_linear will return false.
 *
 * - time dependent linear problem \$ \partial_t u + A*u = f \$. Using a l-step time stepping, the defect will be
 * 	 \[
 * 		d(u^k) = \sum_{i=0}^{l-1} s_{m,i} M*u^{k-i} + s_{a,i} \{A*u^{k-i} - f \}
 * 		J(u^k) = s_{m,0}*M + s_{a,0} A
 *   \]
 *   and assemble_linear will compute A = s_{m,0}*M + s_{a,0} A, b = - \{ \sum_{i=1}^{l-1} s_{m,i} M*u^{k-i} + s_{a,i} \{A*u^{k-i} - f\} \} + s_{a,0} * f
 *
 * - static nonlinear Problem \$ A(u) = f \$. Then, the defect will be
 *   \[
 *   	d(u) = A(u) - f
 *   \]
 *   and assemble_linear will return false.
 *
 * - static linear Problem \$ A*u = f \$. Then, the defect will be
 *   \[
 *   	d(u) = A*u - f
 *   	J(u) = A
 *   \]
 *   and assemble_linear will compute A = A and b = f.
 *
 *
 *
 */
template <	typename TDiscreteFunction,
			typename TAlgebra = typename TDiscreteFunction::algebra_type >
class IAssemble {
	public:
		// forward types and constants

		// algebra type
		typedef TAlgebra algebra_type;

		// type of algebra matrix
		typedef typename TAlgebra::matrix_type matrix_type;

		// type of algebra vector
		typedef typename TAlgebra::vector_type vector_type;

		// type of discrete function
		typedef TDiscreteFunction discrete_function_type;

	public:
		/// assembles Jacobian (or Approximation of Jacobian)
		/**
		 * Assembles Jacobian at a given Solution u.
		 * The size of matrix_type has to match the DoFPattern of the NumericalSolution
		 *
		 * \param[out] J Jacobian J(u) (or Precondition) matrix_type to be filled
		 * \param[in]  u Numerical solution
		 */
		virtual IAssembleReturn assemble_jacobian(matrix_type& J, const discrete_function_type& u)
		{return IAssemble_NOT_IMPLEMENTED;}

		/// assembles Defect
		/**
		 * Assembles Defect at a given Solution u.
		 * The size of vector_type has to match the DoFPattern of the NumericalSolution
		 *
		 * \param[out] d Defect d(u) to be filled
		 * \param[in]  u Numerical solution
		 */
		virtual IAssembleReturn assemble_defect(vector_type& d, const discrete_function_type& u)
		{return IAssemble_NOT_IMPLEMENTED;}

		/// Assembles matrix_type and Right-Hand-Side for a linear problem
		/**
		 * Assembles matrix_type and Right-Hand-Side for a linear problem
		 * The size of matrix_type and vector_type have to match the DoFPattern of the Numerical Solution
		 *
		 * \param[out] A Mass-/Stiffness- matrix_type of the discretization
		 * \param[out] b Right-Hand-Side of the discretization
		 *
		 * \return 	IAssemble_OK 			if problem is linear and assembling successful
		 * 			IAssemble_ERROR 		if problem is linear and an error occurred during assembling
		 * 			IAssemble_NONLINEAR 	if problem is non-linear
		 */
		virtual IAssembleReturn assemble_linear(matrix_type& A, vector_type& b, const discrete_function_type& u)
		{return IAssemble_NOT_IMPLEMENTED;}

		/// sets dirichlet values in solution vector
		/**
		 * Sets dirichlet values of the NumericalSolution u when components are dirichlet
		 *
		 * \param[out] u Numerical Solution where dirichlet values are to be set.
		 *
		 * \return 	IAssemble_OK 				if function is implemented and assembling successful
		 * 			IAssemble_NOT_IMPLEMENTED 	if function has not been implemented
		 * 			IAssemble_ERROR 			if function is implemented and an error occurred during assembling
		 */
		virtual IAssembleReturn assemble_solution(discrete_function_type& u)
		{return IAssemble_NOT_IMPLEMENTED;}

		virtual size_t num_fct() const = 0;

		virtual bool is_dirichlet(int si, size_t fct) = 0;

		/// virtual Destructor
		virtual ~IAssemble(){};
};

}; // name space ug

#endif /* __H__LIB_DISCRETIZATION__ASSEMBLE__ */
