



#ifndef __H__LIB_DISCRETIZATION__ASSEMBLE__
#define __H__LIB_DISCRETIZATION__ASSEMBLE__

#include "lib_algebra/lib_algebra.h"
#include "numericalsolution.h"


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
 * Interface to generate Jacobi - Matrices and Defect - Vectors for a nonlinear problem
 * and to compute Matrix and right-hand side for a linear problem
 *
 * The Interface can be used directly to compute Jacobian and Defect (resp. Mass-Stiffness-Matrix and right-hand side)
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
template <int dim>
class IAssemble {

	public:
		/// assembles Jacobian (or Approximation of Jacobian) and Defect at a given Solution u
		/**
		 * Assembles Jacobian and Defect at a given Solution u on level 'level'.
		 * The size of Matrix and Vector have to match the DoFPattern of the NumericalSolution
		 *
		 * \param[out] J Jacobian J(u) (or Precondition) Matrix to be filled
		 * \param[out] d Defect d(u) to be filled
		 * \param[in]  u Numerical solution
		 */
		virtual IAssembleReturn assemble_jacobian_defect(Matrix& J, Vector& d, NumericalSolution<dim>& u, uint level = 0)
		{return IAssemble_NOT_IMPLEMENTED;}

		/// assembles Jacobian (or Approximation of Jacobian)
		/**
		 * Assembles Jacobian at a given Solution u.
		 * The size of Matrix has to match the DoFPattern of the NumericalSolution
		 *
		 * \param[out] J Jacobian J(u) (or Precondition) Matrix to be filled
		 * \param[in]  u Numerical solution
		 */
		virtual IAssembleReturn assemble_jacobian(Matrix& J, NumericalSolution<dim>& u, uint level = 0)
		{return IAssemble_NOT_IMPLEMENTED;}

		/// assembles Defect
		/**
		 * Assembles Defect at a given Solution u.
		 * The size of Vector has to match the DoFPattern of the NumericalSolution
		 *
		 * \param[out] d Defect d(u) to be filled
		 * \param[in]  u Numerical solution
		 */
		virtual IAssembleReturn assemble_defect(Vector& d, NumericalSolution<dim>& u, uint level = 0)
		{return IAssemble_NOT_IMPLEMENTED;}

		/// Assembles Matrix and Right-Hand-Side for a linear problem
		/**
		 * Assembles Matrix and Right-Hand-Side for a linear problem
		 * The size of Matrix and Vector have to match the DoFPattern of the Numerical Solution
		 *
		 * \param[out] A Mass-/Stiffness- Matrix of the discretization
		 * \param[out] b Right-Hand-Side of the discretization
		 *
		 * \return 	IAssemble_OK 			if problem is linear and assembling successful
		 * 			IAssemble_ERROR 		if problem is linear and an error occurred during assembling
		 * 			IAssemble_NONLINEAR 	if problem is non-linear
		 */
		virtual IAssembleReturn assemble_linear(Matrix& A, Vector& b, uint level = 0)
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
		virtual IAssembleReturn assemble_solution(NumericalSolution<dim>& u, uint level = 0)
		{return IAssemble_NOT_IMPLEMENTED;}

		/// virtual Destructor
		virtual ~IAssemble(){};
};


/**
 *
 * This class can be used in a time solver, if the member functions are specified.
 * Furthermore, if the member functions inherited from IAssemble<dim> are specified, then
 * it can be used in a Newton Solver and MultiGridSolver
 *
 * By its structure it is convenient to implement elementwise Mass-Matrix and Stiffness-Matrix. Then
 * the time - independent member functions can call only the Stiffness-Matrix assembling, while the
 * time - dependent part can call Mass- and Stiffness - Matrix assembling.
 *
 */
template <int dim>
class ISpacialDiscretization : public IAssemble<dim>{

	public:
		/// assembles Jacobian (or Approximation of Jacobian) and Defect at a given Solution u for a time dependent problem
		/**
		 * Assembles Jacobian and Defect at a given Solution u on level 'level'.
		 * The size of Matrix and Vector have to match the DoFPattern of the NumericalSolution
		 *
		 * \param[out] J Jacobian J(u) (or Precondition) Matrix to be filled
		 * \param[out] d Defect d(u) to be filled
		 * \param[in]  u Numerical solution
		 *
		 * \return 	IAssemble_OK  				if problem is time dependent and assembling successful
		 * 			IAssemble_ERROR 			if problem is time dependent and an error occurred
		 * 			IAssemble_TIMEINDEPENDENT 	if problem is time independent
		 *
		 */
		virtual IAssembleReturn assemble_jacobian_defect(Matrix& J, Vector& d, NumericalSolution<dim>& u, number time, number s_m, number s_a, uint level = 0)
		{return IAssemble_NOT_IMPLEMENTED;}

		/// assembles Jacobian (or Approximation of Jacobian)
		/**
		 * Assembles Jacobian at a given Solution u.
		 * The size of Matrix has to match the DoFPattern of the NumericalSolution
		 *
		 * \param[out] J Jacobian J(u) (or Precondition) Matrix to be filled
		 * \param[in]  u Numerical solution
		 *
		 * \return 	IAssemble_OK  				if problem is time dependent and assembling successful
		 * 			IAssemble_ERROR 			if problem is time dependent and an error occurred
		 * 			IAssemble_TIMEINDEPENDENT 	if problem is time independent
		 */
		virtual IAssembleReturn assemble_jacobian(Matrix& J, NumericalSolution<dim>& u, number time, number s_m, number s_a, uint level = 0)
		{return IAssemble_NOT_IMPLEMENTED;}

		/// assembles Defect
		/**
		 * Assembles Defect at a given Solution u.
		 * The size of Vector has to match the DoFPattern of the NumericalSolution
		 *
		 * \param[out] d Defect d(u) to be filled
		 * \param[in]  u Numerical solution
		 *
		 * \return 	IAssemble_OK  				if problem is time dependent and assembling successful
		 * 			IAssemble_ERROR 			if problem is time dependent and an error occurred
		 * 			IAssemble_TIMEINDEPENDENT 	if problem is time independent
		 */
		virtual IAssembleReturn assemble_defect(Vector& d, NumericalSolution<dim>& u, number time, number s_m, number s_a, uint level = 0)
		{return IAssemble_NOT_IMPLEMENTED;}

		/// Assembles Matrix and Right-Hand-Side for a linear problem
		/**
		 * Assembles Matrix and Right-Hand-Side for a linear problem
		 * The size of Matrix and Vector have to match the DoFPattern of the Numerical Solution
		 *
		 * \param[out] A Mass-/Stiffness- Matrix of the discretization
		 * \param[out] b Right-Hand-Side of the discretization
		 *
		 * \return 	IAssemble_OK  				if problem is time dependent and linear and assembling successful
		 * 			IAssemble_ERROR 			if problem is time dependent and linear an error occurred
		 * 			IAssemble_TIMEINDEPENDENT 	if problem is time independent and linear
		 * 			IAssemble_NONLINEAR			if problem is time dependent, but nonlinear
		 */
		virtual IAssembleReturn assemble_linear(Matrix& A, Vector& b, number time, number s_m, number s_a, uint level = 0)
		{return IAssemble_NOT_IMPLEMENTED;}

		/// sets dirichlet values in solution vector
		/**
		 * Sets dirichlet values of the NumericalSolution u when components are dirichlet
		 *
		 * \param[out] u 		Numerical Solution where dirichlet values are to be set.
		 * \param[in]  time 	current time value
		 *
		 * \return 	IAssemble_OK 				if function is implemented and assembling successful
		 * 			IAssemble_NOT_IMPLEMENTED 	if function has not been implemented
		 * 			IAssemble_ERROR 			if function is implemented and an error occurred during assembling
		 */
		virtual IAssembleReturn assemble_solution(NumericalSolution<dim>& u, number time, uint level = 0)
		{return IAssemble_NOT_IMPLEMENTED;}

};


};


#endif /* __H__LIB_DISCRETIZATION__ASSEMBLE__ */
