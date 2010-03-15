/*
 * subset_discretization_interface.h
 *
 *  Created on: 01.03.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__DOMAIN_DISCRETIZATION__SUBSET_DISCRETIZATION__SUBSET_DISCRETIZATION_INTERFACE__
#define __H__LIB_DISCRETIZATION__DOMAIN_DISCRETIZATION__SUBSET_DISCRETIZATION__SUBSET_DISCRETIZATION_INTERFACE__

namespace ug{


template <typename TDomain, typename TAlgebra>
class ISubsetDiscretization : public IAssemble<TDomain, TAlgebra>{
	public:
		// forward types and constants

		// domain type
		typedef TDomain domain_type;

		// algebra type
		typedef TAlgebra algebra_type;

		// type of algebra matrix
		typedef typename TAlgebra::matrix_type matrix_type;

		// type of algebra vector
		typedef typename TAlgebra::vector_type vector_type;

		// type of discrete function
		typedef DiscreteGridFunction<TDomain, TAlgebra> discrete_function_type;

	public:
		/// assembles Jacobian (or Approximation of Jacobian) and Defect at a given Solution u for a time dependent problem
		/**
		 * Assembles Jacobian and Defect at a given Solution u on level 'level'.
		 * The size of matrix_type and vector_type have to match the DoFPattern of the NumericalSolution
		 *
		 * \param[out] J Jacobian J(u) (or Precondition) matrix_type to be filled
		 * \param[out] d Defect d(u) to be filled
		 * \param[in]  u Numerical solution
		 *
		 * \return 	IAssemble_OK  				if problem is time dependent and assembling successful
		 * 			IAssemble_ERROR 			if problem is time dependent and an error occurred
		 * 			IAssemble_TIMEINDEPENDENT 	if problem is time independent
		 *
		 */
		virtual IAssembleReturn assemble_jacobian_defect(matrix_type& J, vector_type& d, discrete_function_type& u, number time, number s_m, number s_a, uint level = 0)
		{return IAssemble_NOT_IMPLEMENTED;}

		/// assembles Jacobian (or Approximation of Jacobian)
		/**
		 * Assembles Jacobian at a given Solution u.
		 * The size of matrix_type has to match the DoFPattern of the NumericalSolution
		 *
		 * \param[out] J Jacobian J(u) (or Precondition) matrix_type to be filled
		 * \param[in]  u Numerical solution
		 *
		 * \return 	IAssemble_OK  				if problem is time dependent and assembling successful
		 * 			IAssemble_ERROR 			if problem is time dependent and an error occurred
		 * 			IAssemble_TIMEINDEPENDENT 	if problem is time independent
		 */
		virtual IAssembleReturn assemble_jacobian(matrix_type& J, discrete_function_type& u, number time, number s_m, number s_a, uint level = 0)
		{return IAssemble_NOT_IMPLEMENTED;}

		/// assembles Defect
		/**
		 * Assembles Defect at a given Solution u.
		 * The size of vector_type has to match the DoFPattern of the NumericalSolution
		 *
		 * \param[out] d Defect d(u) to be filled
		 * \param[in]  u Numerical solution
		 *
		 * \return 	IAssemble_OK  				if problem is time dependent and assembling successful
		 * 			IAssemble_ERROR 			if problem is time dependent and an error occurred
		 * 			IAssemble_TIMEINDEPENDENT 	if problem is time independent
		 */
		virtual IAssembleReturn assemble_defect(vector_type& d, discrete_function_type& u, number time, number s_m, number s_a, uint level = 0)
		{return IAssemble_NOT_IMPLEMENTED;}

		/// Assembles matrix_type and Right-Hand-Side for a linear problem
		/**
		 * Assembles matrix_type and Right-Hand-Side for a linear problem
		 * The size of matrix_type and vector_type have to match the DoFPattern of the Numerical Solution
		 *
		 * \param[out] A Mass-/Stiffness- matrix_type of the discretization
		 * \param[out] b Right-Hand-Side of the discretization
		 *
		 * \return 	IAssemble_OK  				if problem is time dependent and linear and assembling successful
		 * 			IAssemble_ERROR 			if problem is time dependent and linear an error occurred
		 * 			IAssemble_TIMEINDEPENDENT 	if problem is time independent and linear
		 * 			IAssemble_NONLINEAR			if problem is time dependent, but nonlinear
		 */
		virtual IAssembleReturn assemble_linear(matrix_type& A, vector_type& b, number time, number s_m, number s_a, uint level = 0)
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
		virtual IAssembleReturn assemble_solution(discrete_function_type& u, number time, uint level = 0)
		{return IAssemble_NOT_IMPLEMENTED;}

};



} // namespace ug


#endif /* __H__LIB_DISCRETIZATION__DOMAIN_DISCRETIZATION__SUBSET_DISCRETIZATION__SUBSET_DISCRETIZATION_INTERFACE__ */
