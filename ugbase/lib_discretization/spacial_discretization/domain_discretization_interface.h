/*
 * domain_discretization_interface.h
 *
 *  Created on: 01.03.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__DOMAIN_DISCRETIZATION_INTERFACE__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__DOMAIN_DISCRETIZATION_INTERFACE__

#include "lib_discretization/assemble.h"

namespace ug {

/**
 *
 * This class can be used in a time solver, if the member functions are specified.
 * Furthermore, if the member functions inherited from IAssemble<dim> are specified, then
 * it can be used in a Newton Solver and MultiGridSolver
 *
 * By its structure it is convenient to implement elementwise Mass-matrix_type and Stiffness-matrix_type. Then
 * the time - independent member functions can call only the Stiffness-matrix_type assembling, while the
 * time - dependent part can call Mass- and Stiffness - matrix_type assembling.
 *
 */
template <	typename TDiscreteFunction,
			typename TAlgebra = typename TDiscreteFunction::algebra_type >
class IDomainDiscretization : public IAssemble<TDiscreteFunction, TAlgebra>{
	public:
		// forward types and constants

		// discrete function type
		typedef TDiscreteFunction discrete_function_type;

		// algebra type
		typedef TAlgebra algebra_type;

		// type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

		// type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

	public:
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
		virtual IAssembleReturn assemble_jacobian(matrix_type& J, const discrete_function_type& u, number time, number s_m, number s_a)
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
		virtual IAssembleReturn assemble_defect(vector_type& d, const discrete_function_type& u, number time, number s_m, number s_a)
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
		virtual IAssembleReturn assemble_linear(matrix_type& A, vector_type& b, const discrete_function_type& u, number time, number s_m, number s_a)
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
		virtual IAssembleReturn assemble_solution(discrete_function_type& u, number time)
		{return IAssemble_NOT_IMPLEMENTED;}

		/// returns if the number of functions of this assembling
		virtual size_t num_fct() const = 0;

		/// returns true if the subset is dirchlet for the function fct
		virtual bool is_dirichlet(int si, size_t fct) = 0;
};

template <	typename TDiscreteFunction,
			typename TAlgebra = typename TDiscreteFunction::algebra_type >
class IDimensionDomainDiscretization{
	public:
		// forward types and constants

		// discrete function type
		typedef TDiscreteFunction discrete_function_type;

		// algebra type
		typedef TAlgebra algebra_type;

		// type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

		// type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

	public:
		virtual IAssembleReturn assemble_jacobian(matrix_type& J, const discrete_function_type& u, int si)
		{return IAssemble_NOT_IMPLEMENTED;}
		virtual IAssembleReturn assemble_defect(vector_type& d, const discrete_function_type& u, int si)
		{return IAssemble_NOT_IMPLEMENTED;}
		virtual IAssembleReturn assemble_linear(matrix_type& mat, vector_type& rhs, const discrete_function_type& u, int si)
		{return IAssemble_NOT_IMPLEMENTED;}


		virtual IAssembleReturn assemble_jacobian(matrix_type& J, const discrete_function_type& u, int si, number time, number s_m, number s_a)
		{return IAssemble_NOT_IMPLEMENTED;}
		virtual IAssembleReturn assemble_defect(vector_type& d, const discrete_function_type& u, int si, number time, number s_m, number s_a)
		{return IAssemble_NOT_IMPLEMENTED;}
		virtual IAssembleReturn assemble_linear(matrix_type& A, vector_type& b, const discrete_function_type& u, int si, number time, number s_m, number s_a)
		{return IAssemble_NOT_IMPLEMENTED;}

		/// returns if the number of functions of this assembling
		virtual size_t num_fct() const = 0;

		virtual ~IDimensionDomainDiscretization() {};
};

template <	typename TDiscreteFunction,
			typename TAlgebra = typename TDiscreteFunction::algebra_type >
class IDirichletBoundaryValues{
	public:
		// forward types and constants

		// discrete function type
		typedef TDiscreteFunction discrete_function_type;

		// algebra type
		typedef TAlgebra algebra_type;

		// type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

		// type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

	public:
		virtual IAssembleReturn clear_dirichlet_jacobian(matrix_type& J, const discrete_function_type& u, int si, number time = 0.0)
		{return IAssemble_NOT_IMPLEMENTED;}
		virtual IAssembleReturn clear_dirichlet_defect(vector_type& d, const discrete_function_type& u, int si,number time = 0.0)
		{return IAssemble_NOT_IMPLEMENTED;}
		virtual IAssembleReturn set_dirichlet_linear(matrix_type& mat, vector_type& rhs, const discrete_function_type& u, int si, number time = 0.0)
		{return IAssemble_NOT_IMPLEMENTED;}

		virtual IAssembleReturn set_dirichlet_solution(discrete_function_type& u, int si, number time = 0.0)
		{return IAssemble_NOT_IMPLEMENTED;}

		/// returns true if the subset is dirchlet for the function fct
		virtual bool is_dirichlet(size_t fct) = 0;

		virtual ~IDirichletBoundaryValues() {};
};


}; // namespace ug



#endif /* __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__DOMAIN_DISCRETIZATION_INTERFACE__ */
