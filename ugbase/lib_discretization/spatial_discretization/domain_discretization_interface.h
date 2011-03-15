/*
 * domain_discretization_interface.h
 *
 *  Created on: 01.03.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DOMAIN_DISCRETIZATION_INTERFACE__
#define __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DOMAIN_DISCRETIZATION_INTERFACE__

#include "lib_discretization/assemble_interface.h"
#include "./post_process/post_process_interface.h"

namespace ug {

/// \ingroup lib_disc_domain_assemble
/// @{

/// Interface for domain disretization
/**
 * This class is the interface for spatial discretizations. It can be used
 * in the stationary case as well as for the domain dependent part of an
 * instationary problem (i.e. inside ITimeDiscretization).
 *
 * By its structure it is convenient to implement element-wise Mass-Matrix and
 * Stiffness-Matrix. Then the time-independent member functions can call only
 * the Stiffness-Matrix assembling, while the time-dependent part can call Mass-
 * and Stiffness-Matrix assembling.
 * \tparam		TDoFDistribution 	DoF Distribution Type
 * \tparam		TAlgebra			Algebra Type
 */
template <	typename TDoFDistribution,
			typename TAlgebra>
class IDomainDiscretization : public IAssemble<TDoFDistribution, TAlgebra>{
	public:
	// 	DoF Distribution Type
		typedef IDoFDistribution<TDoFDistribution> dof_distribution_type;

	// 	Algebra type
		typedef TAlgebra algebra_type;

	// 	Type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

	// 	Type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

	public:
	/// assembles Jacobian (or Approximation of Jacobian)
	/**
	 * Assembles Jacobian at a given Solution u.
	 *
	 * \param[out] J 		Jacobian J(u) Matrix to be filled
	 * \param[in]  u 		Current solution
	 * \param[in]  dofDistr DoF Distribution
	 * \param[in]  time		current time
	 * \param[in]  s_m		scaling for mass matrix
	 * \param[in]  s_a		scaling for stiffness matrix
	 *
	 * \return 	IAssemble_OK  				if time dependent and successful
	 * 			IAssemble_ERROR 			if time dependent and error occurred
	 * 			IAssemble_TIMEINDEPENDENT 	if problem is time independent
	 */
		virtual
		IAssembleReturn assemble_jacobian(matrix_type& J,
		                                  const vector_type& u,
		                                  const dof_distribution_type& dofDistr,
		                                  number time, number s_m, number s_a)
		{return IAssemble_NOT_IMPLEMENTED;}

	/// assembles Defect
	/**
	 * Assembles Defect at a given Solution u.
	 *
	 * \param[out] d Defect d(u) to be filled
	 * \param[in]  u 		Current solution
	 * \param[in]  dofDistr DoF Distribution
	 * \param[in]  time		current time
	 * \param[in]  s_m		scaling for mass matrix
	 * \param[in]  s_a		scaling for stiffness matrix
	 *
	 * \return 	IAssemble_OK  				if time dependent and successful
	 * 			IAssemble_ERROR 			if time dependent and error occurred
	 * 			IAssemble_TIMEINDEPENDENT 	if problem is time independent
	 */
		virtual
		IAssembleReturn assemble_defect(vector_type& d,
		                                const vector_type& u,
		                                const dof_distribution_type& dofDistr,
		                                number time, number s_m, number s_a)
		{return IAssemble_NOT_IMPLEMENTED;}

	/// Assembles matrix_type and Right-Hand-Side for a linear problem
	/**
	 * Assembles matrix_type and Right-Hand-Side for a linear problem
	 *
	 * \param[out] A Mass-/Stiffness- matrix_type of the discretization
	 * \param[out] b Right-Hand-Side of the discretization
	 * \param[in]  u 		Current solution
	 * \param[in]  dofDistr DoF Distribution
	 * \param[in]  time		current time
	 * \param[in]  s_m		scaling for mass matrix
	 * \param[in]  s_a		scaling for stiffness matrix
	 *
	 * \return 	IAssemble_OK  				if time dependent and linear and successful
	 * 			IAssemble_ERROR 			if time dependent and linear and error occurred
	 * 			IAssemble_TIMEINDEPENDENT 	if problem is time independent and linear
	 * 			IAssemble_NONLINEAR			if problem is time dependent, but nonlinear
	 */
		virtual
		IAssembleReturn assemble_linear(matrix_type& A,
		                                vector_type& b,
		                                const vector_type& u,
		                                const dof_distribution_type& dofDistr,
		                                number time, number s_m, number s_a)
		{return IAssemble_NOT_IMPLEMENTED;}

	/// sets dirichlet values in solution vector
	/**
	 * Sets dirichlet values of the Solution u when components are dirichlet
	 *
	 * \param[in]  u 		Solution to set values at
	 * \param[in]  dofDistr DoF Distribution
	 * \param[in]  time		current time
	 *
	 * \return 	IAssemble_OK 				if successful
	 * 			IAssemble_NOT_IMPLEMENTED 	if function has not been implemented
	 * 			IAssemble_ERROR 			if implemented but error occurred
	 */
		virtual
		IAssembleReturn assemble_solution(vector_type& u,
		                                  const dof_distribution_type& dofDistr,
		                                  number time)
		{return IAssemble_NOT_IMPLEMENTED;}

	///	returns the number of post processes
		virtual size_t num_post_process() const = 0;

	///	returns the i'th post process
		virtual IPostProcess<TDoFDistribution, TAlgebra>* get_post_process(size_t i) = 0;
};

/// @}

}; // namespace ug

#endif /* __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DOMAIN_DISCRETIZATION_INTERFACE__ */
