/*
 * domain_disc_interface.h
 *
 *  Created on: 01.03.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__DOMAIN_DISC_INTERFACE__
#define __H__UG__LIB_DISC__SPATIAL_DISC__DOMAIN_DISC_INTERFACE__

#include "lib_disc/assemble_interface.h"
#include "lib_disc/time_disc/solution_time_series.h"
#include "lib_disc/spatial_disc/constraints/constraint_interface.h"

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
 *
 * \tparam		TAlgebra			Algebra Type
 */
template <typename TAlgebra>
class IDomainDiscretization : public IAssemble<TAlgebra>
{
	public:
	/// Algebra type
		typedef TAlgebra algebra_type;

	/// Type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

	/// Type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

	public:
		/// assembles Jacobian (or Approximation of Jacobian)
		/**
		 * Assembles Jacobian at a given iterate u.
		 *
		 * \param[out] 	J 	Jacobian J(u) matrix to be filled
		 * \param[in]  	u 	Current iterate
		 * \param[in]	dd	DoF Distribution
		 */
		virtual void assemble_jacobian(matrix_type& J, const vector_type& u, GridLevel gl) = 0;

		/// assembles Defect
		/**
		 * Assembles Defect at a given Solution u.
		 *
		 * \param[out] 	d 	Defect d(u) to be filled
		 * \param[in] 	u 	Current iterate
		 * \param[in]	dd	DoF Distribution
		 */
		virtual void assemble_defect(vector_type& d, const vector_type& u, GridLevel gl) = 0;

		/// Assembles Matrix and Right-Hand-Side for a linear problem
		/**
		 * Assembles matrix_type and Right-Hand-Side for a linear problem
		 *
		 * \param[out] 	A 	Mass-/Stiffness- Matrix
		 * \param[out] 	b 	Right-Hand-Side
		 * \param[in]	dd	DoF Distribution
		 */
		virtual void assemble_linear(matrix_type& A, vector_type& b, GridLevel gl) = 0;

		/// sets dirichlet values in solution vector
		/**
		 * Sets dirichlet values of the NumericalSolution u when components
		 * are dirichlet
		 *
		 * \param[out] 	u	Numerical Solution
		 * \param[in]	dd	DoF Distribution
		 */
		virtual void adjust_solution(vector_type& u, GridLevel gl) = 0;


	public:
	/// prepares Timestep
	/**
	 * prepares timestep at a given Solution u.
	 *
	 * \param[in]  vSol			vector of previous and current (iterated) solution
	 * \param[in]  dd 			DoF Distribution
	 */
	virtual	void prepare_timestep(ConstSmartPtr<VectorTimeSeries<vector_type> > vSol, GridLevel gl) = 0;

	///	prepares timestep on surface level
	void prepare_timestep(ConstSmartPtr<VectorTimeSeries<vector_type> > vSol)
		{prepare_timestep(vSol, GridLevel());}

	/// assembles Jacobian (or Approximation of Jacobian)
	/**
	 * Assembles Jacobian at a given Solution u.
	 *
	 * \param[out] J 			Jacobian J(u) Matrix to be filled
	 * \param[in]  vSol			vector of previous and current (iterated) solution
	 * \param[in]  s_a			scaling factors for mass matrix
	 * \param[in]  dd 			DoF Distribution
	 */
		virtual void assemble_jacobian(matrix_type& J,
		                               ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
		                               const number s_a, GridLevel gl) = 0;

	///	assembles jacobian on surface level
		void assemble_jacobian(matrix_type& J,
		                       ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
		                       const number s_a)
		{assemble_jacobian(J, vSol, s_a, GridLevel());}


	/// assembles Defect
	/**
	 * Assembles Defect at a given Solution u.
	 *
	 * \param[out] d 			Defect d(u) to be filled
	 * \param[in]  vSol			vector of previous and current (iterated) solution
	 * \param[in]  vScaleMass	scaling factors for mass matrix
	 * \param[in]  vScaleStiff	scaling factors for stiffness matrix
	 * \param[in]  dd 			DoF Distribution
	 */
		virtual	void assemble_defect(vector_type& d,
		       	                     ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
		       	                     const std::vector<number>& vScaleMass,
		       	                     const std::vector<number>& vScaleStiff,
		       	                     GridLevel gl) = 0;

	///	assembles defect on surface level
		void assemble_defect(vector_type& d,
		                     ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
		                     const std::vector<number>& vScaleMass,
		                     const std::vector<number>& vScaleStiff)
		{assemble_defect(d, vSol, vScaleMass, vScaleStiff, GridLevel());}

	/// Assembles matrix_type and Right-Hand-Side for a linear problem
	/**
	 * Assembles matrix_type and Right-Hand-Side for a linear problem
	 *
	 * \param[out] A 			Mass-/Stiffness- matrix_type of the discretization
	 * \param[out] b 			Right-Hand-Side of the discretization
	 * \param[in]  vSol			vector of previous and current (iterated) solution
	 * \param[in]  vScaleMass	scaling factors for mass matrix
	 * \param[in]  vScaleStiff	scaling factors for stiffness matrix
	 * \param[in]  dd 			DoF Distribution
	 */
		virtual void assemble_linear(matrix_type& A, vector_type& b,
		                             ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
		                             const std::vector<number>& vScaleMass,
		                             const std::vector<number>& vScaleStiff,
		                             GridLevel gl) = 0;

	///	assembles linear on surface level
		void assemble_linear(matrix_type& A, vector_type& b,
		                     ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
		                     const std::vector<number>& vScaleMass,
		                     const std::vector<number>& vScaleStiff)
		{assemble_linear(A,b,vSol,vScaleMass,vScaleStiff, GridLevel());}

	/// Assembles Right-Hand-Side for a linear problem
	/**
	 * Assembles Right-Hand-Side for a linear problem
	 *
	 * \param[out] b 			Right-Hand-Side of the discretization
	 * \param[in]  vSol			vector of previous and current (iterated) solution
	 * \param[in]  vScaleMass	scaling factors for mass matrix
	 * \param[in]  vScaleStiff	scaling factors for stiffness matrix
	 * \param[in]  dd 			DoF Distribution
	 */
		virtual void assemble_rhs(	 vector_type& b,
									 ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
									 const std::vector<number>& vScaleMass,
									 const std::vector<number>& vScaleStiff,
									 GridLevel gl) = 0;

	///	assembles rhs on surface level
		void assemble_rhs(vector_type& b,
							 ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
							 const std::vector<number>& vScaleMass,
							 const std::vector<number>& vScaleStiff)
		{assemble_rhs(b,vSol,vScaleMass,vScaleStiff, GridLevel());}


	/// sets dirichlet values in solution vector
	/**
	 * Sets dirichlet values of the Solution u when components are dirichlet
	 *
	 * \param[in]  u 		Solution to set values at
	 * \param[in]  time		time of next (to be computed) timestep
	 * \param[in]  dd 		DoF Distribution
	 */
		virtual void adjust_solution(vector_type& u, number time, GridLevel gl) = 0;

	///	adjust solution on surface level
		void adjust_solution(vector_type& u, number time)
		{adjust_solution(u,time, GridLevel());}

	/// finishes timestep
	/**
	 * finishes timestep at a given Solution u.
	 *
	 * \param[in]  vSol			vector of previous and current (iterated) solution
	 * \param[in]  dd 			DoF Distribution
	 */
		virtual	void finish_timestep(ConstSmartPtr<VectorTimeSeries<vector_type> > vSol, GridLevel gl) = 0;

	///	prepares timestep on surface level
		void finish_timestep(ConstSmartPtr<VectorTimeSeries<vector_type> > vSol)
			{finish_timestep(vSol, GridLevel());}

	///	returns the number of post processes
		virtual size_t num_dirichlet_constraints() const = 0;

	///	returns the i'th post process
		virtual SmartPtr<IConstraint<TAlgebra> > dirichlet_constraint(size_t i) = 0;
};

/// @}

}; // namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__DOMAIN_DISC_INTERFACE__ */
