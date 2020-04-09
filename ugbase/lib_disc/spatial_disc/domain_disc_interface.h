/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__DOMAIN_DISC_INTERFACE__
#define __H__UG__LIB_DISC__SPATIAL_DISC__DOMAIN_DISC_INTERFACE__

#include "lib_disc/assemble_interface.h"
#include "lib_disc/time_disc/solution_time_series.h"
#include "lib_disc/spatial_disc/constraints/constraint_interface.h"
#include "lib_grid/refinement/refiner_interface.h"
#include "lib_disc/function_spaces/error_elem_marking_strategy.h"

namespace ug {

/// \ingroup lib_disc_domain_assemble
/// @{

/// Interface for an object that can estimate the (global) error
template <typename TAlgebra>
class IDomainErrorIndicator
{
public:
	/// Type of algebra vector
	typedef typename TAlgebra::vector_type vector_type;

	// (virtual) destructor
	virtual ~IDomainErrorIndicator() {};

		/**
		 * Computes the error estimator.
		 *
		 * \param[in]  u			vector of the solution to estimate the error for
		 * \param[in]  dd 			DoF Distribution
		 */
			virtual	void calc_error
			(	const vector_type& u,
				const GridLevel& gl,
				vector_type* u_vtk = NULL
			) = 0;
			virtual void calc_error
			(	const vector_type& u,
				ConstSmartPtr<DoFDistribution> dd,
				vector_type* u_vtk = NULL
			) = 0;

			//! Transient version
			virtual void calc_error
			(	ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
				ConstSmartPtr<DoFDistribution> dd,
				const std::vector<number>& vScaleMass,
				const std::vector<number>& vScaleStiff,
				vector_type* u_vtk
			) = 0;

			//! Transient version
			virtual void calc_error
			(	ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
				const std::vector<number>& vScaleMass,
				const std::vector<number>& vScaleStiff,
				const GridLevel& gl,
				vector_type* u_vtk
			) = 0;

	/// marks error indicators as invalid
	virtual void invalidate_error() = 0;

	/// returns whether current error values are valid
	virtual bool is_error_valid() = 0;
};

/// Interface for an object that can mark elements based on a strategy
template <typename TDomain>
class IDomainMarker
{
public:
	/// Type of algebra vector
	typedef IElementMarkingStrategy<TDomain> element_marking_strategy_type;

	// (virtual) destructor
	virtual ~IDomainMarker() {};

	virtual void mark_with_strategy
	(	IRefiner& refiner,
		SmartPtr<element_marking_strategy_type> spMarkingStrategy
	) = 0 ;
};

/// Interface for domain discretization
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
class IDomainDiscretization : public IAssemble<TAlgebra>, public IDomainErrorIndicator<TAlgebra>
{
	using IAssemble<TAlgebra>::assemble_rhs;
	
	public:
	/// Algebra type
		typedef TAlgebra algebra_type;

	/// Type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

	/// Type of algebra vector
		typedef typename algebra_type::vector_type vector_type;
		virtual ~IDomainDiscretization() {};

	public:
		/// assembles Jacobian (or Approximation of Jacobian)
		/**
		 * Assembles Jacobian at a given iterate u.
		 *
		 * \param[out] 	J 	Jacobian J(u) matrix to be filled
		 * \param[in]  	u 	Current iterate
		 * \param[in]	dd	DoF Distribution
		 */
		virtual void assemble_jacobian(matrix_type& J, const vector_type& u, const GridLevel& gl) = 0;
		virtual void assemble_jacobian(matrix_type& J, const vector_type& u, ConstSmartPtr<DoFDistribution> dd) = 0;

		/// assembles Defect
		/**
		 * Assembles Defect at a given Solution u.
		 *
		 * \param[out] 	d 	Defect d(u) to be filled
		 * \param[in] 	u 	Current iterate
		 * \param[in]	dd	DoF Distribution
		 */
		virtual void assemble_defect(vector_type& d, const vector_type& u, const GridLevel& gl) = 0;
		virtual void assemble_defect(vector_type& d, const vector_type& u, ConstSmartPtr<DoFDistribution> dd) = 0;

		/// Assembles Matrix and Right-Hand-Side for a linear problem
		/**
		 * Assembles matrix_type and Right-Hand-Side for a linear problem
		 *
		 * \param[out] 	A 	Mass-/Stiffness- Matrix
		 * \param[out] 	b 	Right-Hand-Side
		 * \param[in]	dd	DoF Distribution
		 */
		virtual void assemble_linear(matrix_type& A, vector_type& b, const GridLevel& gl) = 0;
		virtual void assemble_linear(matrix_type& A, vector_type& b, ConstSmartPtr<DoFDistribution> dd) = 0;

		/// assembles the rhs
		virtual void assemble_rhs(vector_type& rhs, const vector_type& u, const GridLevel& gl) = 0;
		virtual void assemble_rhs(vector_type& rhs, const vector_type& u, ConstSmartPtr<DoFDistribution> dd) = 0;

		/// assembles the rhs
		virtual void assemble_rhs(vector_type& b, const GridLevel& gl) = 0;
		virtual void assemble_rhs(vector_type& b, ConstSmartPtr<DoFDistribution> dd) = 0;

		/// sets dirichlet values in solution vector
		/**
		 * Sets dirichlet values of the NumericalSolution u when components
		 * are dirichlet
		 *
		 * \param[out] 	u	Numerical Solution
		 * \param[in]	dd	DoF Distribution
		 */
		virtual void adjust_solution(vector_type& u, const GridLevel& gl) = 0;
		virtual void adjust_solution(vector_type& u, ConstSmartPtr<DoFDistribution> dd) = 0;

		/// assembles the mass matrix
		virtual void assemble_mass_matrix(matrix_type& M, const vector_type& u, const GridLevel& gl) = 0;
		virtual void assemble_mass_matrix(matrix_type& M, const vector_type& u, ConstSmartPtr<DoFDistribution> dd) = 0;

		/// assembles the stiffness matrix
		virtual void assemble_stiffness_matrix(matrix_type& A, const vector_type& u, const GridLevel& gl) = 0;
		virtual void assemble_stiffness_matrix(matrix_type& A, const vector_type& u, ConstSmartPtr<DoFDistribution> dd) = 0;


	public:
	/// prepares time step
	/**
	 * Prepares time step at a given solution u.
	 * This method is called only once before any time step.
	 *
	 * \param[in]  vSol			vector of previous and current (iterated) solution
	 * \param[in]  dd 			DoF distribution
	 */
		virtual void prepare_timestep(ConstSmartPtr<VectorTimeSeries<vector_type> > vSol, number future_time, ConstSmartPtr<DoFDistribution> dd) = 0;
		virtual void prepare_timestep(ConstSmartPtr<VectorTimeSeries<vector_type> > vSol, number future_time, const GridLevel& gl) = 0;
		virtual void prepare_timestep(ConstSmartPtr<VectorTimeSeries<vector_type> > vSol, number future_time)
		{prepare_timestep(vSol, future_time, GridLevel());}

	/// prepares time step element-wise
	/**
	 * prepares time step element-wise at a given solution u.
	 *
	 * \param[in]  vSol			vector of previous and current (iterated) solution
	 * \param[in]  dd 			DoF Distribution
	 */
		virtual void prepare_timestep_elem(ConstSmartPtr<VectorTimeSeries<vector_type> > vSol, ConstSmartPtr<DoFDistribution> dd) = 0;
		virtual void prepare_timestep_elem(ConstSmartPtr<VectorTimeSeries<vector_type> > vSol, const GridLevel& gl) = 0;

	///	prepares timestep on surface level
		void prepare_timestep_elem(ConstSmartPtr<VectorTimeSeries<vector_type> > vSol)
		{prepare_timestep_elem(vSol, GridLevel());}

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
		                               const number s_a, const GridLevel& gl) = 0;
		virtual void assemble_jacobian(matrix_type& J,
		                               ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
		                               const number s_a0,
		                               ConstSmartPtr<DoFDistribution> dd) = 0;

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
		       	                     const GridLevel& gl) = 0;
		virtual void assemble_defect(vector_type& d,
		                             ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
		                             const std::vector<number>& vScaleMass,
		                             const std::vector<number>& vScaleStiff,
		                             ConstSmartPtr<DoFDistribution> dd) = 0;

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
		                             const GridLevel& gl) = 0;
		virtual void assemble_linear(matrix_type& A, vector_type& b,
									 ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
									 const std::vector<number>& vScaleMass,
									 const std::vector<number>& vScaleStiff,
									 ConstSmartPtr<DoFDistribution> dd) = 0;

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
									 const GridLevel& gl) = 0;
		virtual void assemble_rhs(	 vector_type& b,
									 ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
									 const std::vector<number>& vScaleMass,
									 const std::vector<number>& vScaleStiff,
									 ConstSmartPtr<DoFDistribution> dd) = 0;

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
		virtual void adjust_solution(vector_type& u, number time, const GridLevel& gl) = 0;
		virtual void adjust_solution(vector_type& u, number time, ConstSmartPtr<DoFDistribution> dd) = 0;

	///	adjust solution on surface level
		void adjust_solution(vector_type& u, number time)
		{adjust_solution(u,time, GridLevel());}

	/// finishes time step
	/**
	 * Finishes time step at a given solution u.
	 * This method is called only once after any time step.
	 *
	 * \param[in]  vSol			vector of previous and current (iterated) solution
	 * \param[in]  dd 			DoF distribution
	 */
		virtual void finish_timestep(ConstSmartPtr<VectorTimeSeries<vector_type> > vSol, ConstSmartPtr<DoFDistribution> dd) = 0;
		virtual void finish_timestep(ConstSmartPtr<VectorTimeSeries<vector_type> > vSol, const GridLevel& gl) = 0;
		virtual void finish_timestep(ConstSmartPtr<VectorTimeSeries<vector_type> > vSol)
		{finish_timestep(vSol, GridLevel());}

	/// finishes timestep
	/**
	 * finishes timestep at a given Solution u.
	 *
	 * \param[in]  vSol			vector of previous and current (iterated) solution
	 * \param[in]  dd 			DoF Distribution
	 */
		virtual void finish_timestep_elem(ConstSmartPtr<VectorTimeSeries<vector_type> > vSol, const GridLevel& gl) = 0;
		virtual void finish_timestep_elem(ConstSmartPtr<VectorTimeSeries<vector_type> > vSol, ConstSmartPtr<DoFDistribution> dd) = 0;

	///	finishes timestep on surface level
		void finish_timestep_elem(ConstSmartPtr<VectorTimeSeries<vector_type> > vSol)
			{finish_timestep_elem(vSol, GridLevel());}




	///	returns the number of post processes
		virtual size_t num_constraints() const = 0;

	///	returns the i'th post process
		virtual SmartPtr<IConstraint<TAlgebra> > constraint(size_t i) = 0;
};

/// @}

}; // namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__DOMAIN_DISC_INTERFACE__ */
