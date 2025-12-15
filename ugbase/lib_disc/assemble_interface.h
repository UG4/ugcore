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

#ifndef __H__UG__LIB_DISC__ASSEMBLE__
#define __H__UG__LIB_DISC__ASSEMBLE__

// #include "lib_grid/tools/selector_grid.h"
#include "lib_grid/tools/grid_level.h"
#include "lib_disc/spatial_disc/ass_tuner.h"

namespace ug {

//predeclaration
template <typename TAlgebra>
class IConstraint;


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
 *	\tparam		TAlgebra				Algebra type
 */
template <typename TAlgebra>
class IAssemble
{
	public:
	///	Algebra type
	using algebra_type = TAlgebra;

	///	Type of algebra matrix
	using matrix_type = typename TAlgebra::matrix_type;

	///	Type of algebra vector
	using vector_type = typename TAlgebra::vector_type;

	public:
		/// assembles Jacobian (or Approximation of Jacobian)
		/**
		 * Assembles Jacobian at a given iterate u.
		 *
		 * \param[out] 	J 	Jacobian J(u) matrix to be filled
		 * \param[in]  	u 	Current iterate
		 * \param[in]	gl	Grid Level
		 */
		virtual void assemble_jacobian(matrix_type& J, const vector_type& u, const GridLevel& gl) = 0;
		void assemble_jacobian(matrix_type& J, const vector_type& u)
		{assemble_jacobian(J,u,GridLevel());}

		/// assembles Defect
		/**
		 * Assembles Defect at a given Solution u.
		 *
		 * \param[out] 	d 	Defect d(u) to be filled
		 * \param[in] 	u 	Current iterate
		 * \param[in]	gl	Grid Level
		 */
		virtual void assemble_defect(vector_type& d, const vector_type& u, const GridLevel& gl) = 0;
		void assemble_defect(vector_type& d, const vector_type& u)
		{assemble_defect(d,u, GridLevel());}

		/// Assembles Matrix and Right-Hand-Side for a linear problem
		/**
		 * Assembles matrix_type and Right-Hand-Side for a linear problem
		 *
		 * \param[out] 	A 	Mass-/Stiffness- Matrix
		 * \param[out] 	b 	Right-Hand-Side
		 * \param[in]	gl	Grid Level
		 */
		virtual void assemble_linear(matrix_type& A, vector_type& b, const GridLevel& gl) = 0;
		void assemble_linear(matrix_type& A, vector_type& b)
		{assemble_linear(A,b, GridLevel());}

	///	assembles rhs
		virtual void assemble_rhs(vector_type& rhs, const vector_type& u, const GridLevel& gl) = 0;
		virtual void assemble_rhs(vector_type& rhs, const vector_type& u)
		{assemble_rhs(rhs, u, GridLevel());}

		/// Assembles Right-Hand-Side for a linear problem
		/**
		 * Assembles Right-Hand-Side for a linear problem
		 *
		 * \param[out] 	b 	Right-Hand-Side
		 * \param[in]	gl	Grid Level
		 */
		virtual void assemble_rhs(vector_type& b, const GridLevel& gl) = 0;
		void assemble_rhs(vector_type& b)
		{assemble_rhs(b, GridLevel());}

		/// sets dirichlet values in solution vector
		/**
		 * Sets dirichlet values of the NumericalSolution u when components
		 * are dirichlet
		 *
		 * \param[out] 	u	Numerical Solution
		 * \param[in]	gl	Grid Level
		 */
		virtual void adjust_solution(vector_type& u, const GridLevel& gl) = 0;
		void adjust_solution(vector_type& u)
		{adjust_solution(u, GridLevel());}

	///	assembles mass matrix
		virtual void assemble_mass_matrix(matrix_type& M, const vector_type& u, const GridLevel& gl)
		{UG_THROW("IAssemble: assemble_mass_matrix not implemented.");}
		void assemble_mass_matrix(matrix_type& M, const vector_type& u)
		{assemble_mass_matrix(M,u,GridLevel());}

	///	assembles stiffness matrix
		virtual void assemble_stiffness_matrix(matrix_type& A, const vector_type& u, const GridLevel& gl)
		{UG_THROW("IAssemble: assemble_stiffness_matrix not implemented.");}
		void assemble_stiffness_matrix(matrix_type& A, const vector_type& u)
		{assemble_stiffness_matrix(A,u,GridLevel());}

	/// \{
		virtual SmartPtr<AssemblingTuner<TAlgebra> > ass_tuner() = 0;
		[[nodiscard]] virtual ConstSmartPtr<AssemblingTuner<TAlgebra> > ass_tuner() const = 0;
	/// \}

	///	returns the number of constraints
		[[nodiscard]] virtual size_t num_constraints() const = 0;

	///	returns the i'th constraint
		virtual SmartPtr<IConstraint<TAlgebra> > constraint(size_t i) = 0;

	/// Virtual Destructor
		virtual ~IAssemble() = default;

};

/// @}

} // name space ug

#endif
