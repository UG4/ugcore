/*
 * assemble.h
 *
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__ASSEMBLE__
#define __H__UG__LIB_DISC__ASSEMBLE__

#include "lib_grid/tools/bool_marker.h"
#include "lib_grid/tools/selector_grid.h"
#include "lib_disc/dof_manager/grid_level.h"

namespace ug{

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
		typedef TAlgebra algebra_type;

	///	Type of algebra matrix
		typedef typename TAlgebra::matrix_type matrix_type;

	///	Type of algebra vector
		typedef typename TAlgebra::vector_type vector_type;

	public:
		/// assembles Jacobian (or Approximation of Jacobian)
		/**
		 * Assembles Jacobian at a given iterate u.
		 *
		 * \param[out] 	J 	Jacobian J(u) matrix to be filled
		 * \param[in]  	u 	Current iterate
		 * \param[in]	dd	DoF Distribution
		 */
		virtual void assemble_jacobian(matrix_type& J, const vector_type& u,
		                               GridLevel gl) = 0;

		///	assembles jacobian on surface grid
		void assemble_jacobian(matrix_type& J, const vector_type& u)
			{assemble_jacobian(J,u,GridLevel());}

		/// assembles Defect
		/**
		 * Assembles Defect at a given Solution u.
		 *
		 * \param[out] 	d 	Defect d(u) to be filled
		 * \param[in] 	u 	Current iterate
		 * \param[in]	dd	DoF Distribution
		 */
		virtual void assemble_defect(vector_type& d, const vector_type& u,
		                             GridLevel gl) = 0;

		/// assembles jacobian on surface grid
		void assemble_defect(vector_type& d, const vector_type& u)
			{assemble_defect(d,u, GridLevel());}

		/// Assembles Matrix and Right-Hand-Side for a linear problem
		/**
		 * Assembles matrix_type and Right-Hand-Side for a linear problem
		 *
		 * \param[out] 	A 	Mass-/Stiffness- Matrix
		 * \param[out] 	b 	Right-Hand-Side
		 * \param[in]	dd	DoF Distribution
		 */
		virtual void assemble_linear(matrix_type& A, vector_type& b,
		                             GridLevel gl) = 0;

		/// assembles linear on the surface grid
		void assemble_linear(matrix_type& A, vector_type& b)
			{assemble_linear(A,b, GridLevel());}

	///	assembles rhs
		virtual void assemble_rhs(vector_type& rhs, const vector_type& u,
								   GridLevel gl) = 0;

	///	assembles rhs on surface grid
		virtual void assemble_rhs(vector_type& rhs, const vector_type& u)
			{assemble_rhs(rhs, u, GridLevel());}

		/// Assembles Right-Hand-Side for a linear problem
		/**
		 * Assembles Right-Hand-Side for a linear problem
		 *
		 * \param[out] 	b 	Right-Hand-Side
		 * \param[in]	dd	DoF Distribution
		 */
		virtual void assemble_rhs(vector_type& b, GridLevel gl) = 0;

		/// assembles linear on the surface grid
		void assemble_rhs(vector_type& b)
			{assemble_rhs(b, GridLevel());}

		/// sets dirichlet values in solution vector
		/**
		 * Sets dirichlet values of the NumericalSolution u when components
		 * are dirichlet
		 *
		 * \param[out] 	u	Numerical Solution
		 * \param[in]	dd	DoF Distribution
		 */
		virtual void adjust_solution(vector_type& u,
		                             GridLevel gl) = 0;

	/// adjust solution on surface grid
		void adjust_solution(vector_type& u)
			{adjust_solution(u, GridLevel());}

	///	assembles mass matrix
		virtual void assemble_mass_matrix(matrix_type& M, const vector_type& u,
		                               GridLevel gl)
		{UG_THROW("IAssemble: assemble_mass_matrix not implemented.");}

	///	assembles mass matrix on surface grid
		void assemble_mass_matrix(matrix_type& M, const vector_type& u)
			{assemble_mass_matrix(M,u,GridLevel());}

	///	assembles stiffness matrix
		virtual void assemble_stiffness_matrix(matrix_type& A, const vector_type& u,
		                                       GridLevel gl)
		{UG_THROW("IAssemble: assemble_stiffness_matrix not implemented.");}

	///	assembles stiffness matrix on surface grid
		void assemble_stiffness_matrix(matrix_type& A, const vector_type& u)
			{assemble_stiffness_matrix(A,u,GridLevel());}

	/// forces the assembling to consider the grid as regular
		virtual void force_regular_grid(bool bForce) = 0;

	///	returns type of constraints enabled
		virtual int constraints_enabled() const = 0;

	///	enables constraints
		virtual void enable_constraints(int TypesEnable) = 0;

	///	returns type of boundary elem discs enabled
		virtual int elem_discs_enabled() const = 0;

	///	enables boundary elem discs
		virtual void enable_elem_discs(int TypesEnable) = 0;

	///	sets a marker to exclude elements from assembling
	/**
	 * This method sets a marker. Only elements that are marked will be
	 * assembled during assembling process. If no marker is set, this
	 * corresponds to a marker where all elements have been marked.
	 *
	 * \param[in]	mark	BoolMarker
	 */
		virtual void set_marker(BoolMarker* mark = NULL) = 0;

	///	sets a selector of elements for assembling
	/**
	 * This method sets a selector. Only elements of this selector will be
	 * assembled during assembling process. Especially the selector defines the begin
	 * and end of the element-iterator in the element assembling-loop.
	 * If no selector is set, this corresponds to a assembling where the loop is
	 * carried out over all elements of a given subset.
	 *
	 * \param[in]	sel		Selector
	 */
		virtual void set_selector(Selector* sel = NULL) = 0;
	
		virtual void ass_index(size_t ind, bool index_set = true) = 0;

	///	returns the number of constraints
		virtual size_t num_constraints() const = 0;

	///	returns the i'th constraint
		virtual SmartPtr<IConstraint<TAlgebra> > constraint(size_t i) = 0;

	/// Virtual Destructor
		virtual ~IAssemble(){};
};

/// @}

}; // name space ug

#endif /* __H__UG__LIB_DISC__ASSEMBLE__ */
