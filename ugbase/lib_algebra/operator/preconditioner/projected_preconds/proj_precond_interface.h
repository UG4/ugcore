/*
 * proj_precond_interface.h
 *
 *  Created on: 13.11.2013
 *      Author: raphaelprohl
 */

#ifndef __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_PRECONDS__PROJ_PRECOND_INTERFACE__
#define __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_PRECONDS__PROJ_PRECOND_INTERFACE__

#include "lib_algebra/operator/interface/linear_iterator.h"

#include "obstacles/obstacle_constraint_interface.h"

namespace ug{

/// Interface for Projected Preconditioners
/**
 * 	This class provides an interface to define a preconditioner which can be applied to solve
 * 	problems of the form
 *
 * 		A * u >= b				(I)
 * 		c(u) >= 0				(II)
 * 		c(u)^T * [A*u - b] = 0,	(III)
 *
 * 	where u, b are vectors and A is a matrix. '*' denotes componentwise multiplication.
 * 	c(u) denotes an obstacle-function, which depends on the solution vector u. One possible
 * 	example for such an obstacle-function could be the scalar obstacle function
 *
 * 		u >= 0.
 *
 * 	The obstacle function c(u) is defined by creating an instance of IObstacleConstraint, which is
 * 	passed to the projected preconditioner by the method 'set_obstacle_constraint'.
 *
 *	Similar problems, which e.g. only differ in the sign in (I) and/or (II) can be
 * 	equivalently treated by these preconditioners.
 *
 * 	Note: Due to (II) the old solution needs to be stored within this method.
 *	This is a difference to the classical smoothers/preconditioners, which usually work
 *	on the correction and the defect only.
 *
 * 	Since the problem formulation (I)-(III) consists of inequalities, the projected preconditioner
 * 	performs a projection on a constraint c(u) in every preconditioner-step.
 *
 *  \tparam 	TAlgebra		Algebra type
 */
template <typename TAlgebra>
class IProjPreconditioner:
	public ILinearIterator<typename TAlgebra::vector_type>
{
	public:
	///	Algebra type
		typedef TAlgebra algebra_type;

	///	Matrix type
		typedef typename algebra_type::matrix_type matrix_type;

	///	Vector type
		typedef typename algebra_type::vector_type vector_type;

	///	Value type
		typedef typename vector_type::value_type value_type;

	public:
	/// constructor
		IProjPreconditioner(): ILinearIterator<typename TAlgebra::vector_type>(),
			m_bObsCons(false), m_spMat(NULL), m_bInit(false){
			m_relax = 1.0;};

	///	preprocess checks if matrix is diagonal invertible
		bool preprocess(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp);

	///	sets the obstacle constraint function c(u)
		void set_obstacle_constraint(SmartPtr<IObstacleConstraint<TAlgebra> > spObsCons){
			m_spObsConstraint = spObsCons; m_bObsCons = true;}

	///	set relaxation parameter to define a SOR-method
		void set_sor_relax(number relaxFactor){ m_relax = relaxFactor;}

	///	computes a new correction c = B*d and projects on the underlying constraint.
	/**
	 * This method computes a new correction c = B*d. The matrix operator B is defined in the
	 * particular inheritated sub-class. It can only be called, when the preprocess has been done.
	 *
	 * \param[out]	c			correction
	 * \param[out]	sol			solution
	 * \param[in]	mat			underlying matrix (i.e. A in A*u = b)
	 * \param[in]	d			defect
	 */
		virtual void projected_precond_step(vector_type& c, vector_type& sol, const matrix_type& mat, const vector_type& d) = 0;

	///////////////////////////////////////////////////////////////////////////
	//	Linear Solver interface methods
	///////////////////////////////////////////////////////////////////////////

	/// Prepare for Operator J(u) and linearization point u (current solution)
		bool init(SmartPtr<ILinearOperator<vector_type> > J, const vector_type& u);

	///	Prepare for Linear Operartor L
		bool init(SmartPtr<ILinearOperator<vector_type> > L);

	///	Compute new correction c = B*d
		bool apply(vector_type& c, const vector_type& d);

	///	Compute new correction c = B*d and return new defect d := d - A*c
		bool apply_update_defect(vector_type& c, vector_type& d);

	///	Destructor
		virtual ~IProjPreconditioner(){};

	private:
	///	adjust defect of the active indices for the case that a constraint/obstacle is set
		void adjust_defect(vector_type& d);


	protected:
	///	storage for last solution u
		SmartPtr<vector_type> m_lastSol;

	///	relaxation parameter
		number m_relax;

	/// point to obstacle constraint
		SmartPtr<IObstacleConstraint<TAlgebra> > m_spObsConstraint;

	private:
	/// flag indicating if obstacle constraint has been set
		bool m_bObsCons;

	/// operator to invert
		SmartPtr<MatrixOperator<matrix_type, vector_type> > m_spMat;

	/// init flag indicating if init has been called
		bool m_bInit;

};

} // end namespace ug

// include implementation
#include "proj_precond_interface_impl.h"

#endif /* __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_PRECONDS__PROJ_PRECOND_INTERFACE__ */
