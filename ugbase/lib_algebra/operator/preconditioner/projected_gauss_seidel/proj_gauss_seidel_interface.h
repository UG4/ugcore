/*
 * proj_gauss_seidel_interface.h
 *
 *  Created on: 13.11.2013
 *      Author: raphaelprohl
 */

#ifndef __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_GAUSS_SEIDEL__PROJ_GAUSS_SEIDEL_INTERFACE__
#define __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_GAUSS_SEIDEL__PROJ_GAUSS_SEIDEL_INTERFACE__

#include "obstacles/obstacle_constraint_interface.h"
#include "lib_algebra/operator/preconditioner/gauss_seidel.h"

namespace ug{

/// Interface for Projected GaussSeidel Preconditioner
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
template <typename TDomain, typename TAlgebra>
class IProjGaussSeidel:
	public GaussSeidelBase<TAlgebra>
{
	public:
	///	Base class type
		typedef GaussSeidelBase<TAlgebra> base_type;

	///	Algebra type
		typedef TAlgebra algebra_type;

	///	Matrix type
		typedef typename algebra_type::matrix_type matrix_type;

	///	Vector type
		typedef typename algebra_type::vector_type vector_type;

	///	Value type
		typedef typename vector_type::value_type value_type;

	///	Grid Function type
		typedef GridFunction<TDomain, TAlgebra> GF;

	public:
	/// constructor
		IProjGaussSeidel(): GaussSeidelBase<TAlgebra>(){
			m_spvObsConstraint.clear();
			m_bObsCons = false;
		};


	/// clone constructor
		IProjGaussSeidel( const IProjGaussSeidel<TDomain, TAlgebra> &parent )
			: base_type(parent)
		{
			m_spvObsConstraint = parent.m_spvObsConstraint;
			m_bObsCons = parent.m_bObsCons;
		}

	///	adds the obstacle constraint function c(u)
		void add_obstacle_constraint(SmartPtr<IObstacleConstraint<TDomain,TAlgebra> > spObsCons)
		{
			m_spvObsConstraint.push_back(spObsCons);
			m_bObsCons = true;

			//	inits the obstacle constraint
			spObsCons->init();
		}

	///	Destructor
		~IProjGaussSeidel(){};

	///	name
		virtual const char* name() const = 0;

	/// Prepare for Operator J(u) and linearization point u (current solution)
		virtual bool init(SmartPtr<ILinearOperator<vector_type> > J, const vector_type& u);

	///	computes a new correction c = B*d and projects on the underlying constraint
	/**
	 * This method computes a new correction c = B*d. B is here the underlying matrix operator.
	 *
	 * \param[in]	mat			underlying matrix (i.e. A in A*u = b)
	 * \param[out]	c			correction
	 * \param[in]	d			defect
	 * \param[in]	relax		relaxation parameter
	 */
		virtual void step(const matrix_type& mat, vector_type& c, const vector_type& d, const number relax) = 0;

	///	projects the correction on the underlying constraints set by the obstacleConstraints
		void project_correction(value_type& c_i, const size_t i);

	///	Compute new correction c = B*d
		virtual bool apply(vector_type& c, const vector_type& d);

	///	Compute new correction c = B*d and return new defect d := d - A*c
		virtual bool apply_update_defect(vector_type& c, vector_type& d);

	private:
	///	for all indices stored in vInd:
	///	the entry of vec is set to zero
		void truncateVec(vector_type& vec, vector<DoFIndex>& vInd);

	///	for all indices stored in vInd:
	///	all rows and columns of mat are set to zero
		void truncateMat(matrix_type& mat, vector<DoFIndex>& vInd);

	protected:
	///	obstacle constraint
		vector<SmartPtr<IObstacleConstraint<TDomain,TAlgebra> > > m_spvObsConstraint;

	private:
	///	pointer to solution
		SmartPtr<vector_type> m_spSol;

	/// flag indicating if obstacle constraint has been set
		bool m_bObsCons;

	/// init flag indicating if init has been called
		bool m_bInit;

};


} // end namespace ug

// include implementation
#include "proj_gauss_seidel_interface_impl.h"

#endif /* __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_GAUSS_SEIDEL__PROJ_GAUSS_SEIDEL_INTERFACE__ */
