/*
 * proj_precond_interface.h
 *
 *  Created on: 13.11.2013
 *      Author: raphaelprohl
 */

#ifndef __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_PRECONDS__PROJ_PRECOND_INTERFACE__
#define __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_PRECONDS__PROJ_PRECOND_INTERFACE__

#include "lib_algebra/operator/interface/linear_iterator.h"

namespace ug{

/// Interface for Projected Preconditioners
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
			m_bLowerObs(false), m_bUpperObs(false), m_spMat(NULL), m_bInit(false){
			m_vActiveIndicesLow.resize(0); m_vActiveIndicesUp.resize(0);
			m_vInactiveIndices.resize(0);
			m_relax = 1.0;};

	///	preprocess checks if matrix is diagonal invertible
		bool preprocess(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp);

	///	set constraint/obstacle g_low (for u >= g_low)
		void set_lower_obstacle(SmartPtr<vector_type> lowObs){
			m_spVecOfLowObsValues = lowObs; m_bLowerObs = true;}

	///	set constraint/obstacle g_up (for u <= g_up)
		void set_upper_obstacle(SmartPtr<vector_type> upObs){
			m_spVecOfUpObsValues = upObs; m_bUpperObs = true;}

	///	set relaxation parameter to define a SOR-method
		void set_sor_relax(number relaxFactor){ m_relax = relaxFactor;}

	///	computes a new correction c = B*d and projects on the underlying constraint
	/**
	 * This method computes a new correction c = B*d. B is here the underlying matrix operator.
	 * It can only be called, when the preprocess has been done.
	 *
	 * \param[out]	c			correction
	 * \param[in]	mat			underlying matrix (i.e. A in A*u = b)
	 * \param[in]	d			defect
	 */
		virtual void projected_precond_step(vector_type& c, const matrix_type& mat, const vector_type& d) = 0;

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

	protected:
	///	computes the correction for the case that only a lower obstacle is set, i.e. u >= g_low
		void correction_for_lower_obs(vector_type& c, const size_t index, const value_type tmpSol);

	///	computes the correction for the case that only an upper obstacle is set, i.e. u <= g_up
		void correction_for_upper_obs(vector_type& c, const size_t index, const value_type tmpSol);

	///	computes the correction for the case that a lower and an upper obstacle is set
		void correction_for_lower_and_upper_obs(vector_type& c, const size_t index, const value_type tmpSol);

	private:
	///	adjust defect of the active indices for the case that a constraint/obstacle is set
		void adjust_defect(vector_type& d);

	protected:
	/// flag indicating if an obstacle is set
		bool m_bLowerObs, m_bUpperObs;

	///	storage for last solution u
		vector_type m_lastSol;

	///	relaxation parameter
		number m_relax;

	private:
	///	pointer to constraint/obstacle values
		SmartPtr<vector_type> m_spVecOfLowObsValues, m_spVecOfUpObsValues;

	/// operator to invert
		SmartPtr<matrix_type> m_spMat;

	/// init flag indicating if init has been called
		bool m_bInit;

	///	store the indices, which satisfy the constraints with equality in m_vActiveIndices.
	///	Other indices are stored in m_vInactiveIndices.
		std::vector<size_t> m_vActiveIndicesLow, m_vActiveIndicesUp, m_vInactiveIndices;
};

} // end namespace ug

// include implementation
#include "proj_precond_interface_impl.h"

#endif /* __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_PRECONDS__PROJ_PRECOND_INTERFACE__ */
