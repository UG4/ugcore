/*
 * proj_gauss_seidel.h
 *
 *  Created on: 10.10.2013
 *      Author: raphaelprohl
 */

#ifndef __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJ_GAUSS_SEIDEL__PROJ_GAUSS_SEIDEL__
#define __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJ_GAUSS_SEIDEL__PROJ_GAUSS_SEIDEL__

namespace ug{

/// Projected GaussSeidel-method
/**
 * The projected GaussSeidel method can be applied to problems of the form
 *
 * 		A * u <= b				(I)
 * 		u >= 0					(II)
 * 		u^T * [A*u - b] = 0,	(III)
 *
 * where u, b are vectors and A is a matrix. '*' denotes componentwise multiplication.
 * Similar problems, which e.g. only differ in the sign in (I) and/or (II) can be
 * equivalently treated by the method.
 *
 * Note: Due to (II) the old solution needs to stored within this method.
 * This is a difference to the classical smoothers/preconditioners, which usually work
 * on the correction and the defect only.
 *
 * Note: At the moment the method is only implemented for the 'forward' ordering of the dofs
 *
 * References:
 * <ul>
 * <li> A. Brandt and C. W. Cryer. Multigrid algorithms for the solution of linear
 * <li>	complementarity problems arising from free boundary problems. SIAM J. SCI. STAT. COMPUT. Vol. 4, No. 4 (1983)
 * </ul>
 *
 * \tparam 	TAlgebra		Algebra type
 */
template <typename TAlgebra>
class ProjGaussSeidel:
	public ILinearIterator<typename TAlgebra::vector_type>
{

	public:
	///	Algebra type
		typedef TAlgebra algebra_type;

	///	Vector type
		typedef typename algebra_type::vector_type vector_type;

	///	Matrix type
		typedef typename algebra_type::matrix_type matrix_type;

	public:
	/// constructor
		ProjGaussSeidel(): m_spMat(NULL), m_bInit(false){ m_vActiveIndices.resize(0); };

	///	preprocess checks if matrix is diagonal invertible
		bool preprocess(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp);

	///	computes a new correction c = B*d and projects on the underlying constraint
	/**
	 * This method computes a new correction c = B*d. B is here the underlying matrix operator.
	 * It can only be called, when the preprocess has been done.
	 *
	 * \param[out]	c			correction
	 * \param[in]	mat			underlying matrix (i.e. A in A*u = b)
	 * \param[in]	d			defect
	 * \returns		bool		success flag
	 */
		bool gs_step_with_projection(vector_type& c, const matrix_type& mat, const vector_type& d);

	///////////////////////////////////////////////////////////////////////////
	//	Linear Solver interface methods
	///////////////////////////////////////////////////////////////////////////

	///	name
		virtual const char* name() const {return "Projected GaussSeidel";}

	/// Prepare for Operator J(u) and linearization point u (current solution)
		virtual bool init(SmartPtr<ILinearOperator<vector_type> > J, const vector_type& u);

	///	Prepare for Linear Operartor L
		virtual bool init(SmartPtr<ILinearOperator<vector_type> > L);

	///	Compute new correction c = B*d
		virtual bool apply(vector_type& c, const vector_type& d);

	///	Compute new correction c = B*d and return new defect d := d - A*c
		virtual bool apply_update_defect(vector_type& c, vector_type& d);

	///	Clone
		SmartPtr<ILinearIterator<vector_type> > clone();

		///	Destructor
		~ProjGaussSeidel()
		{};

	private:
	/// operator to invert
		SmartPtr<matrix_type> m_spMat;

	///	storage for last solution u
		vector_type m_lastSol;

	/// init flag indicating if init has been called
		bool m_bInit;

	///	storage for the indices, which satisfy the constraints with equality
		std::vector<size_t> m_vActiveIndices;
};

} // end namespace ug

// include implementation
#include "proj_gauss_seidel_impl.h"

#endif /* __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJ_GAUSS_SEIDEL__PROJ_GAUSS_SEIDEL__ */
