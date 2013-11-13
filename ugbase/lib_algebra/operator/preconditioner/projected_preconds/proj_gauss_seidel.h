/*
 * proj_gauss_seidel.h
 *
 *  Created on: 10.10.2013
 *      Author: raphaelprohl
 */

#ifndef __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_PRECONDS__PROJ_GAUSS_SEIDEL__
#define __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_PRECONDS__PROJ_GAUSS_SEIDEL__

#include "proj_precond_interface.h"

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
	public IProjPreconditioner<TAlgebra>
{
	public:
	///	Base class type
		typedef IProjPreconditioner<TAlgebra> base_type;

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
		ProjGaussSeidel():IProjPreconditioner<TAlgebra>(){};

	///	name
		const char* name() const {return "Projected GaussSeidel";}

	///	Clone
		SmartPtr<ILinearIterator<vector_type> > clone()
		{
			SmartPtr<ProjGaussSeidel<TAlgebra> > newInst(
					new ProjGaussSeidel<TAlgebra>());
			newInst->set_damp(this->damping());
			return newInst;
		}

	///	computes a new correction c = B*d and projects on the underlying constraint
	/**
	 * This method computes a new correction c = B*d. B is here the underlying matrix operator.
	 * It can only be called, when the preprocess has been done.
	 *
	 * \param[out]	c			correction
	 * \param[in]	mat			underlying matrix (i.e. A in A*u = b)
	 * \param[in]	d			defect
	 */
		void projected_precond_step(vector_type& c, const matrix_type& mat, const vector_type& d);

	///	Destructor
		~ProjGaussSeidel(){};

	protected:
	/// flag indicating if an obstacle is set
		using base_type::m_bLowerObs;
		using base_type::m_bUpperObs;

	///	storage for last solution u
		using base_type::m_lastSol;

	///	relaxation parameter
		using base_type::m_relax;
};

} // end namespace ug

// include implementation
#include "proj_gauss_seidel_impl.h"

#endif /* __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_PRECONDS__PROJ_GAUSS_SEIDEL__ */
