/*
 * proj_gauss_seidel.h
 *
 *  Created on: 10.10.2013
 *      Author: raphaelprohl
 */

#ifndef __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_GAUSS_SEIDEL__PROJ_GAUSS_SEIDEL__
#define __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_GAUSS_SEIDEL__PROJ_GAUSS_SEIDEL__

#include "proj_gauss_seidel_interface.h"

namespace ug{

/// Projected GaussSeidel (SOR) -method
/**
 * The projected GaussSeidel method can be applied to problems of the form
 *
 * 		A * u >= b				(I)
 * 		c(u) >= 0				(II)
 * 		c(u)^T * [A*u - b] = 0,	(III)
 *
 * where u, b are vectors and A is a matrix. '*' denotes componentwise multiplication.
 * c(u) denotes an obstacle-function, which depends on the solution vector u. One possible
 * example for such an obstacle-function could be the scalar obstacle function
 *
 * 	u >= 0.
 *
 * The obstacle function c(u) is defined by creating an instance of IObstacleConstraint, which is
 * passed to the projected preconditioner by the method 'IProjPreconditioner::set_obstacle_constraint'.
 *
 * Similar problems, which e.g. only differ in the sign in (I) and/or (II) can be
 * equivalently treated by the method.
 *
 * Note: Due to (II) the old solution needs to be stored within this method.
 * This is a difference to the classical smoothers/preconditioners, which usually work
 * on the correction and the defect only.
 *
 * By calling 'set_sor_damp(number)' one gets the successive overrelaxation-version of the
 * projected preconditioners of GaussSeidel type.
 *
 * References:
 * <ul>
 * <li> A. Brandt and C. W. Cryer. Multigrid algorithms for the solution of linear
 * <li>	complementarity problems arising from free boundary problems. SIAM J. SCI. STAT. COMPUT. Vol. 4, No. 4 (1983)
 * </ul>
 *
 * \tparam 	TAlgebra		Algebra type
 */

template <typename TDomain, typename TAlgebra>
class ProjGaussSeidel:
	public IProjGaussSeidel<TDomain,TAlgebra>
{
	public:
	///	Base class type
		typedef IProjGaussSeidel<TDomain,TAlgebra> base_type;

	///	Algebra type
		typedef TAlgebra algebra_type;

	///	Matrix type
		typedef typename algebra_type::matrix_type matrix_type;

	///	Vector type
		typedef typename algebra_type::vector_type vector_type;

	public:
	/// constructor
		ProjGaussSeidel(): IProjGaussSeidel<TDomain,TAlgebra>(){};

	///	Destructor
		~ProjGaussSeidel(){};

	protected:
	///	name
		virtual const char* name() const {return "Projected GaussSeidel";}

	///	Clone
		SmartPtr<ILinearIterator<vector_type> > clone()
		{
			SmartPtr<ProjGaussSeidel<TDomain,TAlgebra> > newInst(
					new ProjGaussSeidel<TDomain,TAlgebra>());
			base_type::copy_config(*newInst);
			return newInst;
		}

	///	computes a new correction c = B*d and projects on the underlying constraint
		virtual void step(const matrix_type& mat, vector_type& c, const vector_type& d, const number relax);
};

template <typename TDomain, typename TAlgebra>
class ProjBackwardGaussSeidel:
	public IProjGaussSeidel<TDomain,TAlgebra>
{
	public:
	///	Base class type
		typedef IProjGaussSeidel<TDomain,TAlgebra> base_type;

	///	Algebra type
		typedef TAlgebra algebra_type;

	///	Matrix type
		typedef typename algebra_type::matrix_type matrix_type;

	///	Vector type
		typedef typename algebra_type::vector_type vector_type;

	public:
	/// constructor
		ProjBackwardGaussSeidel(): IProjGaussSeidel<TDomain,TAlgebra>(){};

	///	Destructor
		~ProjBackwardGaussSeidel(){};

	protected:
	///	name
		virtual const char* name() const {return "Projected Backward GaussSeidel";}

	///	Clone
		SmartPtr<ILinearIterator<vector_type> > clone()
		{
			SmartPtr<ProjBackwardGaussSeidel<TDomain,TAlgebra> > newInst(
					new ProjBackwardGaussSeidel<TDomain,TAlgebra>());
			base_type::copy_config(*newInst);
			return newInst;
		}

	///	computes a new correction c = B*d and projects on the underlying constraint
		virtual void step(const matrix_type& mat, vector_type& c, const vector_type& d, const number relax);
};


template <typename TDomain, typename TAlgebra>
class ProjSymmetricGaussSeidel:
	public IProjGaussSeidel<TDomain,TAlgebra>
{
	public:
	///	Base class type
		typedef IProjGaussSeidel<TDomain,TAlgebra> base_type;

	///	Algebra type
		typedef TAlgebra algebra_type;

	///	Matrix type
		typedef typename algebra_type::matrix_type matrix_type;

	///	Vector type
		typedef typename algebra_type::vector_type vector_type;

	public:
	/// constructor
		ProjSymmetricGaussSeidel(): IProjGaussSeidel<TDomain,TAlgebra>(){};

	///	Destructor
		~ProjSymmetricGaussSeidel(){};

	protected:
	///	name
		virtual const char* name() const {return "Projected Symmetric GaussSeidel";}

	///	Clone
		SmartPtr<ILinearIterator<vector_type> > clone()
		{
			SmartPtr<ProjSymmetricGaussSeidel<TDomain,TAlgebra> > newInst(
					new ProjSymmetricGaussSeidel<TDomain,TAlgebra>());
			base_type::copy_config(*newInst);
			return newInst;
		}

	///	computes a new correction c = B*d and projects on the underlying constraint
		virtual void step(const matrix_type& mat, vector_type& c, const vector_type& d, const number relax);
};

} // end namespace ug

// include implementation
#include "proj_gauss_seidel_impl.h"

#endif /* __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_GAUSS_SEIDEL__PROJ_GAUSS_SEIDEL__ */
