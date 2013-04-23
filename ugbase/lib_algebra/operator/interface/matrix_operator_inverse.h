/*
 * matrix_operator_inverse.h
 *
 *  Created on: 23.04.2013
 *      Author: andreasvogel
 */

#ifndef __H__LIB_ALGEBRA__OPERATOR__INTERFACE__MATRIX_OPERATOR_INVERSE__
#define __H__LIB_ALGEBRA__OPERATOR__INTERFACE__MATRIX_OPERATOR_INVERSE__

#include "linear_operator_inverse.h"
#include "matrix_operator.h"
#include "common/error.h"
#include "common/util/smart_pointer.h"

namespace ug{

///////////////////////////////////////////////////////////////////////////////
// Inverse of a Matrix-based Linear Operator
///////////////////////////////////////////////////////////////////////////////

/// describes an inverse linear mapping X->Y based on a matrix
/**
 * This class is the base class for the inversion of linear matrix-based operator
 * given in form of the IMatrixOperator interface class. Given a operator L,
 * the basic usage of this class is to invert this operator, i.e. to compute
 * the solution u of
 *
 * 		L*u = f     i.e. u := L^{-1} f
 *
 *
 * \tparam	X		domain space (i.e. a vector corresponding to the matrix)
 * \tparam	Y		range space (i.e. a vector corresponding to the matrix)
 * \tparam	M		matrix type used to represent linear mapping
 */
template <typename M, typename X, typename Y = X>
class IMatrixOperatorInverse
	: public virtual ILinearOperatorInverse<X,Y>
{
	public:
	///	Domain space
		typedef X domain_function_type;

	///	Range space
		typedef Y codomain_function_type;

	///	Matrix type
		typedef M matrix_type;

	public:
	///	initializes this inverse operator for a matrix-based operator
	/**
	 * This method passes the operator A that is inverted by this operator. In
	 * addition some preparation step can be made.
	 *
	 * \param[in]	A		linear matrix-basewd operator to invert
	 * \returns		bool	success flag
	 */
		virtual bool init(SmartPtr<MatrixOperator<M,Y,X> > A) = 0;

	/// applies the inverse operator, i.e. returns u = A^{-1} * f
	/**
	 * This method applies the inverse operator.
	 *
	 * \param[out]	u		solution
	 * \param[in]	f		right-hand side
	 * \returns		bool	success flag
	 */
		virtual bool apply(Y& u, const X& f) = 0;

	/// applies the inverse operator and updates the defect
	/**
	 * This method applies the inverse operator and updates the defect, i.e.
	 * returns u = A^{-1} * f and in f the last defect d:= f - A*u is returned.
	 *
	 * \param[out]	u		solution
	 * \param[in]	f		right-hand side on entry, defect on exit
	 * \returns		bool	success flag
	 */
		virtual bool apply_return_defect(Y& u, X& f) = 0;

	/// virtual destructor
		virtual ~IMatrixOperatorInverse() {};

	public:
	///	initializes this inverse operator for a linear operator
	/**
	 * This method implements the ILinearOperatorInverse interface method.
	 * Basically, the request is forwarded to the matrix-based init method,
	 * if the the operator is matrix-based. If the operator is not matrix-based
	 * this inverse can not be used and false is returned
	 *
	 * \param[in]	A		linear matrix-based operator to invert
	 * \param[in]	u		linearization point
	 * \returns		bool	success flag
	 */
		virtual bool init(SmartPtr<ILinearOperator<Y,X> > A, const Y&u)
		{
		//	forget about u and forward request.
			return init(A);
		}

	///	initializes this inverse operator for a linear operator
	/**
	 * This method implements the ILinearOperatorInverse interface method.
	 * Basically, the request is forwarded to the matrix-based init method,
	 * if the the operator is matrix-based. If the operator is not matrix-based
	 * this inverse can not be used and false is returned
	 *
	 * \param[in]	A		linear matrix-based operator to invert
	 * \returns		bool	success flag
	 */
		virtual bool init(SmartPtr<ILinearOperator<Y,X> > A)
		{
		//	cast operator
			SmartPtr<MatrixOperator<M,Y,X> > op =
									A.template cast_dynamic<MatrixOperator<M,Y,X> >();

		//	check if correct types are present
			if(op.invalid())
				UG_THROW("IMatrixOperatorInverse::init:"
						" Passed operator is not matrix-based.");

		//	forward request
			return init(op);
		}
};

}
#endif /* __H__LIB_ALGEBRA__OPERATOR__INTERFACE__MATRIX_OPERATOR_INVERSE__ */
