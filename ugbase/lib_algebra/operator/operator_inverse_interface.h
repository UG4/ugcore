/*
 * operator_inverse_interface.h
 *
 *  Created on: 22.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_ALGEBRA__OPERATOR__OPERATOR_INVERSE_INTERFACE__
#define __H__LIB_ALGEBRA__OPERATOR__OPERATOR_INVERSE_INTERFACE__

#include "operator_interface.h"

namespace ug{

template <typename X, typename Y>
class ILinearOperatorInverse
{
	public:
	// 	Domain space
		typedef X domain_function_type;

	// 	Range space
		typedef Y codomain_function_type;

	public:
	// 	Init for Linear Operator L
		virtual bool init(ILinearOperator<Y,X>& L) = 0;

	// 	Init for Linear Operator J and Linearization point (current solution)
		virtual bool init(ILinearOperator<Y,X>& J, const Y& u) = 0;

	// 	Solve A*u = f, such that u = A^{-1} f
		virtual bool apply(Y& u, const X& f) = 0;

	// 	Solve A*u = f, such that u = A^{-1} f
	// 	This is done by iterating: u := u + B(f - A*u)
	// 	In f the last defect f := f - A*u is returned
		virtual bool apply_return_defect(Y& u, X& f) = 0;

	// 	Destructor
		virtual ~ILinearOperatorInverse() {};
};



template <typename X, typename Y, typename M>
class IMatrixOperatorInverse : public virtual ILinearOperatorInverse<X,Y>
{
	public:
	//	Domain space
		typedef X domain_function_type;

	// 	Range space
		typedef Y codomain_function_type;

	// 	Matrix type
		typedef M matrix_type;

	public:
	// 	Init for Operator A
		virtual bool init(IMatrixOperator<Y,X,M>& A) = 0;

	// 	Apply solver, i.e. return u = A^{-1} * f
		virtual bool apply(Y& u, const X& f) = 0;

	// 	Solve A*u = f, such that u = A^{-1} f
	// 	This is done by iterating: u := u + B(f - A*u)
	// 	In f the last defect f := f - A*u is returned
		virtual bool apply_return_defect(Y& u, X& f) = 0;

	// 	Destructor
		virtual ~IMatrixOperatorInverse() {};

	public:
	//	Implement functions of LinearOperator
		virtual bool init(ILinearOperator<Y,X>& A, const Y&u)
		{
			return init(A);
		}

		virtual bool init(ILinearOperator<Y,X>& A)
		{
		//	cast operator
			IMatrixOperator<Y,X,M>* op = dynamic_cast<IMatrixOperator<Y,X,M>*>(&A);

		//	check if correct types are present
			if(op == NULL)
			{
				UG_LOG("Type mismatch in IMatrixOperatorInverse::init:"
						" Wrong type of operator A detected.\n");
				return false;
			}

			return init(*op);
		}
};


// describes a mapping X->Y
template <typename X, typename Y>
class IOperatorInverse
{
	public:
		// domain space
		typedef X domain_function_type;

		// range space
		typedef Y codomain_function_type;

	public:
		// This operator inverts the Operator N: Y -> X
		virtual bool init(IOperator<Y,X>& N) = 0;

		// prepare Operator
		virtual bool prepare(X& uOut) = 0;

		// apply Operator, i.e. u = N^{-1}(0)
		virtual bool apply(X& u) = 0;

		// destructor
		virtual ~IOperatorInverse() {};
};



} // end namespace ug

#endif /* __H__LIB_ALGEBRA__OPERATOR__OPERATOR_INVERSE_INTERFACE__ */
