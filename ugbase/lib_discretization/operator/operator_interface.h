/*
 * operator_interface.h
 *
 *  Created on: 22.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__OPERATOR__OPERATOR_INTERFACE__
#define __H__LIB_DISCRETIZATION__OPERATOR__OPERATOR_INTERFACE__

#include "operator_base_interface.h"

namespace ug{

///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
// Operator
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////

// describes a mapping X->Y
template <typename X, typename Y>
class IOperator : public IOperatorBase
{
	public:
		// domain space
		typedef X domain_function_type;

		// range space
		typedef Y codomain_function_type;

	public:
		// init
		virtual bool init() = 0;

		// prepare Operator
		virtual bool prepare(Y& dOut, X& uIn) = 0;

		// apply Operator, i.e. d := L(u);
		virtual bool apply(Y& dOut, X& uIn) = 0;

		// destructor
		virtual ~IOperator() {};

	public:
		virtual bool prepare(IFunctionBase& dOut, IFunctionBase& uIn)
		{
		//	cast function to needed class types
			X* u = dynamic_cast<X*>(&uIn);
			Y* d = dynamic_cast<Y*>(&dOut);

		//	check if correct types are present
			if(u == NULL)
			{
				UG_LOG("ERROR in IOperator::prepare:"
						" Wrong type of function u detected.\n");
				return false;
			}
			if(d == NULL)
			{
				UG_LOG("ERROR in IOperator::prepare:"
						" Wrong type of function d detected.\n");
				return false;
			}

		//	types are correct. Forward to type aware function
			return prepare(*d, *u);
		}

		virtual bool apply(IFunctionBase& dOut, IFunctionBase& uIn)
		{
		//	cast function to needed class types
			X* u = dynamic_cast<X*>(&uIn);
			Y* d = dynamic_cast<Y*>(&dOut);

		//	check if correct types are present
			if(u == NULL)
			{
				UG_LOG("Type mismatch in IOperator::apply:"
						" Wrong type of function u detected.\n");
				return false;
			}
			if(d == NULL)
			{
				UG_LOG("Type mismatch in IOperator::apply:"
						" Wrong type of function d detected.\n");
				return false;
			}

		//	types are correct. Forward to type aware function
			return apply(*d, *u);
		}
};

// describes a mapping X->Y
template <typename X, typename Y>
class IOperatorInverse : public IOperatorInverseBase
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
		virtual bool apply(X& uOut) = 0;

		// destructor
		virtual ~IOperatorInverse() {};

	public:
		virtual bool init(IOperatorBase& N)
		{
		//	cast operator
			IOperator<Y,X>* op = dynamic_cast<IOperator<Y,X>*>(&N);

		//	check if correct types are present
			if(op == NULL)
			{
				UG_LOG("Type mismatch in IOperatorInverse::init:"
						" Wrong type of operator N detected.\n");
				return false;
			}

			return init(*op);
		}

		virtual bool prepare(IFunctionBase& uOut)
		{
		//	cast function to needed class types
			X* u = dynamic_cast<X*>(&uOut);

		//	check if correct types are present
			if(u == NULL)
			{
				UG_LOG("Type mismatch in IOperatorInverse::prepare:"
						" Wrong type of function u detected.\n");
				return false;
			}

		//	types are correct. Forward to type aware function
			return prepare(*u);
		}

		virtual bool apply(IFunctionBase& uOut)
		{
		//	cast function to needed class types
			X* u = dynamic_cast<X*>(&uOut);

		//	check if correct types are present
			if(u == NULL)
			{
				UG_LOG("Type mismatch in IOperatorInverse::apply:"
						" Wrong type of function u detected.\n");
				return false;
			}

		//	types are correct. Forward to type aware function
			return apply(*u);
		}
};


///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
// Linearized Operator
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////


// describes a mapping X->Y
template <typename X, typename Y>
class ILinearizedOperator : public ILinearizedOperatorBase
{
	public:
		// domain space
		typedef X domain_function_type;

		// range space
		typedef Y codomain_function_type;

	public:
		// init
		virtual bool init() = 0;

		// prepare J(u) for application of d = J(u)*c
		virtual bool prepare(Y& dOut, X& uIn, X& cIn) = 0;

		// apply Operator, i.e. d = J(u)*c
		virtual bool apply(Y& dOut, X& cIn) = 0;

		// apply Operator, i.e. d = d - J(u)*c;
		virtual bool apply_sub(Y& dOut, X& cIn) = 0;

		// destructor
		virtual ~ILinearizedOperator() {};

	public:
		virtual bool prepare(IFunctionBase& dOut, IFunctionBase& uIn, IFunctionBase& cIn)
		{
		//	cast function to needed class types
			X* c = dynamic_cast<X*>(&cIn);
			X* u = dynamic_cast<X*>(&uIn);
			Y* d = dynamic_cast<Y*>(&dOut);

		//	check if correct types are present
			if(c == NULL)
			{
				UG_LOG("Type mismatch in ILinearizedOperator::prepare:"
						" Wrong type of function c detected.\n");
				return false;
			}
			if(u == NULL)
			{
				UG_LOG("Type mismatch in ILinearizedOperator::prepare:"
						" Wrong type of function u detected.\n");
				return false;
			}
			if(d == NULL)
			{
				UG_LOG("Type mismatch in ILinearizedOperator::prepare:"
						" Wrong type of function d detected.\n");
				return false;
			}

		//	types are correct. Forward to type aware function
			return prepare(*d, *u, *c);
		}

		virtual bool apply(IFunctionBase& dOut, IFunctionBase& cIn)
		{
		//	cast function to needed class types
			X* c = dynamic_cast<X*>(&cIn);
			Y* d = dynamic_cast<Y*>(&dOut);

		//	check if correct types are present
			if(c == NULL)
			{
				UG_LOG("Type mismatch in ILinearizedOperator::apply:"
						" Wrong type of function c detected.\n");
				return false;
			}
			if(d == NULL)
			{
				UG_LOG("Type mismatch in ILinearizedOperator::apply:"
						" Wrong type of function d detected.\n");
				return false;
			}

		//	types are correct. Forward to type aware function
			return apply(*d, *c);
		}

		virtual bool apply_sub(IFunctionBase& dOut, IFunctionBase& cIn)
		{
		//	cast function to needed class types
			X* c = dynamic_cast<X*>(&cIn);
			Y* d = dynamic_cast<Y*>(&dOut);

		//	check if correct types are present
			if(c == NULL)
			{
				UG_LOG("Type mismatch in ILinearizedOperator::apply_sub:"
						" Wrong type of function c detected.\n");
				return false;
			}
			if(d == NULL)
			{
				UG_LOG("Type mismatch in ILinearizedOperator::apply_sub:"
						" Wrong type of function d detected.\n");
				return false;
			}

		//	types are correct. Forward to type aware function
			return apply_sub(*d, *c);
		}
};

/* This Operator type behaves different on application. It not only computes c = B(u)*d, but also changes d. */
/* It is used in iterative schemes. */
template <typename X, typename Y>
class ILinearizedIteratorOperator : public ILinearizedIteratorOperatorBase
{
	public:
		// domain space
		typedef X domain_function_type;

		// range space
		typedef Y codomain_function_type;

	public:
		// prepare for Linearized Operator J(u)
		virtual bool init(ILinearizedOperator<Y, X>& J) = 0;

		// prepare B(u) for application of B(u)*d = c
		virtual bool prepare(Y& cOut, X& uIn, X& dIn) = 0;

		/** apply
		 *
		 * This function computes a new correction c = B(u)*d by
		 * applying the Operator. The defect is updated, if
		 * updateDefect is true.
		 *
		 * \param[in] 	d 				Defect
		 * \param[out]	c				Correction, c = B(u)*d
		 * \param[in]	updateDefect	if true, the defect is updated, d:= d - J(u)*c
		 */
		virtual bool apply(Y& cOut, X& dIn, bool updateDefect) = 0;

		// clone
		virtual ILinearizedIteratorOperator<X,Y>* clone() = 0;

		// destructor
		virtual ~ILinearizedIteratorOperator() {};

	public:
		virtual bool init(ILinearizedOperatorBase& N)
		{
		//	cast operator
			ILinearizedOperator<Y,X>* op = dynamic_cast<ILinearizedOperator<Y,X>*>(&N);

		//	check if correct types are present
			if(op == NULL)
			{
				UG_LOG("Type mismatch in ILinearizedIteratorOperator::init:"
						" Wrong type of operator N detected.\n");
				return false;
			}

			return init(*op);
		}

		virtual bool prepare(IFunctionBase& cOut, IFunctionBase& uIn, IFunctionBase& dIn)
		{
		//	cast function to needed class types
			X* u = dynamic_cast<X*>(&uIn);
			X* d = dynamic_cast<X*>(&dIn);
			Y* c = dynamic_cast<Y*>(&cOut);

		//	check if correct types are present
			if(u == NULL)
			{
				UG_LOG("Type mismatch in ILinearizedIteratorOperator::prepare:"
						" Wrong type of function u detected.\n");
				return false;
			}
			if(d == NULL)
			{
				UG_LOG("Type mismatch in ILinearizedIteratorOperator::prepare:"
						" Wrong type of function d detected.\n");
				return false;
			}
			if(c == NULL)
			{
				UG_LOG("Type mismatch in ILinearizedIteratorOperator::prepare:"
						" Wrong type of function c detected.\n");
				return false;
			}

		//	types are correct. Forward to type aware function
			return prepare(*c, *u, *d);
		}

		virtual bool apply(IFunctionBase& cOut, IFunctionBase& dIn, bool updateDefect)
		{
		//	cast function to needed class types
			X* d = dynamic_cast<X*>(&dIn);
			Y* c = dynamic_cast<Y*>(&cOut);

		//	check if correct types are present
			if(d == NULL)
			{
				UG_LOG("Type mismatch in ILinearizedIteratorOperator::apply:"
						" Wrong type of function d detected.\n");
				return false;
			}
			if(c == NULL)
			{
				UG_LOG("Type mismatch in ILinearizedIteratorOperator::apply:"
						" Wrong type of function c detected.\n");
				return false;
			}

		//	types are correct. Forward to type aware function
			return apply(*c, *d, updateDefect);
		}
};

template <typename X, typename Y>
class ILinearizedOperatorInverse : public ILinearizedOperatorInverseBase
{
	public:
		// domain space
		typedef X domain_function_type;

		// range space
		typedef Y codomain_function_type;

	public:
		// init for Linearized Operator A
		virtual bool init(ILinearizedOperator<Y,X>& A) = 0;

		// prepare Operator
		virtual bool prepare(Y& cOut, X& uIn, X& dIn) = 0;

		// solve c = J(u)^{-1}*d
		// The updated defect is returned in d
		virtual bool apply(Y& cOut, X& dInOut) = 0;

		// destructor
		virtual ~ILinearizedOperatorInverse() {};

	public:
		virtual bool init(ILinearizedOperatorBase& N)
		{
		//	cast operator
			ILinearizedOperator<Y,X>* op = dynamic_cast<ILinearizedOperator<Y,X>*>(&N);

		//	check if correct types are present
			if(op == NULL)
			{
				UG_LOG("Type mismatch in ILinearizedOperatorInverse::init:"
						" Wrong type of operator N detected.\n");
				return false;
			}

			return init(*op);
		}

		virtual bool prepare(IFunctionBase& cOut, IFunctionBase& uIn, IFunctionBase& dIn)
		{
		//	cast function to needed class types
			X* u = dynamic_cast<X*>(&uIn);
			X* d = dynamic_cast<X*>(&dIn);
			Y* c = dynamic_cast<Y*>(&cOut);

		//	check if correct types are present
			if(u == NULL)
			{
				UG_LOG("Type mismatch in ILinearizedOperatorInverse::prepare:"
						" Wrong type of function u detected.\n");
				return false;
			}
			if(d == NULL)
			{
				UG_LOG("Type mismatch in ILinearizedOperatorInverse::prepare:"
						" Wrong type of function d detected.\n");
				return false;
			}
			if(c == NULL)
			{
				UG_LOG("Type mismatch in ILinearizedOperatorInverse::prepare:"
						" Wrong type of function c detected.\n");
				return false;
			}

		//	types are correct. Forward to type aware function
			return prepare(*c, *u, *d);
		}

		virtual bool apply(IFunctionBase& cOut, IFunctionBase& dInOut)
		{
		//	cast function to needed class types
			X* d = dynamic_cast<X*>(&dInOut);
			Y* c = dynamic_cast<Y*>(&cOut);

		//	check if correct types are present
			if(d == NULL)
			{
				UG_LOG("Type mismatch in ILinearizedOperatorInverse::apply:"
						" Wrong type of function d detected.\n");
				return false;
			}
			if(c == NULL)
			{
				UG_LOG("Type mismatch in ILinearizedOperatorInverse::apply:"
						" Wrong type of function c detected.\n");
				return false;
			}

		//	types are correct. Forward to type aware function
			return apply(*c, *d);
		}

};

///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
// Linear Operator
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////

template <typename X, typename Y>
class ILinearOperator : public ILinearOperatorBase, public ILinearizedOperator<X,Y>
{
	public:
		// domain space
		typedef X domain_function_type;

		// range space
		typedef Y codomain_function_type;

	public:
		// init
		virtual bool init() = 0;

		// prepare Operator
		virtual bool prepare(Y& fOut, X& uIn) = 0;

		// implement the interface for Linearized Operator
		virtual bool prepare(Y& dOut, X& cIn, X& uIn) {return prepare(dOut,cIn);};

		// apply Operator, i.e. f = L*u; (or d = J*c in case of Linearized Operator, i.e. u = c, f = d)
		virtual bool apply(Y& fOut, X& uIn) = 0;

		// apply Operator, i.e. f = f - L*u;
		virtual bool apply_sub(Y& fOut, X& uIn) = 0;

		// destructor
		virtual ~ILinearOperator() {};

	public:
		virtual bool prepare(IFunctionBase& fOut, IFunctionBase& uIn)
		{
		//	cast function to needed class types
			X* u = dynamic_cast<X*>(&uIn);
			Y* f = dynamic_cast<Y*>(&fOut);

		//	check if correct types are present
			if(u == NULL)
			{
				UG_LOG("Type mismatch in ILinearOperator::prepare:"
						" Wrong type of function u detected.\n");
				return false;
			}
			if(f == NULL)
			{
				UG_LOG("Type mismatch in ILinearOperator::prepare:"
						" Wrong type of function f detected.\n");
				return false;
			}

		//	types are correct. Forward to type aware function
			return prepare(*f, *u);
		}

		virtual bool apply(IFunctionBase& fOut, IFunctionBase& uIn)
		{
		//	cast function to needed class types
			X* u = dynamic_cast<X*>(&uIn);
			Y* f = dynamic_cast<Y*>(&fOut);

		//	check if correct types are present
			if(u == NULL)
			{
				UG_LOG("Type mismatch in ILinearOperator::apply:"
						" Wrong type of function u detected.\n");
				return false;
			}
			if(f == NULL)
			{
				UG_LOG("Type mismatch in ILinearOperator::apply:"
						" Wrong type of function f detected.\n");
				return false;
			}

		//	types are correct. Forward to type aware function
			return apply(*f, *u);
		}

		virtual bool apply_sub(IFunctionBase& fOut, IFunctionBase& uIn)
		{
		//	cast function to needed class types
			X* u = dynamic_cast<X*>(&uIn);
			Y* f = dynamic_cast<Y*>(&fOut);

		//	check if correct types are present
			if(u == NULL)
			{
				UG_LOG("Type mismatch in ILinearOperator::apply_sub:"
						" Wrong type of function u detected.\n");
				return false;
			}
			if(f == NULL)
			{
				UG_LOG("Type mismatch in ILinearOperator::apply_sub:"
						" Wrong type of function f detected.\n");
				return false;
			}

		//	types are correct. Forward to type aware function
			return apply_sub(*f, *u);
		}

};

/* This Operator type behaves different on application. It not only computes c = B*d, but also changes d. */
/* It is used in iterative schemes. */
template <typename X, typename Y>
class ILinearIteratorOperator : public ILinearIteratorOperatorBase,
								public ILinearizedIteratorOperator<X,Y>
{
	public:
		// domain space
		typedef X domain_function_type;

		// range space
		typedef Y codomain_function_type;

	public:
		// prepare for Operator
		virtual bool init(ILinearOperator<Y,X>& A) = 0;

		// prepare for correction and defect
		virtual bool prepare(Y& cOut, X& dIn) = 0;

		// Implement Interface for Linearized Operator
		virtual bool prepare(Y& cOut, X& uIn, X& dIn){return prepare(cOut, dIn);}

		/** apply
		 *
		 * This function computes a new correction c = B*d by applying the Operator.
		 * The defect is updated, if updateDefect is true
		 *
		 * \param[in] 	d 				Defect
		 * \param[out]	c				Correction, c = B*d
		 * \param[in]	updateDefect	if true, the defect is updated, d:= d - A*c
		 */
		virtual bool apply(Y& cOut, X& dIn, bool updateDefect) = 0;

		// destructor
		virtual ~ILinearIteratorOperator() {};

	public:
		virtual bool init(ILinearOperatorBase& N)
		{
		//	cast operator
			ILinearOperator<Y,X>* op = dynamic_cast<ILinearOperator<Y,X>*>(&N);

		//	check if correct types are present
			if(op == NULL)
			{
				UG_LOG("Type mismatch in ILinearIteratorOperator::init:"
						" Wrong type of operator N detected.\n");
				return false;
			}

			return init(*op);
		}

		virtual bool prepare(IFunctionBase& cOut, IFunctionBase& dIn)
		{
		//	cast function to needed class types
			X* c = dynamic_cast<X*>(&cOut);
			Y* d = dynamic_cast<Y*>(&dIn);

		//	check if correct types are present
			if(c == NULL)
			{
				UG_LOG("Type mismatch in ILinearIteratorOperator::apply:"
						" Wrong type of function c detected.\n");
				return false;
			}
			if(d == NULL)
			{
				UG_LOG("Type mismatch in ILinearIteratorOperator::apply:"
						" Wrong type of function d detected.\n");
				return false;
			}

		//	types are correct. Forward to type aware function
			return prepare(*c, *d);
		}

		virtual bool apply(IFunctionBase& cOut, IFunctionBase& dIn, bool updateDefect)
		{
		//	cast function to needed class types
			X* c = dynamic_cast<X*>(&cOut);
			Y* d = dynamic_cast<Y*>(&dIn);

		//	check if correct types are present
			if(c == NULL)
			{
				UG_LOG("Type mismatch in ILinearIteratorOperator::apply:"
						" Wrong type of function c detected.\n");
				return false;
			}
			if(d == NULL)
			{
				UG_LOG("Type mismatch in ILinearIteratorOperator::apply:"
						" Wrong type of function d detected.\n");
				return false;
			}

		//	types are correct. Forward to type aware function
			return apply(*c, *d, updateDefect);
		}

};

template <typename X, typename Y>
class ILinearOperatorInverse : 	public ILinearOperatorInverseBase,
								public ILinearizedOperatorInverse<X,Y>
{
	public:
		// domain space
		typedef X domain_function_type;

		// range space
		typedef Y codomain_function_type;

	public:
		// init for Linearized Operator A
		virtual bool init(ILinearOperator<Y,X>& A) = 0;

		// prepare Operator
		virtual bool prepare(Y& uOut, X& fIn) = 0;

		// prepare Operator
		virtual bool prepare(Y& cOut, X& uIn, X& dIn) {return prepare(cOut,dIn);}

		// Solve A*u = f, such that u = A^{-1} f
		// This is done by iterating: u := u + B(f - A*u)
		// In f the last defect f := f - A*u is returned
		virtual bool apply(Y& uOut, X& fInOut) = 0;

		// destructor
		virtual ~ILinearOperatorInverse() {};

	public:
		virtual bool init(ILinearOperatorBase& N)
		{
		//	cast operator
			ILinearOperator<Y,X>* op = dynamic_cast<ILinearOperator<Y,X>*>(&N);

		//	check if correct types are present
			if(op == NULL)
			{
				UG_LOG("Type mismatch in ILinearOperatorInverse::init:"
						" Wrong type of operator N detected.\n");
				return false;
			}

			return init(*op);
		}

		virtual bool prepare(IFunctionBase& uOut, IFunctionBase& fIn)
		{
		//	cast function to needed class types
			X* u = dynamic_cast<X*>(&uOut);
			Y* f = dynamic_cast<Y*>(&fIn);

		//	check if correct types are present
			if(u == NULL)
			{
				UG_LOG("Type mismatch in ILinearOperatorInverse::prepare:"
						" Wrong type of function u detected.\n");
				return false;
			}
			if(f == NULL)
			{
				UG_LOG("Type mismatch in ILinearOperatorInverse::prepare:"
						" Wrong type of function f detected.\n");
				return false;
			}

		//	types are correct. Forward to type aware function
			return prepare(*u, *f);
		}

		virtual bool apply(IFunctionBase& uOut, IFunctionBase& fInOut)
		{
		//	cast function to needed class types
			X* u = dynamic_cast<X*>(&uOut);
			Y* f = dynamic_cast<Y*>(&fInOut);

		//	check if correct types are present
			if(u == NULL)
			{
				UG_LOG("Type mismatch in ILinearOperatorInverse::apply:"
						" Wrong type of function u detected.\n");
				return false;
			}
			if(f == NULL)
			{
				UG_LOG("Type mismatch in ILinearOperatorInverse::apply:"
						" Wrong type of function f detected.\n");
				return false;
			}

		//	types are correct. Forward to type aware function
			return apply(*u, *f);
		}
};


} // end namespace ug

#endif
