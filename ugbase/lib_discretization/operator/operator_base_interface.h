/*
 * operator_base_interface.h
 *
 *  Created on: 30.09.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__OPERATOR__OPERATOR_BASE_INTERFACE__
#define __H__LIB_DISCRETIZATION__OPERATOR__OPERATOR_BASE_INTERFACE__

#include "operator_base_interface.h"

namespace ug{


class IFunctionBase
{
	public:
		virtual number two_norm() = 0;
};


///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
// (Nonlinear) Operator
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////

// describes a mapping N
class IOperatorBase
{
	public:
		// prepare Operator
		virtual bool prepare(IFunctionBase& dOut, IFunctionBase& uIn) = 0;

		// apply Operator, i.e. d := N(u);
		virtual bool apply(IFunctionBase& dOut, IFunctionBase& uIn) = 0;

		// destructor
		virtual ~IOperatorBase() {};
};

// describes an inverse mapping
class IOperatorInverseBase
{
	public:
		// This operator inverts the Operator N
		virtual bool init(IOperatorBase& N) = 0;

		// prepare Operator
		virtual bool prepare(IFunctionBase& uOut) = 0;

		// apply Operator, i.e. u = N^{-1}(0)
		virtual bool apply(IFunctionBase& uOut) = 0;

		// destructor
		virtual ~IOperatorInverseBase() {};
};


///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
// Linearized Operator
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////

// describes a linearized mapping X->Y
class ILinearizedOperatorBase
{
	public:
		// prepare J(u) for application of d = J(u)*c
		virtual bool prepare(IFunctionBase& dOut, IFunctionBase& uIn, IFunctionBase& cIn) = 0;

		// apply Operator, i.e. d = J(u)*c
		virtual bool apply(IFunctionBase& dOut, IFunctionBase& cIn) = 0;

		// apply Operator, i.e. d = d - J(u)*c;
		virtual bool apply_sub(IFunctionBase& dOut, IFunctionBase& cIn) = 0;

		// destructor
		virtual ~ILinearizedOperatorBase() {};
};

// Stepping operator to compute new correction
class ILinearizedIteratorOperatorBase
{
	public:
		// prepare for Linearized Operator J(u)
		virtual bool init(ILinearizedOperatorBase& J) = 0;

		// prepare B(u) for application of B(u)*d = c
		virtual bool prepare(IFunctionBase& cOut, IFunctionBase& uIn, IFunctionBase& dIn) = 0;

		/** apply
		 *
		 * This function computes a new correction c = B(u)*d by applying the Operator.
		 * The defect is updated, if updateDefect is true
		 *
		 * \param[in] 	d 				Defect
		 * \param[out]	c				Correction, c = B(u)*d
		 * \param[in]	updateDefect	if true, the defect is updated, d:= d - J(u)*c
		 */
		virtual bool apply(IFunctionBase& cOut, IFunctionBase& dIn, bool updateDefect) = 0;

		// destructor
		virtual ~ILinearizedIteratorOperatorBase() {};
};

class ILinearizedOperatorInverseBase
{
	public:
		// init for Linearized Operator A
		virtual bool init(ILinearizedOperatorBase& A) = 0;

		// prepare Operator
		virtual bool prepare(IFunctionBase& cOut, IFunctionBase& uIn, IFunctionBase& dIn) = 0;

		// solve c = J(u)^{-1}*d
		virtual bool apply(IFunctionBase& cOut, IFunctionBase& dInOut) = 0;

		// destructor
		virtual ~ILinearizedOperatorInverseBase() {};
};

///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
// Linear Operator
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////

class ILinearOperatorBase
{
	public:
		// init
		virtual bool init() = 0;

		// prepare Operator
		virtual bool prepare(IFunctionBase& fOut, IFunctionBase& uIn) = 0;

		// apply Operator, i.e. f = L*u
		virtual bool apply(IFunctionBase& fOut, IFunctionBase& uIn) = 0;

		// apply Operator, i.e. f = f - L*u;
		virtual bool apply_sub(IFunctionBase& fOut, IFunctionBase& uIn) = 0;

		// destructor
		virtual ~ILinearOperatorBase() {};
};

// stepping operator
class ILinearIteratorOperatorBase
{
	public:
		// prepare for Operator
		virtual bool init(ILinearOperatorBase& A) = 0;

		// prepare for correction and defect
		virtual bool prepare(IFunctionBase& cOut, IFunctionBase& dIn) = 0;

		// compute c = B*d
		virtual bool apply(IFunctionBase& cOut, IFunctionBase& dIn, bool updateDefect) = 0;

		// destructor
		virtual ~ILinearIteratorOperatorBase() {};

	public:
};

class ILinearOperatorInverseBase
{
	public:
		// init for Linear Operator A
		virtual bool init(ILinearOperatorBase& A) = 0;

		// prepare Operator
		virtual bool prepare(IFunctionBase& uOut, IFunctionBase& fIn) = 0;

		// Solve u = A^{-1} f
		virtual bool apply(IFunctionBase& uOut, IFunctionBase& fInOut) = 0;

		// destructor
		virtual ~ILinearOperatorInverseBase() {};
};


} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__OPERATOR__OPERATOR_BASE_INTERFACE__ */
