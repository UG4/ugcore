/*
 * operator.h
 *
 *  Created on: 22.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__OPERATOR__OPERATOR__
#define __H__LIBDISCRETIZATION__OPERATOR__OPERATOR__

#include <iomanip>
#include <cmath>

#include "common/common.h"
#include "lib_discretization/assemble.h"
#include "lib_discretization/io/vtkoutput.h"

namespace ug{

///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
// Operator
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////

// describes a mapping X->Y
template <typename X, typename Y>
class IOperator
{
	public:
		// export types:

		// domain space
		typedef X domain_function_type;

		// range space
		typedef Y codomain_function_type;

	public:
		// init
		virtual bool init() = 0;

		// prepare Operator
		virtual bool prepare(domain_function_type& u, codomain_function_type& d) = 0;

		// apply Operator, i.e. f := L(u);
		virtual bool apply(domain_function_type& u, codomain_function_type& d) = 0;

		// destructor
		virtual ~IOperator() {};
};

// describes a mapping X->Y
template <typename X, typename Y>
class IOperatorInverse
{
	public:
		// export types:

		// domain space
		typedef X domain_function_type;

		// range space
		typedef Y codomain_function_type;

	public:
		// init: This operator inverts the Operator N: Y -> X
		virtual bool init(IOperator<Y,X>& N) = 0;

		// prepare Operator
		virtual bool prepare(codomain_function_type& u) = 0;

		// apply Operator, i.e. N^{-1}(0) = u
		virtual bool apply(codomain_function_type& u) = 0;

		// destructor
		virtual ~IOperatorInverse() {};
};


///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
// Linearized Operator
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////


// describes a mapping X->Y
template <typename X, typename Y>
class ILinearizedOperator
{
	public:
		// export types:

		// domain space
		typedef X domain_function_type;

		// range space
		typedef Y codomain_function_type;

	public:
		// init
		virtual bool init() = 0;

		// prepare J(u) for application of J(u)*c = d
		virtual bool prepare(domain_function_type& u, domain_function_type& c, codomain_function_type& d) = 0;

		// apply Operator, i.e. d = J(u)*c
		virtual bool apply(domain_function_type& c, codomain_function_type& d) = 0;

		// apply Operator, i.e. f = f - L*u;
		virtual bool apply_sub(domain_function_type& u, codomain_function_type& f) = 0;

		// destructor
		virtual ~ILinearizedOperator() {};
};

/* This Operator type behaves different on application. It not only computes c = B(u)*d, but also changes d. */
/* It is used in iterative schemes. */
template <typename X, typename Y>
class ILinearizedIteratorOperator
{
	public:
		// export types:

		// domain space
		typedef X domain_function_type;

		// range space
		typedef Y codomain_function_type;

	public:
		// prepare for Linearized Operator J(u)
		virtual bool init(ILinearizedOperator<Y, X>& J) = 0;

		// prepare B(u) for application of B(u)*d = c
		virtual bool prepare(domain_function_type& u, domain_function_type& d, codomain_function_type& c) = 0;

		// compute new correction c = B(u)*d
		//    AND
		// update defect: d := d - J(u)*c
		virtual bool apply(domain_function_type& d, codomain_function_type& c) = 0;

		// clone
		virtual ILinearizedIteratorOperator<X,Y>* clone() = 0;

		// destructor
		virtual ~ILinearizedIteratorOperator() {};
};

template <typename X, typename Y>
class ILinearizedOperatorInverse
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
		virtual bool prepare(codomain_function_type& u, domain_function_type& d_nl, codomain_function_type& c_nl) = 0;

		// Solve J(u)*c_nl = d_nl, such that c_nl = J(u)^{-1} d_nl
		// This is done by iterating: c_nl := c_nl + B(u)(d_nl - J(u)*c_nl)
		// In d_nl the last defect d := d_nl - J(u)*c_nl is returned
		// In the following:
		// c_nl, d_nl refer to the non-linear defect and correction as e.g. in J(u) * c_nl = d_nl as it appears in Newton scheme
		// c, d are the correction and defect for solving that linear equation iteratively.
		virtual bool apply(domain_function_type& d_nl, codomain_function_type& c_nl) = 0;

		// destructor
		virtual ~ILinearizedOperatorInverse() {};
};

///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
// Linear Operator
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////

template <typename X, typename Y>
class ILinearOperator : public ILinearizedOperator<X,Y>
{
	public:
		// export types:

		// domain space
		typedef X domain_function_type;

		// range space
		typedef Y codomain_function_type;

	public:
		// init
		virtual bool init() = 0;

		// prepare Operator
		virtual bool prepare(domain_function_type& u, codomain_function_type& f) = 0;

		// implement the interface for Linearized Operator
		virtual bool prepare(domain_function_type& u, domain_function_type& c, codomain_function_type& d) {return prepare(c,d);};

		// apply Operator, i.e. f = L*u; (or d = J*c in case of Linearized Operator, i.e. u = c, f = d)
		virtual bool apply(domain_function_type& u, codomain_function_type& f) = 0;

		// apply Operator, i.e. f = f - L*u;
		virtual bool apply_sub(domain_function_type& u, codomain_function_type& f) = 0;

		// destructor
		virtual ~ILinearOperator() {};
};

/* This Operator type behaves different on application. It not only computes c = B*d, but also changes d. */
/* It is used in iterative schemes. */
template <typename X, typename Y>
class ILinearIteratorOperator : public ILinearizedIteratorOperator<X,Y>
{
	public:
		// export types:

		// domain space
		typedef X domain_function_type;

		// range space
		typedef Y codomain_function_type;

	public:
		// prepare for Operator
		virtual bool init(ILinearizedOperator<Y,X>& A) = 0;

		// prepare for correction and defect
		virtual bool prepare(domain_function_type& d, codomain_function_type& c) = 0;

		// Implement Interface for Linearized Operator
		virtual bool prepare(domain_function_type& u, domain_function_type& d, codomain_function_type& c){return prepare(d,c);}

		// compute new correction c = B*d
		//    AND
		// update defect: d := d - A*c
		virtual bool apply(domain_function_type& d, codomain_function_type& c) = 0;

		// destructor
		virtual ~ILinearIteratorOperator() {};
};

}

#include "linear_operator/interpolation_operator.h"
#include "linear_operator/projection_operator.h"
#include "linear_operator/transfer_operator.h"
#include "linear_operator/assembled_linear_operator.h"
#include "linear_operator/multi_grid_solver/mg_solver.h"
#include "linear_operator/linear_solver.h"
#include "linear_operator/cg_solver.h"

#ifdef LAPACK_AVAILABLE
#ifdef BLAS_AVAILABLE
#include "linear_operator/lapack_lu_operator.h" 
#endif
#endif

#include "non_linear_operator/assembled_non_linear_operator.h"
#include "non_linear_operator/newton_solver/newton.h"

#endif /* __H__LIBDISCRETIZATION__OPERATOR__OPERATOR__ */
