/*
 * linear_operator.h
 *
 *  Created on: 22.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__LINEAR_OPERATOR__LINEAR_OPERATOR__
#define __H__LIBDISCRETIZATION__LINEAR_OPERATOR__LINEAR_OPERATOR__

namespace ug{

// describes a mapping X->Y
template <typename X, typename Y>
class Operator
{
	public:
		// export types:

		// domain space
		typedef X domain_type;

		// range space
		typedef Y codomain_type;

		// domain function type
		typedef typename X::function_type domain_function_type;

		// codomain function type
		typedef typename Y::function_type codomain_function_type;

	public:
		// prepare Operator
		virtual bool prepare() = 0;

		// apply Operator, i.e. v = L(u);
		virtual bool apply(const domain_function_type& u, codomain_function_type& v) = 0;

		// destructor
		virtual ~Operator() {};
};


template <typename X, typename Y>
class LinearOperator :public Operator<X,Y>
{
	public:
		// export types:

		// domain space
		typedef X domain_type;

		// range space
		typedef Y codomain_type;

		// domain function type
		typedef typename X::function_type domain_function_type;

		// codomain function type
		typedef typename Y::function_type codomain_function_type;

	public:
		// prepare Operator
		virtual bool prepare() = 0;

		// apply Operator, i.e. v = L*u;
		virtual bool apply(const domain_function_type& u, codomain_function_type& v) = 0;

		// destructor
		virtual ~LinearOperator() {};
};


template <typename X, typename Y>
class DiscreteLinearOperator : public LinearOperator<X,Y>
{
	public:
		// export types:

		// domain space
		typedef X domain_type;

		// range space
		typedef Y codomain_type;

		// domain function type
		typedef typename X::function_type domain_function_type;

		// codomain function type
		typedef typename Y::function_type codomain_function_type;

	public:
		// prepare Operator
		virtual bool prepare() = 0;

		// apply Operator, i.e. v = L*u;
		virtual bool apply(const domain_function_type& u, codomain_function_type& v) = 0;

		// destructor
		virtual ~DiscreteLinearOperator() {};
};


template <typename ApproximationSpaceX, typename ApproximationSpaceY, typename TAlgebra>
class MatrixBasedLinearOperator : public DiscreteLinearOperator<ApproximationSpaceX, ApproximationSpaceY>
{
	public:
		// export types:

		// domain space
		typedef ApproximationSpaceX domain_type;

		// range space
		typedef ApproximationSpaceY codomain_type;

		// domain function type
		typedef typename ApproximationSpaceX::function_type domain_function_type;

		// codomain function type
		typedef typename ApproximationSpaceY::function_type codomain_function_type;

	public:
		// prepare the operator for application (e.g. compute an intern Matrix L)
		virtual bool prepare() = 0;

		// compute v = L*u (here, L is a Matrix)
		virtual bool apply(const domain_function_type& u, codomain_function_type& v)
		{
			const typename domain_function_type::vector_type& x = u.get_vector();
			typename codomain_function_type::vector_type& y = v.get_vector();

			return m_Matrix.apply(u,v);
		}

		// destructor
		virtual ~MatrixBasedLinearOperator() {};

	protected:
		// matrix type used
		typedef typename TAlgebra::matrix_type matrix_type;

	protected:
		// matrix storage
		matrix_type m_Matrix;
};


}

#include "interpolation_operator.h"
#include "transfer_operator.h"

#endif /* __H__LIBDISCRETIZATION__LINEAR_OPERATOR__LINEAR_OPERATOR__ */
