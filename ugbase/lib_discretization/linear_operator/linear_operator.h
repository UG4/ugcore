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
		// Init Operator
		virtual bool init() = 0;

		// prepare Operator
		virtual bool prepare() = 0;

		// apply Operator, i.e. v = L(u);
		virtual bool apply(const domain_function_type& u, codomain_function_type v) = 0;
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
		// Init Operator
		virtual bool init() = 0;

		// prepare Operator
		virtual bool prepare() = 0;

		// apply Operator, i.e. v = L*u;
		virtual bool apply(const domain_function_type& u, codomain_function_type& v) = 0;
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
		// Init Operator
		virtual bool init() = 0;

		// prepare Operator
		virtual bool prepare() = 0;

		// apply Operator, i.e. v = L*u;
		virtual bool apply(const domain_function_type& u, codomain_function_type& v) = 0;

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
		// initialze the operator (e.g. reserve memory for the Matrix L)
		virtual bool init() = 0;

		// prepare the operator for application (e.g. compute an intern Matrix L)
		virtual bool prepare() = 0;

		// compute v = L*u (here, L is a Matrix)
		virtual bool apply(const domain_function_type& u, codomain_function_type& v)
		{
			const domain_function_type::vector_type& x = u.get_vector();
			codomain_function::vector_type& y = v.get_vector();

			return m_Matrix.apply(u,v);
		}

	protected:
		// matrix type used
		typedef typename TAlgebra::matrix_type matrix_type;

	protected:
		// matrix storage
		matrix_type m_Matrix;
};


}


#endif /* __H__LIBDISCRETIZATION__LINEAR_OPERATOR__LINEAR_OPERATOR__ */
