/*
 * operator_interface.h
 *
 *  Created on: 22.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_ALGEBRA__OPERATOR__OPERATOR_INTERFACE__
#define __H__LIB_ALGEBRA__OPERATOR__OPERATOR_INTERFACE__

#include "operator_base_interface.h"

namespace ug{

///////////////////////////////////////////////////////////
// Operator
///////////////////////////////////////////////////////////

// describes a mapping X->Y
template <typename X, typename Y>
class IOperator
{
	public:
	// 	Domain space
		typedef X domain_function_type;

	// 	Range space
		typedef Y codomain_function_type;

	public:
	// 	Init Operator
		virtual bool init() = 0;

	// 	Prepare out function
		virtual bool prepare(Y& d, X& u) = 0;

	// 	Apply Operator, i.e. d := N(u);
		virtual bool apply(Y& d, const X& u) = 0;

	// 	Destructor
		virtual ~IOperator() {};
};

///////////////////////////////////////////////////////////
// Linear Operator
///////////////////////////////////////////////////////////

// describes a mapping X->Y
template <typename X, typename Y>
class ILinearOperator
{
	public:
	// 	Domain space
		typedef X domain_function_type;

	// 	Range space
		typedef Y codomain_function_type;

	public:
	// 	Init Operator J(u)
		virtual bool init(const X& u) = 0;

	// 	Init Operator L
		virtual bool init() = 0;

	// 	Apply Operator f = L*u (e.g. d = J(u)*c in iterative scheme)
		virtual bool apply(Y& f, const X& u) = 0;

	// 	Apply Operator, i.e. f = f - L*u;
		virtual bool apply_sub(Y& f, const X& u) = 0;

	//  apply_axpy(dest, alpha, v, beta, w) dest = alpha*v + beta*A*w
	//	energy_prod(v, w) dest = v^T(Aw)

	// 	Destructor
		virtual ~ILinearOperator() {};
};

///////////////////////////////////////////////////////////
// Matrix Operator
///////////////////////////////////////////////////////////

template <typename X, typename Y, typename M>
class IMatrixOperator :	public virtual ILinearOperator<X,Y>
{
	public:
	// 	Domain space
		typedef X domain_function_type;

	// 	Range space
		typedef Y codomain_function_type;

	// 	Matrix type
		typedef M matrix_type;

	public:
	// 	Access to matrix
		virtual M& get_matrix() = 0;
};

template <typename X, typename Y, typename M>
class PureMatrixOperator :	public virtual IMatrixOperator<X,Y,M>
{
	public:
	// 	Domain space
		typedef X domain_function_type;

	// 	Range space
		typedef Y codomain_function_type;

	// 	Matrix type
		typedef M matrix_type;

	public:
	// 	Init Operator J(u)
		virtual bool init(const X& u) {return true;}

	// 	Init Operator L
		virtual bool init() {return true;}

	// 	Apply Operator f = L*u (e.g. d = J(u)*c in iterative scheme)
		virtual bool apply(Y& f, const X& u) {return m_Matrix.apply(f,u);}

	// 	Apply Operator, i.e. f = f - L*u;
		virtual bool apply_sub(Y& f, const X& u) {return m_Matrix.matmul_minus(f,u);}

	// 	Access to matrix
		virtual M& get_matrix() {return m_Matrix;};

	protected:
	//	memory
		M m_Matrix;
};

///////////////////////////////////////////////////////////
// Prolongation Operator
///////////////////////////////////////////////////////////

template <typename X, typename Y>
class IProlongationOperator :	public virtual ILinearOperator<X,Y>
{
	public:
	// 	Domain space
		typedef X domain_function_type;

	// 	Range space
		typedef Y codomain_function_type;

	public:
	// 	Apply Transposed Operator u = L^T*f
		virtual bool apply_transposed(X& u, const Y& f) = 0;

	// 	Set Levels for Prolongation coarse -> fine
		virtual bool set_levels(size_t coarseLevel, size_t fineLevel) = 0;

	//	Clone
		virtual IProlongationOperator<X,Y>* clone() = 0;
};

///////////////////////////////////////////////////////////
// Projection Operator
///////////////////////////////////////////////////////////

template <typename X, typename Y>
class IProjectionOperator :	public virtual ILinearOperator<X,Y>
{
	public:
	// 	Domain space
		typedef X domain_function_type;

	// 	Range space
		typedef Y codomain_function_type;

	public:
	// 	Apply Transposed Operator u = L^T*f
		virtual bool apply_transposed(X& u, const Y& f) = 0;

	// 	Set Levels for Prolongation coarse -> fine
		virtual bool set_levels(size_t coarseLevel, size_t fineLevel) = 0;

	//	Clone
		virtual IProjectionOperator<X,Y>* clone() = 0;
};

} // end namespace ug

#endif
