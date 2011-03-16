/*
 * operator_inverse_interface.h
 *
 *  Created on: 22.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_ALGEBRA__OPERATOR__OPERATOR_INVERSE_INTERFACE__
#define __H__LIB_ALGEBRA__OPERATOR__OPERATOR_INVERSE_INTERFACE__

#include "operator_interface.h"
#include "convergence_check.h"

namespace ug{

///////////////////////////////////////////////////////////////////////////////
// Inverse of a (Nonlinear-) Operator
///////////////////////////////////////////////////////////////////////////////

// \todo: The prepare methods seems obsolet, since it could be handled in the
//		  apply method. It depends on how expensive the prepare is, and how
//		  often we would apply the operator to an already apply one (??)
/// describes an inverse mapping X->Y
/**
 * This class is the base class for the inversion of the operator given in form
 * of the IOperator interface class. Given a operator N(u), the basic usage
 * of this class is to invert this operator, i.e. to compute
 *
 *     N(u) = 0, i.e. u := N^{-1}(0).
 *
 * This application has been split up into three steps:
 *
 * 1. init(N): This method initializes the inverse operator. The inverse operator
 * 			   is initialized the way that, its application will be the inverse
 * 			   application of the operator N passed in by this function. The
 * 			   prepare method can only be called, when this method has been
 * 			   called once.
 *
 * 2. prepare(u): This method prepares the function u, before it can be used in
 * 				  the apply method. Here, typically dirichlet values are set.
 *
 * 3. apply(u):	This method performs the inversion. Before this method is called
 * 				the init and prepare methods have to be called.
 *
 * This splitting has been made, since initialization and preparation may be
 * computationally expansive. Thus, the user of this class has the choice
 * when to call this initialization/preparation. E.g. when the operator is
 * applied several times on the same vectors, those have only to be prepared
 * once and the init of the operator is only needed once.
 */
template <typename X, typename Y>
class IOperatorInverse
{
	public:
	///	Domain space
		typedef X domain_function_type;

	/// Range space
		typedef Y codomain_function_type;

	public:
	/// initializes the operator for the inversion of the operator N: Y -> X
	/**
	 * This method sets the (nonlinear) operator that should be inverted by
	 * this inverse operator. In addition preparations can be made to
	 * facilitate the application of the inverse.
	 *
	 * \param[in]	N		operator that is to be inverted
	 * \returns 	bool	success flag
	 */
		virtual bool init(IOperator<Y,X>& N) = 0;

	/// prepares the function u for application
	/**
	 * This method prepares the function u before it can be used to find the
	 * solution of N(u) = 0. Typically, dirichlet values are set here.
	 *
	 * \param[in]	u		domain function
	 * \returns		bool	success flag
	 */
		virtual bool prepare(X& u) = 0;

	/// apply the operator, i.e. u = N^{-1}(0)
	/**
	 * This method inverts the operator N and returns the solution u = N^{-1}(0).
	 *
	 * \param[in,out]	u		domain function with solution at output
	 * \returns			bool	success flag
	 */
		virtual bool apply(X& u) = 0;

	/// virtual destructor
		virtual ~IOperatorInverse() {};
};


///////////////////////////////////////////////////////////////////////////////
// Inverse of a Linear Operator
///////////////////////////////////////////////////////////////////////////////

/// describes an inverse linear mapping X->Y
/**
 * This class is the base class for the inversion of linear soperator given in
 * form of the ILinearOperator interface class. Given a operator L, the basic
 * usage of this class is to invert this operator, i.e. to compute the solution
 * u of
 *
 * 		L*u = f     i.e. u := L^{-1} f
 *
 * This application has been split up into three steps:
 *
 * 1. init(N): This method initializes the inverse operator. The inverse operator
 * 			   is initialized the way that, its application will be the inverse
 * 			   application of the operator N passed in by this function. The
 * 			   prepare method can only be called, when this method has been
 * 			   called once.
 *
 * 2. prepare(u): This method prepares the function u, before it can be used in
 * 				  the apply method. Here, typically dirichlet values are set.
 *
 * 3. apply(u):	This method performs the inversion. Before this method is called
 * 				the init and prepare methods have to be called.
 *
 * This splitting has been made, since initialization and preparation may be
 * computationally expansive. Thus, the user of this class has the choice
 * when to call this initialization/preparation. E.g. when the operator is
 * applied several times on the same vectors, those have only to be prepared
 * once and the init of the operator is only needed once.
 */

template <typename X, typename Y>
class ILinearOperatorInverse
{
	public:
	///	Domain space
		typedef X domain_function_type;

	///	Range space
		typedef Y codomain_function_type;

	public:
	///	returns the name of the operator inverse
	/**
	 * This method returns the name of the inverse operator. This function is
	 * typically needed, when the inverse operator is used inside of another and
	 * some debug output should be printed
	 *
	 * \returns 	const char* 	name of inverse operator
	 */
		virtual const char* name() const = 0;

	/// initializes for the inverse for a linear operator
	/**
	 * This method passes the operator L that is inverted by this operator. In
	 * addition some preparation step can be made.
	 *
	 * \param[in]	L		linear operator to invert
	 * \returns		bool	success flag
	 */
		virtual bool init(ILinearOperator<Y,X>& L) = 0;

	/// initializes for the inverse for a linearized operator at linearization point u
	/**
	 * This method passes the linear operator J(u) that should be inverted by
	 * this operator. As second argument the linearization point is passed.
	 * This is needed e.g. for the geometric multigrid method, that inverts
	 * a linearized operator based on coarser grid operators, that have to be
	 * initialized based on the linearization point.
	 *
	 * \param[in]	J		linearized operator to invert
	 * \param[in]	u		linearization point
	 * \returns		bool	success flag
	 */
		virtual bool init(ILinearOperator<Y,X>& J, const Y& u) = 0;

	///	applies inverse operator, i.e. returns u = A^{-1} f
	/**
	 * This method applies the inverse operator, i.e. u = A^{-1} f. The
	 * domain function f remains unchanged.
	 * Note, that this method can always be implemented by creating a copy of
	 * f and calling apply_return_defect with this copy.
	 *
	 * \param[in]	f		right-hand side
	 * \param[out]	u		solution
	 * \returns		bool	success flag
	 */
		virtual bool apply(Y& u, const X& f) = 0;

	///	applies inverse operator, i.e. returns u = A^{-1} f and returns defect d := f - A*u
	/**
	 * This method applies the inverse operator, i.e. u = A^{-1} f. The
	 * domain function f is changed in the way, that the defect d := f - A*u
	 * is returned in the function. This is always useful, when the inverting
	 * algorithm can (or must) update the defect during computation (this is
	 * e.g. the case for the geometric multigrid method).
	 * Note, that this method can always be implemented by calling apply and
	 * then computing d := f - A*u.
	 *
	 * \param[in,out]	f		right-hand side
	 * \param[out]		u		solution
	 * \returns			bool	success flag
	 */
		virtual bool apply_return_defect(Y& u, X& f) = 0;

	///	set the convergence check
		virtual void set_convergence_check(IConvergenceCheck& convCheck) = 0;

	///	returns the convergence check
		virtual IConvergenceCheck* get_convergence_check() = 0;

	/// virtual destructor
		virtual ~ILinearOperatorInverse() {};
};


///////////////////////////////////////////////////////////////////////////////
// Inverse of a Matrix Linear Operator
///////////////////////////////////////////////////////////////////////////////

template <typename X, typename Y, typename M>
class IMatrixOperatorInverse
	: public virtual ILinearOperatorInverse<X,Y>
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


} // end namespace ug

#endif /* __H__LIB_ALGEBRA__OPERATOR__OPERATOR_INVERSE_INTERFACE__ */
