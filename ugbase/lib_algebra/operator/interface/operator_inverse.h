/*
 * operator_inverse.h
 *
 *  Created on: 22.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_ALGEBRA__OPERATOR__INTERFACE__OPERATOR_INVERSE__
#define __H__LIB_ALGEBRA__OPERATOR__INTERFACE__OPERATOR_INVERSE__

#include "common/util/smart_pointer.h"
#include "operator.h"
#include <string>

namespace ug{

///////////////////////////////////////////////////////////////////////////////
// Inverse of a (Nonlinear-) Operator
///////////////////////////////////////////////////////////////////////////////

// \todo: The prepare methods seems obsolet, since it could be handled in the
//		  apply method. It depends on how expensive the prepare is, and how
//		  often we would apply the operator to an already apply one (?)
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
template <typename X, typename Y = X>
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
		virtual bool init(SmartPtr<IOperator<Y,X> > N) = 0;

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

	///	returns information about configuration parameters
	/**
	 * this should return necessary information about parameters and possibly
	 * calling config_string of subcomponents.
	 *
	 * \returns std::string	necessary information about configuration parameters
	 */
		virtual std::string config_string() const = 0;
};




} // end namespace ug

#endif /* __H__LIB_ALGEBRA__OPERATOR__INTERFACE__OPERATOR_INVERSE__ */
