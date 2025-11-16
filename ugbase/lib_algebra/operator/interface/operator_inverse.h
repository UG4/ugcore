/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
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
		using domain_function_type = X;

	/// Range space
		using codomain_function_type = Y;

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
