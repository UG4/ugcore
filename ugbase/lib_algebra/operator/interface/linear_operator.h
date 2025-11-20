/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__LIB_ALGEBRA__OPERATOR__INTERFACE__LINEAR_OPERATOR__
#define __H__LIB_ALGEBRA__OPERATOR__INTERFACE__LINEAR_OPERATOR__

#include "operator.h"

namespace ug{

///////////////////////////////////////////////////////////////////////////////
// Linearized Operator
///////////////////////////////////////////////////////////////////////////////

/// describes a linear mapping X->Y
/**
 * This class is the base class for all linear mappings between two spaces.
 * The domain space and the codomain space are passed as template parameters.
 * The mapping must be linear. For nonlinear mappings see IOperator. The basic
 * usage of this class is to provide the computation of:
 *
 * 		f := L*u,    (resp.  d := J(u) * c in iterative schemes),
 *
 * where f (resp. d) is from the codomain space, u (resp. c) a function of
 * the domain space and L is a linear mapping (resp. J(u) a linearized mapping)
 *
 * This application is splitted into two steps, that have to be called in the
 * correct order:
 *
 * 1. init() or init(u):
 * 		Theses methods initialize the operator for application. One of these
 * 		methods has to be called once before one of the other two methods can be
 * 		invoked. There is no need to init the operator more than once, but
 * 		sometimes - due to parameter change - this is desirable and can be done.
 *
 * 2. apply() or apply_sub():
 * 		These methods can be called when the operator has been initialized
 * 		by a call of init. These function perform the linear mapping, where in
 * 		the case of apply_sub() the result is subtracted from the input function.
 *
 * This splitting has been made, since initialization may be computationally
 * expansive. Thus, the user of this class has the choice when to call this
 * initialization. E.g. when the operator is applied several times, the init
 * of the operator is only needed once.
 *
 * \tparam	X 	Domain space function
 * \tparam	Y	Range space function
 */
template <typename X, typename Y = X>
class ILinearOperator : public IOperator<X,Y>
{
	public:
	///	Domain space
	using domain_function_type = X;

	///	Range space
	using codomain_function_type = Y;

	public:
	///	init operator depending on a function u
	/**
	 * This method initializes the operator. Once initialized the 'apply'-method
	 * can be called. The function u is passed here, since the linear operator
	 * may be the linearization of some non-linear operator. Thus, the operator
	 * depends on the linearization point.
	 * If the operator is not a linearization, this method can be implemented
	 * by simply calling init() and forgetting about the linearization point.
	 *
	 * \param[in]	u		function (linearization point)
	 * \returns 	bool	success flag
	 */
		virtual void init(const X& u) = 0;

	///	init operator
	/**
	 * This method initializes the operator. Once initialized the 'apply'-method
	 * can be called.
	 * If the operator is a linearization this function returns false.
	 *
	 * \returns 	bool	success flag
	 */
		void init() override = 0;

	///	default implementation for IOperator interface
		void prepare(X& u) override {}

	// 	applies the operator
	/**
	 * This method applies the operator, i.e. f = L*u (or d = J(u)*c in
	 * iterative schemes). Note, that the operator must have been initialized
	 * once before this method can be used.
	 *
	 * \param[in]	u		domain function
	 * \param[out]	f		codomain function
	 * \returns		bool	success flag
	 */
		void apply(Y& f, const X& u) override = 0;

	// 	applies the operator and subtracts the result from the input
	/**
	 * This method applies the operator and subracts the result from the input
	 * codomain function, i.e. f -= L*u (or d -= J(u)*c in iterative schemes).
	 * Note, that the operator must have been initialized once before this
	 * method can be used.
	 *
	 * \param[in]		u		domain function
	 * \param[in,out]	f		codomain function
	 * \returns			bool	success flag
	 */
		virtual void apply_sub(Y& f, const X& u) = 0;

	/// virtual	destructor
		~ILinearOperator() override = default;
};

}
#endif