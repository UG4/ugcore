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

#ifndef __H__LIB_ALGEBRA__OPERATOR__INTERFACE__OPERATOR__
#define __H__LIB_ALGEBRA__OPERATOR__INTERFACE__OPERATOR__


namespace ug {

///////////////////////////////////////////////////////////////////////////////
// (Nonlinear-) Operator
///////////////////////////////////////////////////////////////////////////////

/// describes a mapping X->Y
/**
 * This class is the base class for all mappings between two spaces. The domain
 * space and the codomain space are passed as template parameters. In particular,
 * the mapping can be nonlinear. For linear (or linearized) mappings see
 * ILinearizedOperator. The basic usage of this class is to provide the
 * computation of:
 *
 * 		d := N(u),
 *
 * where d is from the codomain space, u a function of the domain space and
 * N() is a (nonlinear-) mapping.
 *
 * This application is splitted into three methods, that have to be called
 * in the correct order:
 *
 * 1. init(): This method initializes the operator for application. It has
 * 			  to be called once before one of the other two methods can be
 * 			  invoked. There is no need to call this method more than once, but
 * 			  sometimes - due to parameter change - this is desirable and
 * 			  can be done.
 *
 * 2. prepare(): This method is used to prepare the in- and output vector used
 * 				 later in apply. It can be called, after the init method has
 * 				 been called at least once. The prepare method is e.g. used
 * 				 to set dirichlet values.
 *
 * 3. apply(): This method can be called when the operator has been initialized
 * 			   by a call of init and with two functions, that have been prepare
 * 			   using the prepare method. It maps the function from the domain
 * 			   space to the range space.
 *
 * This splitting has been made, since initialization and preparation may be
 * computationally expansive. Thus, the user of this class has the choice
 * when to call this initialization/preparation. E.g. when the operator is
 * applied several times on the same vectors, those have only to be prepared
 * once and the init of the operator is only needed once.
 *
 * \tparam	X 	Domain space function
 * \tparam	Y	Range space function
 */
template <typename X, typename Y = X>
class IOperator
{
	public:
	///	Domain space
		using domain_function_type = X;

	///	Range space
		using codomain_function_type = Y;

	public:
	/// initializes the operator
	/**
	 * This method initializes the operator. It must be called before any of
	 * the other methods are called.
	 *
	 * \returns 	bool	success flag
	 */
		virtual void init() = 0;

	/// prepares domain and codomain functions for application
	/**
	 * This method prepares the in- and output functions for the application
	 * and has to be called once before the apply method can be invoked with
	 * the functions used here.
	 *
	 * \param[in]	u		domain function
	 * \returns 	bool	flag if preparation successful
	 */
		virtual void prepare(X& u) = 0;

	///	computes the nonlinear mapping d := N(u)
	/**
	 * This method maps a function from the domain space to the range space.
	 * Note, that is must be called with functions, that have previously been
	 * prepared using the 'prepare'-method and that the operator must have been
	 * initialized using the 'init'-method
	 *
	 * \param[in]	u		domain function
	 * \param[out]	d		codomain function
	 * \returns 	bool	flag if application successful
	 */
		virtual void apply(Y& d, const X& u) = 0;

	///	virtual destructor
		virtual ~IOperator() = default;
};

} // end namespace ug

#endif