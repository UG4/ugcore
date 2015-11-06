/*
 * Copyright (c) 2009-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
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

#ifndef __H__COMMON__STATIC_ASSERT__
#define __H__COMMON__STATIC_ASSERT__

/// \addtogroup ugbase_common
/// \{

//	STATIC_ASSERT is only active during debug-mode.
#ifndef NDEBUG

template <bool> struct CompileTimeAssertion;
template <> struct CompileTimeAssertion<true>
{
	void do_assert()	{};
};

////////////////////////////////////////////////////////////////////////
//	UG_STATIC_ASSERT
///	Checks an expression at compile-time and raises a compile-error if the expression equals 0.
/**
 * UG_STATIC_ASSERT is only active if NDEBUG is not defined.
 *
 * \param expr: an arbitrary expression that can be converted to bool.
 * \param msg: this message will be part of the error-message that the
 * 				compiler prints to the screen.
 * 				Only intended to help the user to identify the error-source.
 * 				The format of the message has to adhere to the rules
 * 				that apply for class-names or variable-names.
 *
 * (bad) example:	UG_STATIC_ASSERT(sizeof(int) == 4, size_of_int_has_to_be_4);
 */
#define UG_STATIC_ASSERT(expr, msg) \
	{CompileTimeAssertion<expr ? true : false> UG_STATIC_ASSERT_ERROR_##msg;\
	UG_STATIC_ASSERT_ERROR_##msg.do_assert();}

#else
//	define an empty STATIC_ASSERT.
#define UG_STATIC_ASSERT(expr, msg)

#endif

// end group ugbase_common
/// \}

#endif
