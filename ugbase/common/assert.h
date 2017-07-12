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

#ifndef __H__COMMON__ASSERT__
#define __H__COMMON__ASSERT__

#include <cassert>
#include "common/log.h"
#include "ug_config.h"

/// \addtogroup ugbase_common
/// \{

#define UG_TO_STRING1(x)  #x
#define UG_TO_STRING(x)  UG_TO_STRING1(x)

////////////////////////////////////////////////////////////////////////
//	UG_ASSERT
///	Checks an expression at runtime and raises a runtime-error if the expression equals 0.
/**
 * UG_ASSERT is only active if NDEBUG is not defined.
 *
 * \param expr: an arbitrary expression that can be converted to bool.
 * \param msg: this message will be part of the error-message. Intended to give
 * 				more information about the checked condition and possible errors.
 *
 * example:		UG_ASSERT( num_dof >= 0, "Number of Degrees of Freedom can not be negative");
 *
 * output:
 *   UG_ASSERT failed:
 *     Condition:   num_dof >= 0
 *     Description: Number of Degrees of Freedom can not be negative
 *     File:        debug.cpp
 *     Line:        130
 *
 */
UG_API void ug_assert_failed();

#ifndef NDEBUG

#define UG_ASSERT(expr, msg)  {if(!(expr)) \
								{ \
									UG_LOG_ALL_PROCS(	"\n  UG_ASSERT failed:\n"\
														"    Condition:   " << UG_TO_STRING(expr) << "\n"\
														"    Description: " << msg << "\n"\
														"    File:        " << __FILE__ << "\n"\
														"    Line:        " << __LINE__ << "\n\n"); \
									ug_assert_failed(); \
									UG_LOG_ALL_PROCS(	"\n  UG_ASSERT failed:\n"\
														"    Condition:   " << UG_TO_STRING(expr) << "\n"\
														"    Description: " << msg << "\n"\
														"    File:        " << __FILE__ << "\n"\
														"    Line:        " << __LINE__ << "\n\n"); \
									assert(expr);\
								}}
#else /* NDEBUG */
#define UG_ASSERT(expr, msg) {}
#endif

// end group ugbase_common
/// \}

#endif /* __H__COMMON__ASSERT__ */
