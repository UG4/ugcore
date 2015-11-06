/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__FUNCTION_CAST__
#define __H__FUNCTION_CAST__

/// \addtogroup ugbase_common_util
/// \{

///	Casts a function of one type to a function of another type
/**	ATTENTION: Use with care!
 * This method should only be used on functions which have similar parameters,
 * and where different parameters in TTarget and TSrc are pointers or references
 * to types which are connected through inheritance.
 *
 * The main purpose of this method is to make code readable and to allow
 * more secure compiler dependent implementations.
 *
 * \remark	The cast is currently performed through a c-style cast, since a
 *			static_cast or reinterpret_cast for the described purpose is not
 *			raises a compile error in Microsoft Visual Studio 2010.*/
template <typename TTarget, typename TSrc>
TTarget function_cast(TSrc src)
{
	return (TTarget)src;
}

// end group ugbase_common_util
/// \}

#endif
