/*
 * Copyright (c) 2010-2013:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
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

#ifndef __H__UG__CPU_ALGEBRA__ALGEBRA_MISC__
#define __H__UG__CPU_ALGEBRA__ALGEBRA_MISC__
////////////////////////////////////////////////////////////////////////////////////////////////

namespace ug
{

/// \addtogroup lib_algebra
/// \{

#ifndef NDEBUG
//!
//! use this to force the creation of prsize_t routines or similar for use in gdb.
#define FORCE_CREATION volatile size_t ___never_happens___ = 0; if(___never_happens___)
#else
#define FORCE_CREATION if(0)
#endif

//! prevent unused variable-warnings
#define UNUSED_VARIABLE(var) ((void) var);


//!
//! template struct for sorting some keys after values
//! for example, sorting a vector of ints and know original pos
template<typename T>
struct sortStruct
{
	size_t index; // for example "original" position.
	T sortValue;

	bool operator < (const sortStruct<T> &other) const
	{
		return sortValue < other.sortValue;
	}
};

// end group lib_algebra
/// \}

}


#endif
