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

#ifndef __H__LIB_ALGEBRA__PARALLELIZATION__PARALLEL_STORAGE_TYPE__
#define __H__LIB_ALGEBRA__PARALLELIZATION__PARALLEL_STORAGE_TYPE__

#include <iostream>
//#include "common/types.h"

namespace ug {

///\ingroup lib_algebra_parallelization
/// Parallel Storage type
/**
 * The storage type of vector is used in parallel applications.
 * We assume that the dofs are distributed to the processes in the way that
 * each dof is master on exactly one process and can be a slave (i.e. a local
 * copy) on several other processes. Given the real values of the dofs the
 * different storage type are defined as follows:
 *  - PST_UNDEFINED: no information given
 *  - PST_CONSISTENT: The real value is saved in the master and every slave
 *  - PST_ADDITIVE: The sum over the values in the master and all slaves gives the exact value
 *  - PST_UNIQUE: Same as PST_ADDITIV, but value is zero in all slaves (i.e. master has exact value)
 *
 *  Note, that a Vector can have more than one type. E.g. every unique Vector
 *  is additive. Moreover, the Vector being zero everywhere is consistent,
 *  additive and unique at the same time. Therefore, the information is given
 *  bitwise.
 *
 *  To check, whether a Vector v is in a type one may ask:
 *  	 -- v.has_storage_type(PST_CONSISTENT), etc.
 *  To change into another type
 *       -- v.change_storage_type(PST_ADDITIVE)
 */
enum ParallelStorageType : uint
{
	PST_UNDEFINED = 0,
	PST_CONSISTENT = 1 << 0,
	PST_ADDITIVE = 1 << 1,
	PST_UNIQUE = 1 << 2
};

// bitwise and for Parallel Storage Type
inline ParallelStorageType operator & (const ParallelStorageType &a, const ParallelStorageType &b)
{
	return static_cast<ParallelStorageType>(static_cast<uint>(a) & static_cast<uint>(b));
}

inline std::ostream& operator << (std::ostream& outStream, const ParallelStorageType& type)
{
	if(!type) outStream << "undefined";
	if(type & PST_CONSISTENT) outStream << "consistent";
	if(type & PST_UNIQUE) outStream << "unique";
	else if (type & PST_ADDITIVE) outStream << "additive";
	return outStream;
}

} // end namespace ug

#endif