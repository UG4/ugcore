/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__LIB_GRID__GEOMETRIC_OBJECT_ID__
#define __H__LIB_GRID__GEOMETRIC_OBJECT_ID__

//#include <algorithm>
#include <ostream>
// #include "lib_grid/grid/grid_base_objects.h"
#include "common/util/hash_function.h"

namespace ug{

/// \addtogroup lib_grid_parallelization
/// @{

////////////////////////////////////////////////////////////////////////
///	The global id can be used to uniquely identify distributed objects.
/**	Note that normally no global IDs are associated with geometric objects.
 * However, methods exist, which can assign global ids based on the
 * current layouts.
 */
using GeomObjID = std::pair<int, size_t>;

UG_API std::ostream& operator << (std::ostream& out, const GeomObjID& goId);

UG_API bool operator < (const GeomObjID& gid1, const GeomObjID& gid2);

////////////////////////////////////////////////////////////////////////
///	Can be used to construct a GeomObjID from a proc-rank and a local id.
inline GeomObjID MakeGeomObjID(int procRank, size_t localGeomObjID)
{
	return std::make_pair(procRank, localGeomObjID);
}

////////////////////////////////////////////////////////////////////////
///	An attachment which can store GeomObjIDs
using AGeomObjID = Attachment<GeomObjID>;

///	This attachment instance should be used to store global ids
extern AGeomObjID aGeomObjID;

///	generates a hash key for a GeomObjID.
/**	\todo Check distribution quality.*/
template <>
size_t hash_key<GeomObjID>(const GeomObjID& key);

}// end of namespace

#endif
