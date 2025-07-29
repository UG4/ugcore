/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__LIB_GRID__PARALLELIZATION_UTIL__
#define __H__LIB_GRID__PARALLELIZATION_UTIL__

#include "distributed_grid.h"
//#include <boost/function.hpp>

#define PROFILE_GRID_DISTRIBUTION
#ifdef PROFILE_GRID_DISTRIBUTION
	#define GDIST_PROFILE_FUNC()	PROFILE_FUNC_GROUP("gdist")
	#define GDIST_PROFILE(name)	PROFILE_BEGIN_GROUP(name, "gdist")
	#define GDIST_PROFILE_END()	PROFILE_END()
	#define GDIST_PROFILE_END_(name)	PROFILE_END_(name)
#else
	#define GDIST_PROFILE_FUNC()
	#define GDIST_PROFILE(name)
	#define GDIST_PROFILE_END()
#endif

namespace ug
{
/// \addtogroup lib_grid_parallelization
///	@{
////////////////////////////////////////////////////////////////////////
//	utility methods
///	Returns the type of associated interfaces
/**
 * \param interfaceType: a constant enumerated in ug::InterfaceNodeTypes.
 * \return the associated interface type or INT_NONE, if no associated
 *		type is known.
 */
int GetAssociatedInterfaceType(int interfaceType);


////////////////////////////////////////////////////////////////////////
///	Creates and distributes global ids for the given element type.
/**	IDs are written to the given attachment (aGeomObjID by default).
 */
template <class TGeomObj>
void CreateAndDistributeGlobalIDs(Grid& g, GridLayoutMap& glm,
								  AGeomObjID& aID = aGeomObjID);

////////////////////////////////////////////////////////////////////////
///	Checks whether the grid-layout-map on this proc is consistent with connected ones.
bool TestGridLayoutMap(MultiGrid& mg, GridLayoutMap& glm, bool verbose = true);

///	@}
}//	end of namespace


////////////////////////////////
//	include implementation
#include "parallelization_util_impl.hpp"

#endif
