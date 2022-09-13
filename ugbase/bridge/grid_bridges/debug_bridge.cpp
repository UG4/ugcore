/*
 * Copyright (c) 2015:  G-CSC, Goethe University Frankfurt
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

#include "grid_bridges.h"
#include "lib_grid/algorithms/debug_util.h"
#include "lib_grid/algorithms/problem_detection_util.h"

using namespace std;

namespace ug{
namespace bridge{

template <typename TRegistry=Registry>
void RegisterGridBridge_Debug_(TRegistry& reg, string parentGroup)
{
	string grp = parentGroup;

	reg.add_function("CheckHangingNodeConsistency", static_cast<bool (*)(MultiGrid&)>(&CheckHangingNodeConsistency), grp)
		.add_function("CheckMultiGridConsistency", &CheckMultiGridConsistency, grp)
		.add_function("CheckDistributedObjectConstraintTypes", &CheckDistributedObjectConstraintTypes, grp)
		.add_function("CheckDistributedParentTypes", &CheckDistributedParentTypes, grp)
		.add_function("CheckElementConsistency", static_cast<bool (*)(MultiGrid&, Vertex*)>(&CheckElementConsistency), grp)
		.add_function("CheckElementConsistency", static_cast<bool (*)(MultiGrid&, Edge*)>(&CheckElementConsistency), grp)
		.add_function("CheckElementConsistency", static_cast<bool (*)(MultiGrid&, Face*)>(&CheckElementConsistency), grp);

	reg.add_function("CheckForUnconnectedSides", &CheckForUnconnectedSides,
					 grp, "foundUnconnectedSides", "grid",
					 "Checks whether unconnected sides exist in the given grid.");
}

}//	end of namespace

UG_REGISTRY_DEFINE(RegisterGridBridge_Debug);

}//	end of namespace
