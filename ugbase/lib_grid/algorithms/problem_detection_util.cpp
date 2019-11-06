/*
 * Copyright (c) 2014-2015:  G-CSC, Goethe University Frankfurt
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

#include "problem_detection_util.h"
#include "lib_grid/grid_objects/tetrahedron_rules.h"
#include "common/math/misc/math_util.h"
#include "lib_grid/grid/grid.h"
#include "debug_util.h"
#include "isolated_elements.h"
#include "lib_disc/domain.h"

namespace ug{

int IsSliver(const vector3& v0, const vector3& v1, const vector3& v2,
			  const vector3& v3, number thresholdRatio)
{
	using namespace tet_rules;
	vector3 v[] = {v0, v1, v2, v3};

	number maxLenSq = 0;
	for(int iedge = 0; iedge < NUM_EDGES; ++iedge){
		maxLenSq = VecDistanceSq(v[EDGE_VRT_INDS[iedge][0]], v[EDGE_VRT_INDS[iedge][1]]);
	}

	number thresholdDist = sqrt(maxLenSq) * thresholdRatio;

	for(int iedge = 0; iedge + 1 < NUM_EDGES; ++iedge){
		int iop = OPPOSED_EDGE[iedge];
		if(iop > iedge){
			number dist = DistanceLineToLine(v[EDGE_VRT_INDS[iedge][0]], v[EDGE_VRT_INDS[iedge][1]],
											 v[EDGE_VRT_INDS[iop][0]], v[EDGE_VRT_INDS[iop][1]]);
			if(dist < thresholdDist)
				return iedge;
		}
	}

	return -1;
}


template <class TSide>
static bool CheckForUnconnectedSidesIMPL(Grid& grid, const ISubsetHandler& sh)
{
	bool gotOne = false;
	std::vector<TSide*> sides;
	if(CollectUnconnectedSides(sides,
								 grid,
								 grid.template begin<TSide>(),
								 grid.template end<TSide>()))
	{
		gotOne = true;
		const size_t numSides = sides.size();
		UG_LOG("WARNING: Found " << numSides << " unconnected sides (those may lead to solver issues!):" << std::endl);
		UG_ERR_LOG("Found " << numSides << " unconnected sides (those may lead to solver issues!):" << std::endl);
		for(size_t i = 0; i < numSides; ++i){
			UG_LOG("  - " << ElementDebugInfo(grid, sides[i]) << std::endl);
			UG_ERR_LOG("  - " << ElementDebugInfo(grid, sides[i]) << std::endl);
			UG_LOG("  - " << ElementSubsetInfo(sh, sides[i]) << std::endl);
			UG_ERR_LOG("  - " << ElementSubsetInfo(sh, sides[i]) << std::endl);
		}
	}
	return gotOne;
}

template <typename TDomain>
bool CheckForUnconnectedSides(TDomain& dom)
{
	if(dom.grid().get()->template num<Edge>() > 0 && CheckForUnconnectedSidesIMPL<Vertex>(*dom.grid().get(), *dom.subset_handler().get()))
		return true;
	if(dom.grid().get()->template num<Face>() > 0 && CheckForUnconnectedSidesIMPL<Edge>(*dom.grid().get(), *dom.subset_handler().get()))
		return true;
	if(dom.grid().get()->template num<Volume>() > 0 && CheckForUnconnectedSidesIMPL<Face>(*dom.grid().get(), *dom.subset_handler().get()))
		return true;
	return false;
}

/// explicit template instantiations
template bool CheckForUnconnectedSides<Domain1d>(Domain1d& domain);
template bool CheckForUnconnectedSides<Domain2d>(Domain2d& domain);
template bool CheckForUnconnectedSides<Domain3d>(Domain3d& domain);

}//	end of namespace
