/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
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

#include "common/common.h"
#include "common/cuthill_mckee.h"
#include "cuthill_mckee.h"
#include <algorithm>
#include <vector>
#include <queue>
#include "common/profiler/profiler.h"
#include "lib_disc/domain.h"

namespace ug{

void OrderCuthillMcKee(DoFDistribution& dofDistr, bool bReverse)
{
	PROFILE_FUNC();
//	get adjacency graph
	std::vector<std::vector<size_t> > vvConnection;
	try{
		dofDistr.get_connections(vvConnection);
	}
	UG_CATCH_THROW("OrderCuthillMcKee: No adjacency graph available.");

//	get mapping for cuthill-mckee order
	std::vector<size_t> vNewIndex;
	ComputeCuthillMcKeeOrder(vNewIndex, vvConnection, bReverse);

//	reorder indices
	dofDistr.permute_indices(vNewIndex);
}

template <typename TDomain>
void OrderCuthillMcKee(ApproximationSpace<TDomain>& approxSpace, bool bReverse)
{
	std::vector<SmartPtr<DoFDistribution> > vDD = approxSpace.dof_distributions();

	for(size_t i = 0; i < vDD.size(); ++i)
		OrderCuthillMcKee(*vDD[i], bReverse);
}

#ifdef UG_DIM_1
template void OrderCuthillMcKee<Domain1d>(ApproximationSpace<Domain1d>& approxSpace, bool bReverse);
#endif
#ifdef UG_DIM_2
template void OrderCuthillMcKee<Domain2d>(ApproximationSpace<Domain2d>& approxSpace, bool bReverse);
#endif
#ifdef UG_DIM_3
template void OrderCuthillMcKee<Domain3d>(ApproximationSpace<Domain3d>& approxSpace, bool bReverse);
#endif

}
