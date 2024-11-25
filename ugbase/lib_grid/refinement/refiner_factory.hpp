/*
 * Copyright (c) 2024: Goethe University Frankfurt
 * Author: Arne Naegel
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

#ifndef __H__UG__DOMAIN_REFINEMENT_BRIDGE_HPP__
#define __H__UG__DOMAIN_REFINEMENT_BRIDGE_HPP__

#include "lib_grid/lib_grid.h"
#include "lib_grid/algorithms/refinement_mark_util.h"
#include "lib_grid/refinement/adaptive_regular_mg_refiner.h"
#include "lib_grid/refinement/global_fractured_media_refiner.h"
#include "lib_grid/refinement/global_multi_grid_refiner.h"
#include "lib_grid/refinement/global_subdivision_multi_grid_refiner.h"
#include "lib_grid/refinement/hanging_node_refiner_multi_grid.h"
#include "lib_grid/refinement/ref_mark_adjusters/horizontal_anisotropy_adjuster.h"
#include "lib_grid/refinement/ref_mark_adjusters/shadow_copy_adjuster.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
///	Creates a global domain refiner.
/**	Automatically chooses whether a parallel refiner is required.*/
template <typename TDomain>
static SmartPtr<IRefiner> GlobalDomainRefiner(TDomain& dom)
{
//todo: support normal grids, too!

	#ifdef UG_PARALLEL
		if(pcl::NumProcs() > 1){
			return SmartPtr<IRefiner>(
					new ParallelGlobalRefiner_MultiGrid(
							*dom.distributed_grid_manager(),
							dom.refinement_projector()));
		}
	#endif

	return SmartPtr<IRefiner>(
				new GlobalMultiGridRefiner(
						*dom.grid(),
						dom.refinement_projector()));
}

}
#endif
