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
#include "lib_grid/parallelization/util/partition_weighting_callbacks.h"
#include "lib_grid/tools/partition_map.h"
#include "lib_grid/algorithms/graph/dual_graph.h"	// DualGraphNeighborCollector

using namespace std;

namespace ug{
namespace bridge{

void RegisterGridBridge_Balancing(Registry& reg, string parentGroup)
{
	string grp = parentGroup;

// partition weighting in metis partitioning
	reg.add_class_<PartitionWeighting>("PartitionWeighting", grp)
		.add_constructor()
		.add_method("set_default_weights", &PartitionWeighting::set_default_weights, "", "hWeight#vWeight")
		.set_construct_as_smart_pointer(true);
	reg.add_class_<InterSubsetPartitionWeighting, PartitionWeighting>("InterSubsetPartitionWeighting", grp)
		.add_constructor()
		.add_method("set_inter_subset_weight", &InterSubsetPartitionWeighting::set_inter_subset_weight, "", "si1#si2#weight")
		.set_construct_as_smart_pointer(true);
	reg.add_class_<ProtectSubsetPartitionWeighting, PartitionWeighting>("ProtectSubsetPartitionWeighting", grp)
		.add_constructor()
		.add_method("set_weight", &ProtectSubsetPartitionWeighting::set_weight, "", "si#weight")
		.set_construct_as_smart_pointer(true);

//	PartitionMap
	reg.add_class_<PartitionMap>("PartitionMap", grp)
		.add_constructor()
		.add_method("clear", &PartitionMap::clear)
		.add_method("get_partition_handler", &PartitionMap::get_partition_handler)
		.add_method("add_target_proc", &PartitionMap::add_target_proc)
		.add_method("add_target_procs", &PartitionMap::add_target_procs)
		.add_method("num_target_procs", &PartitionMap::num_target_procs)
		.add_method("get_target_proc", &PartitionMap::get_target_proc)
		.add_method("shift_target_procs", &PartitionMap::shift_target_procs)
		.set_construct_as_smart_pointer(true);

//	dual graph neighbor collector
	{
		string name = string("IDualGraphNeighborCollector");
		typedef IDualGraphNeighborCollector T;
		reg.add_class_<T>(name, grp);
	}
}

}//	end of namespace
}//	end of namespace
