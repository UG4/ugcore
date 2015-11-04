#include "grid_bridges.h"
#include "lib_grid/parallelization/util/partition_weighting_callbacks.h"
#include "lib_grid/tools/partition_map.h"

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
}

}//	end of namespace
}//	end of namespace
