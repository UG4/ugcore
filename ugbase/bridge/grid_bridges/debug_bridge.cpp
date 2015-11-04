#include "grid_bridges.h"
#include "lib_grid/algorithms/debug_util.h"

using namespace std;

namespace ug{
namespace bridge{

void RegisterGridBridge_Debug(Registry& reg, string parentGroup)
{
	string grp = parentGroup;

	reg.add_function("CheckHangingNodeConsistency", static_cast<bool (*)(MultiGrid&)>(&CheckHangingNodeConsistency), grp)
		.add_function("CheckMultiGridConsistency", &CheckMultiGridConsistency, grp)
		.add_function("CheckDistributedObjectConstraintTypes", &CheckDistributedObjectConstraintTypes, grp)
		.add_function("CheckDistributedParentTypes", &CheckDistributedParentTypes, grp)
		.add_function("CheckElementConsistency", static_cast<bool (*)(MultiGrid&, Vertex*)>(&CheckElementConsistency), grp)
		.add_function("CheckElementConsistency", static_cast<bool (*)(MultiGrid&, Edge*)>(&CheckElementConsistency), grp)
		.add_function("CheckElementConsistency", static_cast<bool (*)(MultiGrid&, Face*)>(&CheckElementConsistency), grp);
}

}//	end of namespace
}//	end of namespace
