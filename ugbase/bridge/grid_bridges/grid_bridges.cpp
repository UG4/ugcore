// created by Sebastian Reiter
// s.b.reiter@gmail.com

#include "../standard_bridges.h"
#include "grid_bridges.h"
#include "bridge/util.h"

using namespace std;

namespace ug{
namespace bridge{

void RegisterBridge_Grid(Registry& reg, string parentGroup)
{
	try{
		parentGroup.append("/Grid");
		RegisterGridBridge_Grid(reg, parentGroup);
		RegisterGridBridge_SubsetHandler(reg, parentGroup);
		RegisterGridBridge_Selector(reg, parentGroup);
		RegisterGridBridge_Refinement(reg, parentGroup);
		RegisterGridBridge_Balancing(reg, parentGroup);
		RegisterGridBridge_FileIO(reg, parentGroup);
		RegisterGridBridge_Layers(reg, parentGroup);
		RegisterGridBridge_Debug(reg, parentGroup);
		RegisterGridBridge_Misc(reg, parentGroup);
	}
	UG_REGISTRY_CATCH_THROW(parentGroup);
}

}//	end of namespace
}//	end of namespace
