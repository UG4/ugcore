#ifndef __H__UG_grid_bridges
#define __H__UG_grid_bridges

#include <string>
#include "registry/registry.h"

namespace ug{
namespace bridge{

void RegisterGridBridge_Grid(Registry& reg, std::string parentGroup);
void RegisterGridBridge_SubsetHandler(Registry& reg, std::string parentGroup);
void RegisterGridBridge_Selector(Registry& reg, std::string parentGroup);
void RegisterGridBridge_Refinement(Registry& reg, std::string parentGroup);
void RegisterGridBridge_Balancing(Registry& reg, std::string parentGroup);
void RegisterGridBridge_FileIO(Registry& reg, std::string parentGroup);
void RegisterGridBridge_Layers(Registry& reg, std::string parentGroup);
void RegisterGridBridge_Debug(Registry& reg, std::string parentGroup);
void RegisterGridBridge_Misc(Registry& reg, std::string parentGroup);

}//	end of namespace	
}//	end of namespace

#endif	//__H__UG_grid_bridges
