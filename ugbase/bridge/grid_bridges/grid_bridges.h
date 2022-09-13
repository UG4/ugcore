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

#ifndef __H__UG_grid_bridges
#define __H__UG_grid_bridges

#include <string>
#include "registry/registry.h"
#include "bridge/bridge.h"

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

namespace pybind{

void RegisterGridBridge_Grid(ug::pybind::RegistryAdapter& reg, std::string parentGroup);
void RegisterGridBridge_SubsetHandler(ug::pybind::RegistryAdapter& reg, std::string parentGroup);
void RegisterGridBridge_Selector(ug::pybind::RegistryAdapter& reg, std::string parentGroup);
void RegisterGridBridge_Refinement(ug::pybind::RegistryAdapter& reg, std::string parentGroup);
void RegisterGridBridge_Balancing(ug::pybind::RegistryAdapter& reg, std::string parentGroup);
void RegisterGridBridge_FileIO(ug::pybind::RegistryAdapter& reg, std::string parentGroup);
void RegisterGridBridge_Layers(ug::pybind::RegistryAdapter& reg, std::string parentGroup);
void RegisterGridBridge_Debug(ug::pybind::RegistryAdapter& reg, std::string parentGroup);
void RegisterGridBridge_Misc(ug::pybind::RegistryAdapter& reg, std::string parentGroup);

}//	end of namespace pybind

UG_REGISTRY_DECL(RegisterGridBridge_Grid);

}//	end of namespace

#endif	//__H__UG_grid_bridges
