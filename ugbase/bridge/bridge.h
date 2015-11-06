/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
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


#ifndef __H__UG_BRIDGE__UG_BRIDGE__
#define __H__UG_BRIDGE__UG_BRIDGE__

#include <string>
#include <sstream>
#include "registry/registry.h"
#include "lib_algebra/algebra_type.h"
#include "common/ug_config.h"

namespace ug{
namespace bridge{

/**
 * \defgroup bridge Bridge
 * \ingroup ugbase
 * \{
 */

/// string for ug4 group
extern const char* UG4_GRP;

///	returns the default registry used in ug
UG_API Registry & GetUGRegistry();

///	Sets the default classes of class-groups based on a tags using default DoFManager
UG_API void InitUG(int dim, const AlgebraType& algebraType);


/// calls RegisterStandardInterfaces and LoadPlugins if UG_PLUGINS is defined
UG_API void InitBridge();

///	registers all standard interfaces.
/**	This method is called by the constructor of Registry automatically.
 *	You don't have to call it yourself!
 */
UG_API void RegisterStandardBridges(Registry& reg, std::string grp = UG4_GRP);

// end group bridge
/// \}

}//	end bridge
}//	end ug

#include "suffix_tag.h"

#endif /* __H__UG_BRIDGE__UG_BRIDGE__ */
