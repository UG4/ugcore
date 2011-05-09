/*
 * lib_disc_bridge_discs.cpp
 *
 *  Created on: 22.09.2010
 *      Author: andreasvogel
 */

// extern headers
#include <iostream>
#include <sstream>

// include bridge
#include "../../../ug_bridge.h"
#include "../../../registry.h"

// other files
#include "navier_stokes_bridge.h"
#include "density_driven_flow_bridge.h"

namespace ug
{

namespace bridge
{

bool RegisterDynamicLibDiscInterfaceDiscs(Registry& reg, int algebra_type, const char* parentGroup)
{
	bool bReturn = true;

	bReturn &= RegisterDynamicNavierStokesDisc(reg, algebra_type, parentGroup);
	bReturn &= RegisterDynamicDensityDrivenFlowDisc(reg, algebra_type, parentGroup);

	return bReturn;
}

}//	end of namespace ug
}//	end of namespace interface
