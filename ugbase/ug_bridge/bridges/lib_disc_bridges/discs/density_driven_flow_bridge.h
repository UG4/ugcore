/*
 * density_driven_flow_bridge.h
 *
 *  Created on: 09.05.2011
 *      Author: andreasvogel
 */

#ifndef __H__UG__BRIDGES__DENSITY_DRIVEN_FLOW_BRIDGE__
#define __H__UG__BRIDGES__DENSITY_DRIVEN_FLOW_BRIDGE__

namespace ug
{
namespace bridge
{

/// registers needed functionality for Density-Driven-Flow Problem
bool RegisterDynamicDensityDrivenFlowDisc(Registry& reg, int algebra_type, const char* parentGroup);

} // end namespace ug
} // end namespace bridge

#endif /* __H__UG__BRIDGES__DENSITY_DRIVEN_FLOW_BRIDGE__ */
