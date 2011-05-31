/*
 * thermohaline_flow_bridge.h
 *
 *  Created on: 20.05.2011
 *      Author: andreasvogel
 */

#ifndef __H__UG__BRIDGES__THERMOHALINE_FLOW_BRIDGE__
#define __H__UG__BRIDGES__THERMOHALINE_FLOW_BRIDGE__

namespace ug
{
namespace bridge
{

/// registers needed functionality for Thermohaline-Flow Problem
bool RegisterDynamicThermohalineFlowDisc(Registry& reg, int algebra_type, const char* parentGroup);

} // end namespace ug
} // end namespace bridge

#endif /* __H__UG__BRIDGES__THERMOHALINE_FLOW_BRIDGE__ */
