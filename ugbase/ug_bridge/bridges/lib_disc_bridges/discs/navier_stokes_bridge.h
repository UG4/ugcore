/*
 * navier_stokes_bridge.h
 *
 *  Created on: 09.05.2011
 *      Author: andreasvogel
 */

#ifndef __H__UG__BRIDGES__NAVIER_STOKES_BRIDGE__
#define __H__UG__BRIDGES__NAVIER_STOKES_BRIDGE__

namespace ug
{
namespace bridge
{

/// registers needed functionality for Navier-Stokes-Problem
bool RegisterDynamicNavierStokesDisc(Registry& reg, int algebra_type, const char* parentGroup);

} // end namespace ug
} // end namespace bridge

#endif /* __H__UG__BRIDGES__NAVIER_STOKES_BRIDGE__ */
