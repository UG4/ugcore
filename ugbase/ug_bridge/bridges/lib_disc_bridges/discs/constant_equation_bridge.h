/*
 * constant_equation_bridge.h
 *
 *  Created on: 12.05.2011
 *      Author: andreasvogel
 */

#ifndef __H__UG__BRIDGES__CONSTANT_EQUATION_BRIDGE__
#define __H__UG__BRIDGES__CONSTANT_EQUATION_BRIDGE__

namespace ug
{
namespace bridge
{

/// registers needed functionality for Convection-Diffusion Problem
bool RegisterDynamicConstantEquationDisc(Registry& reg, int algebra_type, const char* parentGroup);

} // end namespace ug
} // end namespace bridge

#endif /* __H__UG__BRIDGES__CONSTANT_EQUATION_BRIDGE__ */
