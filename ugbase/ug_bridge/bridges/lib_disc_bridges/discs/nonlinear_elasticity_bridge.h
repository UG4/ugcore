/*
 * nonlinear_elasticity_bridge.h
 *
 *  Created on: 08.06.2011
 *      Authors: andreasvogel
 */

#ifndef __H__UG__BRIDGES__NONLINEAR_ELASTICITY_BRIDGE__
#define __H__UG__BRIDGES__NONLINEAR_ELASTICITY_BRIDGE__

namespace ug
{
namespace bridge
{

/// registers needed functionality for the Nonlinear-Elasticity problem
bool RegisterDynamicNonlinearElasticityDisc(Registry& reg, int algebra_type, const char* parentGroup);

} // end namespace ug
} // end namespace bridge

#endif /* __H__UG__BRIDGES__NONLINEAR_ELASTICITY_BRIDGE__ */
