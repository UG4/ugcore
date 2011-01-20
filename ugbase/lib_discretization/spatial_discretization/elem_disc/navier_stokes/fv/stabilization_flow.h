/*
 * stabilization_flow.h
 *
 *  Created on: 20.01.2011
 *      Author: josefdubsky
 */

#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__NAVIER_STOKES__FV__STABILIZATION_FLOW__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__NAVIER_STOKES__FV__STABILIZATION_FLOW__

#include <limits>

#include "lib_discretization/common/geometry_util.h"

namespace ug{

/**
 *
 *
 * \param[in]	geo                         Finite Volume Geometry
 * \param[in]	vCornerValues               Solution values at corners from last iteration
 * \param[in]	vIPVelCurrent               Solution values at IPs from last iteration
 * \param[in]	StabMethod                  Defines which stabilization method should be used
 * \param[in] 	vIPVelUpwindShapesContiEq	Upwind shapes in the IPs
 * \param[in] 	vConvLength                 Convective length corresponding to the Upwind shapes in the IPs
 * \param[in] 	vConvLength                 Convective length corresponding to the Upwind shapes in the IPs
 * \param[in]	dt                          time step size
 * \param[in]	bTimeDependent              flag indicating transient model
 * \param[in]	IPVelOld                    Velocity at ips from old timestep
 * \param[in]	kinematicViscosity          kinematic Viscosity
 * \param[out]	vIPStabVelShapesContiEq     Stabilized velocity shapes in the IPs
 */
template <typename TFVGeometry>
bool GetFlowStabilizedShapes(	const TFVGeometry& geo,
                                    const MathVector<TFVGeometry::world_dim> vCornerValues[TFVGeometry::m_numSCVF],
                                    const MathVector<TFVGeometry::world_dim> vIPVelCurrent[TFVGeometry::m_numSCVF],
                                    const int StabMethod,
                                    const MathVector<TFVGeometry::world_dim> vIPVelUpwindShapesContiEq[TFVGeometry::m_numSCVF][TFVGeometry::m_numSCV][TFVGeometry::world_dim],
                                    const number vConvLength[TFVGeometry::m_numSCVF],
                                    const number dt,
                                    bool bTimeDependent,
                                    const MathVector<TFVGeometry::world_dim> vIPVelOld[TFVGeometry::m_numSCVF],
 									number kinematicViscosity,
                                    MathVector<TFVGeometry::world_dim> vIPStabVelShapesContiEq[TFVGeometry::m_numSCVF][TFVGeometry::m_numSCV][(TFVGeometry::world_dim)+1])
{
    // todo: Implement Flow stabilization
    return false;
}

} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__NAVIER_STOKES__FV__STABILIZATION_FLOW__ */
