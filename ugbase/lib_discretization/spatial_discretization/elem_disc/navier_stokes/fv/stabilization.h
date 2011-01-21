/*
 * stabilization.h
 *
 *  Created on: 28.10.2010
 *      Author: josefdubsky, andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__NAVIER_STOKES__FV__STABILIZATION__
#define __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__NAVIER_STOKES__FV__STABILIZATION__

// other ug4 modules
#include "common/common.h"
#include "lib_discretization/common/geometry_util.h"
#include "stabilization_fields.h"
#include "stabilization_flow.h"

namespace ug{
		enum STABILIZATION_TYPES
		{
            FIELDS = 0,
            FLOW
		};

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
bool GetStabilizedShapes(	const TFVGeometry& geo,
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

    // Compute Upwind Shapes at Ip's and ConvectionLength here
	switch(StabMethod)
	{
		case FIELDS:   if(!GetFieldsStabilizedShapes(geo, vCornerValues, vIPVelCurrent, StabMethod, vIPVelUpwindShapesContiEq, vConvLength,
                                                    dt, bTimeDependent, vIPVelOld, kinematicViscosity, vIPStabVelShapesContiEq))
                                return false;
                            break;

		case FLOW:      if(!GetFlowStabilizedShapes(geo, vCornerValues, vIPVelCurrent, StabMethod, vIPVelUpwindShapesContiEq, vConvLength,
                                                    dt, bTimeDependent, vIPVelOld, kinematicViscosity, vIPStabVelShapesContiEq))
                                return false;
                            break;


        default: 	UG_LOG("Stabilization Type defined incorrecrly.\n");
					return false;
	}

    // Values of Velocities in IPs do not depend in other IPs
	return true;

}

} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__NAVIER_STOKES__FV__STABILIZATION__ */

