/*
 * upwind.h
 *
 *  Created on: 29.10.2010
 *      Author: josefdubsky, andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__NAVIER_STOKES__FV__UPWIND__
#define __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__NAVIER_STOKES__FV__UPWIND__

// other ug4 modules
#include "common/common.h"
#include "lib_discretization/common/geometry_util.h"
#include "lib_discretization/spatial_discretization/disc_helper/upwind_shapes.h"

namespace ug{
		enum UPWIND_TYPES
		{
            NO_UPWIND = 0,
            FULL_UPWIND,
			LPS_UPWIND,
			NUM_UPWIND
		};

/**
 *
 * \param[in]		geo                     Finite Volume Geometry
 * \param[in]		CornerValues            Velocity at element corners
 * \param[out]		vvvIPVelUpwindShapes    Shape functions to compute Upwind Velocity
 * \param[out]		ConvectionLength        Distance from IP to Upwind point
 * \param[out]		IPShape                 factor how much Upwind velocity depend on other IPVel
 * \param[out]		bDependOnOIP            true if and only if IPShape is non-zero
 */
template <typename TFVGeometry>
bool GetUpwindShapes(	const TFVGeometry& geo,
                                        const MathVector<TFVGeometry::world_dim> vCornerVels[TFVGeometry::m_numSCV],
                                        const int UpwindMethod,
                                        number vIPVelUpwindShapes[TFVGeometry::m_numSCVF][TFVGeometry::m_numSCV],
                                        number vIPVelUpwindDependencies[TFVGeometry::m_numSCVF][TFVGeometry::m_numSCVF],
                                        number ConvectionLength[TFVGeometry::m_numSCV])
{
// todo: implement the case whe upwind vels in the ips are mutually dependent
//    MathVector<TFVGeometry::world_dim> vIPVelUpwindCornerShapes[TFVGeometry::m_numSCVF][TFVGeometry::m_numSCV];
//    MathVector<TFVGeometry::world_dim> vIPVelUpwindIPShapes[TFVGeometry::m_numSCVF][TFVGeometry::m_numSCV];



    // Compute Upwind Shapes at Ip's and ConvectionLength here
	switch(UpwindMethod)
	{
		case NO_UPWIND:     if(!GetNoUpwindShapes(geo, vCornerVels, vIPVelUpwindShapes,vIPVelUpwindDependencies, ConvectionLength))
                                return false;
                            break;
        case FULL_UPWIND:   if(!GetFullUpwindShapes(geo, vCornerVels, vIPVelUpwindShapes,vIPVelUpwindDependencies, ConvectionLength))
                                return false;
                            break;
        default: 	UG_LOG("Upwind Type defined incorrecrly.\n");
					return false;
	}

    // Values of Velocities in IPs do not depend in other IPs
	return true;
}

} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__NAVIER_STOKES__FV__UPWIND__ */
