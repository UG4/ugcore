/*
 * upwind.h
 *
 *  Created on: 29.10.2010
 *      Author: josefdubsky, andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__NAVIER_STOKES__FV__UPWIND__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__NAVIER_STOKES__FV__UPWIND__

// other ug4 modules
#include "common/common.h"
#include "lib_discretization/common/geometry_util.h"
#include "lib_discretization/spacial_discretization/disc_helper/upwind_shapes.h"

namespace ug{

/**
 *
 * \param[in]		geo					Finite Volume Geometry
 * \param[in]		IPVel				Velocity at Integration points (ip)
 * \param[out]		CornerShape			Shape functions to compute Upwind Velocity
 * \param[out]		ConvectionLength	Distance from IP to Upwind point
 * \param[out]		IPShape				factor how much Upwind velocity depend on other IPVel
 * \param[out]		bDependOnOIP		true if and only if IPShape is non-zero
 */
template <typename TFVGeometry>
bool GetFullUpwindShapesDependingOnIP(	const TFVGeometry& geo,
										const MathVector<TFVGeometry::world_dim> IPVel[],
										std::vector<std::vector<number> >& CornerShape,
										number ConvectionLength[],
										const number IPScaleNumber[],
										bool &bDependOnIP)
{
	if(!GetFullUpwindShapes(geo, IPVel, CornerShape))
		return false;

	// does not depend in IP Velocities
	bDependOnIP = false;

	return true;
}



} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__NAVIER_STOKES__FV__UPWIND__ */
