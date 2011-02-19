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
#include "diffusion_length.h"

namespace ug{
        enum STABILIZATION_SETUP
        {
            STAB_METHOD = 0,
            DIFF_LENGTH_METHOD,
            PAC_SWITCH,
            PECLET_BLENDING_SWITCH,
            STAB_OPTIONS
        };

		enum STABILIZATION_TYPES
		{
            FIELDS = 0,
            FLOW
		};

        enum DIFFUSION_LENGTH_TYPES
		{
            RAW = 0,
            FIVEPOINT,
            COR
		};

        enum PHYSICAL_ADVECTION_CORRECTION
		{
            NOPAC = 0,
            PAC
		};

        enum PECLET_BLENDING
		{
            NOPEBLEND = 0,
            PEBLEND
		};

/**
 *
 *
 * \param[in]	geo                         Finite Volume Geometry
 * \param[in]	vCornerVels               Solution values at corners from last iteration
 * \param[in]	StabMethod                  Defines which stabilization method should be used
 * \param[in] 	vIPVelUpwindShapesContiEq	Upwind shapes in the IPs
 * \param[in] 	vConvLength                 Convective length corresponding to the Upwind shapes in the IPs
 * \param[in]	dt                          time step size
 * \param[in]	bTimeDependent              flag indicating transient model
 * \param[in]	vCornerVelsOld            Velocity in coners from old timestep
 * \param[in]	kinematicViscosity          kinematic Viscosity
 * \param[out]	vIPStabVelShapesContiEq     Stabilized velocity shapes in the IPs
 */

template <typename TFVGeometry>
bool GetStabilizedShapes(	const TFVGeometry& geo,
                                    const MathVector<TFVGeometry::world_dim> vCornerVels[TFVGeometry::m_numSCV],
                                    number vCornerPress[TFVGeometry::m_numSCV],
                                    const int StabOptions[],
                                    const number vIPVelUpwindShapesContiEq[TFVGeometry::m_numSCVF][TFVGeometry::m_numSCV],
                                    const number vIPVelUpwindDependenciesContiEq[TFVGeometry::m_numSCVF][TFVGeometry::m_numSCVF],
                                    const number vConvLength[TFVGeometry::m_numSCVF],
                                    const number dt,
                                    bool bTimeDependent,
                                    const MathVector<TFVGeometry::world_dim> vCornerVelsOld[TFVGeometry::m_numSCV],
 									number kinematicViscosity,
                                    MathVector<TFVGeometry::world_dim> vIPStabVelShapesContiEq[TFVGeometry::m_numSCVF][TFVGeometry::m_numSCV][(TFVGeometry::world_dim)+1])
{
    number vDiffLengthSqInv[TFVGeometry::m_numSCVF];
    
    // Compute Diffusion Length in corresponding IPs
	switch(StabOptions[DIFF_LENGTH_METHOD])
	{
		case FIVEPOINT: if(!NSDiffLengthFivePoint(vDiffLengthSqInv, geo))
                            return false;
                        break;

		case RAW:       if(!NSDiffLengthRaw(vDiffLengthSqInv, geo))
                            return false;
                        break;

		case COR:       if(!NSDiffLengthCor(vDiffLengthSqInv, geo))
                            return false;
                        break;

        default: 	UG_LOG("Diffusion Length type defined incorrecrly.\n");
					return false;
	}


    // Compute Upwind Shapes at Ip's here
	switch(StabOptions[STAB_METHOD])
	{
		case FIELDS:   if(!GetFieldsStabilizedShapes(geo, vCornerVels, vCornerPress, StabOptions,
                                                    vIPVelUpwindShapesContiEq,vIPVelUpwindDependenciesContiEq, vConvLength,vDiffLengthSqInv,
                                                    dt, bTimeDependent, vCornerVelsOld, kinematicViscosity, vIPStabVelShapesContiEq))
                                return false;
                            break;

		case FLOW:      if(!GetFlowStabilizedShapes(geo, vCornerVels, vCornerPress, StabOptions,
                                                    vIPVelUpwindShapesContiEq,vIPVelUpwindDependenciesContiEq, vConvLength,vDiffLengthSqInv,
                                                    dt, bTimeDependent, vCornerVelsOld, kinematicViscosity, vIPStabVelShapesContiEq))
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

