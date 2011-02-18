/*
 * stabilization_fields.h
 *
 *  Created on: 20.01.2011
 *      Author: josefdubsky
 */

#ifndef __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__NAVIER_STOKES__FV__STABILIZATION_FIELDS__
#define __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__NAVIER_STOKES__FV__STABILIZATION_FIELDS__

#include <limits>

#include "lib_discretization/common/geometry_util.h"

namespace ug{

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
bool GetFieldsStabilizedShapes(	const TFVGeometry& geo,
                                    const MathVector<TFVGeometry::world_dim> vCornerVels[TFVGeometry::m_numSCV],
                                    number vCornerPress[TFVGeometry::m_numSCV],
                                    const int StabMethod,
                                    const MathVector<TFVGeometry::world_dim> vIPVelUpwindShapesContiEq[TFVGeometry::m_numSCVF][TFVGeometry::m_numSCV][TFVGeometry::world_dim],
                                    const number vConvLength[TFVGeometry::m_numSCVF],
                                    const number dt,
                                    bool bTimeDependent,
                                    const MathVector<TFVGeometry::world_dim> vCornerVelsOld[TFVGeometry::m_numSCV],
 									number kinematicViscosity,
                                    MathVector<TFVGeometry::world_dim> vIPStabVelShapesContiEq[TFVGeometry::m_numSCVF][TFVGeometry::m_numSCV][(TFVGeometry::world_dim)+1])
{
    // todo: switch
	// Compute ConvectionLength hereCompute Diffusion Length
	//NSDiffLengthMethod1(vDiffLengthSqInv, geo);


//	todo: Implement time-dependent part
	if(bTimeDependent)
	{
		UG_LOG("ERROR in 'GetFieldsStabilizedVelocitiesDiagonal':"
				"Time-dependent version not implemented.\n");
		return false;
	}

//	Some constants
	static const size_t numIp = TFVGeometry::m_numSCVF;
	static const size_t dim = TFVGeometry::world_dim;
    MathVector<TFVGeometry::world_dim> vIPVelCurrent[TFVGeometry::m_numSCVF];

//	Vector for Diagonal (size = dim * NumIps)
	number Diag[numIp * dim];
	number rhs[numIp * dim];

//	Compute diffusion length
//	todo: compute
	number DiffLengthSq = 1.0;

// Compute velocity in the integration points.
    //	Loop integration points
    for(size_t ip = 0; ip < numIp; ++ip)
	{
    	//	get SubControlVolumeFace
		const typename TFVGeometry::SCVF& scvf = geo.scvf(ip);
		UG_ASSERT(scvf.num_ip() == 1, "Only implemented for first order");

        // reset values to zero
        VecSet(vIPVelCurrent[ip], 0.0);

	// 	Loop components of velocity
       for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
        {
           VecScaleAppend(vIPVelCurrent[ip], scvf.shape(sh, 0), vCornerVels[sh]);
        }
    }

//	Loop integration points
	for(size_t ip = 0; ip < numIp; ++ip)
	{
	// 	Loop components of velocity
		for(size_t d = 0; d < dim; d++)
		{
		//	current component of vector
			const size_t comp = ip * dim + d;

		///////////////////////
		//	assemble diagonal
		////////////////////////

		//	Time part
			if(bTimeDependent)
				Diag[comp] = 1./dt;
			else
				Diag[comp] = 0.0;

		//	Diffusion part
			Diag[comp] += kinematicViscosity/DiffLengthSq;

		//	Convective Term
			Diag[comp] += VecTwoNorm(vIPVelCurrent[ip]) / vConvLength[ip];

		//////////////////////////////
		//	assemble right-hand side
		//////////////////////////////

		//	Source
		//	todo: Implement source/model term
			rhs[comp] = 0.0;

		//	Time
			if(bTimeDependent)
		//		rhs[comp] += vIPVelOld[ip][d] / dt;

		//	Diffusion part
			rhs[comp] += kinematicViscosity * vIPVelCurrent[ip][d] / DiffLengthSq;

		//	Convective part
		//	rhs[comp] += vIPVelUpwindShapesContiEq[ip][d][0][0] * VecTwoNorm(vIPVelCurrent[ip]) / ConvLength;

		//	Pressure term
			//rhs[comp] -= IPPressureGrad[ip][d];
		}
	}

	///////////////////////
	// Solve System
	///////////////////////

	for(size_t ip = 0; ip < numIp; ++ip)
	{
		for(size_t d = 0; d < dim; d++)
		{
		//	current component of vector
			const size_t comp = ip * dim + d;

		//	invert diagonal
			vIPStabVelShapesContiEq[ip][d][0][0] = rhs[comp] / Diag[comp];
		}
	}

//	we're done
	return true;
}
    
} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__NAVIER_STOKES__FV__STABILIZATION_FIELDS__ */
