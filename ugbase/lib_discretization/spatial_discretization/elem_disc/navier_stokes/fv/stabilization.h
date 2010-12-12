/*
 * stabilization.h
 *
 *  Created on: 28.10.2010
 *      Author: josefdubsky, andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__NAVIER_STOKES__FV__STABILIZATION__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__NAVIER_STOKES__FV__STABILIZATION__

// other ug4 modules
#include "common/common.h"

namespace ug{

/**
 *
 *
 * \param[out] 	IPVelStab			Stabilized velocity at integration points (ip)
 * \param[out] 	IPStabVelShape		Derivative of Stab Velocity with respect to Corner Velocities
 * \param[out] 	IPStabPressureShape	Derivative of Stab Velocity with respect to Corner Pressure
 * \param[in]	geo					Finite Volume Geometry
 * \param[in]	CurrentIPVel		Velocity at ips from last iterate
 * \param[in]	bTimeDependent		flag indicating transient model
 * \param[in]	IPVelOld			Velocity at ips from old timestep
 * \param[in]	dt					time step size
 * \param[in]	IPVelUpwind			upwind Velocity at ips
 * \param[in]	bUpwindDependOnIP	flag indicating iff upwind depends in ip values
 * \param[in]	UpwindScalar		upwind scaling factors
 * \param[in]	CornerVel			Velocity at corners from last iterate
 * \param[in]	IPPressureGrad		Derivative of Pressure at ips from last iterate
 * \param[in]	kinematicViscosity	kinematic Viscosity
 */
template <typename TFVGeometry>
bool GetFieldsStabilizedVelocitiesDiagonal
(
		MathVector<TFVGeometry::world_dim> IPVelStab[],
		number IPStabVelShape[][TFVGeometry::world_dim],
		number IPStabPressureShape[][TFVGeometry::world_dim],
		const TFVGeometry& geo,
		MathVector<TFVGeometry::world_dim> CurrentIPVel[],
		bool bTimeDependent,
		MathVector<TFVGeometry::world_dim> IPVelOld[],
		number dt,
		MathVector<TFVGeometry::world_dim> IPVelUpwind[],
		bool bUpwindDependOnIP,
		number UpwindScalar[],
		MathVector<TFVGeometry::world_dim> CornerVel[],
		MathVector<TFVGeometry::world_dim> IPPressureGrad[],
		number kinematicViscosity
)
{
//	Check that upwind is diagonal (i.e. Upwind does only depend on corners)
	if(bUpwindDependOnIP == true)
	{
		UG_LOG("ERROR in 'GetFieldsStabilizedVelocitiesDiagonal':"
				"Using Diagonal computation for non diagonal upwind.\n");
		return false;
	}

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

//	Vector for Diagonal (size = dim * NumIps)
	number Diag[numIp * dim];
	number rhs[numIp * dim];

//	Compute diffusion length
//	todo: compute
	number DiffLengthSq = 1.0;

//	Compute convection length
	number ConvLength = 1.0;


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
			Diag[comp] += VecTwoNorm(CurrentIPVel[ip]) / ConvLength;

			if(bUpwindDependOnIP)
				Diag[comp] -= VecTwoNorm(CurrentIPVel[ip]) / ConvLength *
								UpwindScalar[ip];

		//////////////////////////////
		//	assemble right-hand side
		//////////////////////////////

		//	Source
		//	todo: Implement source/model term
			rhs[comp] = 0.0;

		//	Time
			if(bTimeDependent)
				rhs[comp] += IPVelOld[ip][d] / dt;

		//	Diffusion part
			rhs[comp] += kinematicViscosity * CurrentIPVel[ip][d] / DiffLengthSq;

		//	Convective part
			rhs[comp] += IPVelUpwind[ip][d] * VecTwoNorm(CurrentIPVel[ip]) / ConvLength;

		//	Pressure term
			rhs[comp] -= IPPressureGrad[ip][d];
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
			IPVelStab[ip][d] = rhs[comp] / Diag[comp];
		}
	}

//	we're done
	return true;
}

template <typename TFVGeometry>
bool GetFieldsStabilizedVelocitiesFullMatrix(	MathVector<TFVGeometry::world_dim> IPVelStab[],
									number IPStabVelShape[][TFVGeometry::world_dim],
									number IPStabPressureShape[][TFVGeometry::world_dim],
									const TFVGeometry& geo,
									MathVector<TFVGeometry::world_dim> CurrentIPVel[],
									bool bTimeDependent,
									MathVector<TFVGeometry::world_dim> IPVelOld[],
									number dt,
									MathVector<TFVGeometry::world_dim> IPVelUpwind[],
									bool bUpwindDependOnIP,
									number UpwindScalar[],
									MathVector<TFVGeometry::world_dim> CornerVel[],
									MathVector<TFVGeometry::world_dim> IPPressureGrad[],
									number kinematicViscosity)
{
	return false;
}

template <typename TFVGeometry>
bool GetFieldsStabilizedVelocities(	MathVector<TFVGeometry::world_dim> IPVelStab[],
									number IPStabVelShape[][TFVGeometry::world_dim],
									number IPStabPressureShape[][TFVGeometry::world_dim],
									const TFVGeometry& geo,
									MathVector<TFVGeometry::world_dim> CurrentIPVel[],
									bool bTimeDependent,
									MathVector<TFVGeometry::world_dim> IPVelOld[],
									number dt,
									MathVector<TFVGeometry::world_dim> IPVelUpwind[],
									bool bUpwindDependOnIP,
									number UpwindScalar[],
									MathVector<TFVGeometry::world_dim> CornerVel[],
									MathVector<TFVGeometry::world_dim> IPPressureGrad[],
									number kinematicViscosity)
{
	if(bUpwindDependOnIP)
	return
		GetFieldsStabilizedVelocitiesDiagonal(	IPVelStab,
												IPStabVelShape, IPStabPressureShape,
												geo, CurrentIPVel,
												bTimeDependent, IPVelOld, dt,
												IPVelUpwind, bUpwindDependOnIP,
												UpwindScalar, CornerVel,
												IPPressureGrad, kinematicViscosity);
	else
	return
		GetFieldsStabilizedVelocitiesFullMatrix(IPVelStab,
												IPStabVelShape, IPStabPressureShape,
												geo, CurrentIPVel,
												bTimeDependent, IPVelOld, dt,
												IPVelUpwind, bUpwindDependOnIP,
												UpwindScalar, CornerVel,
												IPPressureGrad, kinematicViscosity);
}




} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__NAVIER_STOKES__FV__STABILIZATION__ */
