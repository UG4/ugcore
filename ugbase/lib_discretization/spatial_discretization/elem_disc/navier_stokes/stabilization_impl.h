/*
 *stabilization_impl.h
 *
 *  Created on: 11.03.2011
 *      Author: andreasvogel
 */

#ifndef NEW_STABILIZATION_IMPL_H___H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__NAVIER_STOKES__FV__STABILIZATION_IMPL__
#define NEW_STABILIZATION_IMPL_H___H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__NAVIER_STOKES__FV__STABILIZATION_IMPL__

#include "stabilization.h"
#include "diffusion_length.h"

#include "common/math/math_vector_matrix/math_vector_functions.h"
#include "common/math/math_vector_matrix/math_matrix_functions.h"

namespace ug{


/////////////////////////////////////////////////////////////////////////////
// Interface for Stabilization
/////////////////////////////////////////////////////////////////////////////

//	register a update function for a Geometry
template <int dim>

template <typename TFVGeom, typename TAssFunc>
void
INavierStokesStabilization<dim>::
register_update_func(TAssFunc func)
{
//	get unique geometry id
	size_t id = GetUniqueFVGeomID<TFVGeom>();

//	make sure that there is enough space
	if((size_t)id >= m_vUpdateFunc.size())
		m_vUpdateFunc.resize(id+1, NULL);

//	set pointer
	m_vUpdateFunc[id] = (UpdateFunc)func;
}

//	set the Geometry type to use for next updates
template <int dim>
template <typename TFVGeom>
bool
INavierStokesStabilization<dim>::
set_geometry_type()
{
//	get unique geometry id
	size_t id = GetUniqueFVGeomID<TFVGeom>();

//	check that function exists
	if(id >= m_vUpdateFunc.size() || m_vUpdateFunc[id] == NULL)
	{
		UG_LOG("ERROR in 'INavierStokesStabilization::set_geometry_type':"
				" No update function registered for this Geometry.\n");
		return false;
	}

//	set current geometry
	m_id = id;

//	set sizes
	TFVGeom& geo = Provider<TFVGeom>::get();
	set_sizes(geo.num_scvf(), geo.num_scv());

//	set sizes in upwind
	if(m_pUpwind != NULL)
		m_pUpwind->set_geometry_type<TFVGeom>();

//	we're done
	return true;
}


//	resize the data arrays
template <int dim>
void
INavierStokesStabilization<dim>::
set_sizes(size_t numScvf, size_t numSh)
{
//	remember sizes
	m_numScvf = numScvf;
	m_numSh = numSh;

//	adjust arrays
	m_vDiffLengthSqInv.resize(m_numScvf, 0);

	m_vStabVel.resize(m_numScvf);
	m_vvvvStabShapeVel.resize(m_numScvf);
	m_vvvvStabShapePressure.resize(m_numScvf);
	for(size_t scvf = 0; scvf < m_numScvf; ++scvf)
	{
		m_vvvvStabShapeVel[scvf].resize(dim);
		m_vvvvStabShapePressure[scvf].resize(dim);
		for(int d1 = 0; d1 < dim; ++d1)
		{
			m_vvvvStabShapeVel[scvf][d1].resize(dim);
			m_vvvvStabShapePressure[scvf][d1].resize(m_numSh, 0);
			for(int d2 = 0; d2 < dim; ++d2)
			{
				m_vvvvStabShapeVel[scvf][d1][d2].resize(m_numSh, 0);
			}
		}
	}
}

template <int dim>
bool
INavierStokesStabilization<dim>::
set_diffusion_length(std::string diffLength)
{
	if      (diffLength == "NS_RAW")        m_diffLengthType = NS_RAW;
	else if (diffLength == "NS_FIVEPOINT")  m_diffLengthType = NS_FIVEPOINT;
	else if (diffLength == "NS_COR")        m_diffLengthType = NS_COR;
	else
	{
		UG_LOG("Diffusion Length calculation method not recognized.\n");
		return false;
	}
	return true;
}

template <int dim>
template <typename TFVGeom>
bool
INavierStokesStabilization<dim>::
compute_diff_length(const TFVGeom& geo)
{
// 	Compute Diffusion Length in corresponding IPs
	switch(m_diffLengthType)
	{
		case NS_FIVEPOINT: return NSDiffLengthFivePoint(&m_vDiffLengthSqInv[0], geo);
		case NS_RAW:       return NSDiffLengthRaw(&m_vDiffLengthSqInv[0], geo);
		case NS_COR:       return NSDiffLengthCor(&m_vDiffLengthSqInv[0], geo);
        default:
        	UG_LOG("ERROR in 'INavierStokesStabilization::compute_diff_length':"
        			" Diffusion Length type defined incorrectly.\n");
        	return false;
	}
	return true;
}

template <int dim>
template <typename TFVGeom>
bool
INavierStokesStabilization<dim>::
compute_upwind(const TFVGeom& geo, const local_vector_type& vCornerValue)
{
//	check, that upwind has been set
	if(m_pUpwind == NULL)
	{
       	UG_LOG("ERROR in 'INavierStokesStabilization::compute_upwind':"
       			" No upwind method has been specified.\n");
       	return false;
 	}

//	compute upwind
	return m_pUpwind->update(geo, vCornerValue);
}

/////////////////////////////////////////////////////////////////////////////
// FIELDS
/////////////////////////////////////////////////////////////////////////////

template <int TDim>
template <typename TElem>
bool
NavierStokesFIELDSStabilization<TDim>::
update(const FV1Geometry<TElem, dim>* geo, const local_vector_type& vCornerValue,
       const DataImport<number, dim>& kinVisco,
       const DataImport<MathVector<dim>, dim>* pSource,
       const local_vector_type* pvCornerValueOldTime, number dt)
{
//	abbreviation for pressure
	static const size_t _P_ = dim;

//	Some constants
	static const size_t numIp = FV1Geometry<TElem, dim>::numSCVF;
	static const size_t numSh = FV1Geometry<TElem, dim>::numSCV;

//	compute upwind
	if(!compute_upwind(geo, vCornerValue))
	{
       	UG_LOG("ERROR in 'NavierStokesFIELDSStabilization::update':"
       			" Cannot compute upwind.\n");
       	return false;
 	}

//	compute diffusion length
	if(!compute_diff_length(*geo))
	{
       	UG_LOG("ERROR in 'NavierStokesFIELDSStabilization::update':"
       			" Cannot compute diffusion length.\n");
       	return false;
 	}

//	Find out if upwinded velocities depend on other ip velocities. In that case
//	we have to solve a matrix system. Else the system is diagonal and we can
//	compute the inverse directly

//	diagonal case (i.e. upwind vel depend only on corner vel)
	if(!non_zero_shape_ip())
	{
	//	Loop integration points
		for(size_t ip = 0; ip < numIp; ++ip)
		{
		//	get SubControlVolumeFace
			const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);

		// 	Compute current (iterated) velocity in the integration points.
			MathVector<dim> vIPVelCurrent;
			VecSet(vIPVelCurrent, 0.0);
			for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
				for(int d = 0; d < dim; ++d)
					vIPVelCurrent[d] += scvf.shape(sh) * vCornerValue(d, sh);

		// 	Loop components of velocity
			for(size_t d = 0; d < (size_t)dim; d++)
			{
			//	First, we compute the contributions to the diagonal
			//	Note: there is no contribution of the upwind vel to the diagonal
			//		  in this case, only for non-diag problems

			//	the diagonal entry
				number diag = 0.0;

			//	Time part
				if(pvCornerValueOldTime != NULL)
					diag = 1./dt;

			//	Diffusion part
				diag += kinVisco[ip] * diff_length_sq_inv(ip);

			//	Convective Term
				diag += VecTwoNorm(vIPVelCurrent) / conv_length(ip);


			//	Now, we can assemble the rhs. This rhs is assembled by all
			//	terms, that are non-dependent on the ip vel.
			//	Note, that we can compute the stab_shapes on the fly when setting
			//	up the system.

			//	Source
				number rhs = 0.0;
				if(pSource != NULL)
					rhs = (*pSource)[ip][d];

			//	Time
				if(pvCornerValueOldTime != NULL)
				{
				//	interpolate old time step
					number oldIPVel = 0.0;
					for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
						oldIPVel += scvf.shape(sh) * (*pvCornerValueOldTime)(d, sh);

				//	add to rhs
					rhs += oldIPVel / dt;
				}

			//	loop shape functions
				for(size_t k = 0; k < scvf.num_sh(); ++k)
				{
				//	Diffusion part
					number sum = kinVisco[ip] * diff_length_sq_inv(ip) * scvf.shape(k);

				//	Convection part
					sum += VecTwoNorm(vIPVelCurrent) * upwind_shape_sh(ip, k) /
																conv_length(ip);

				//	Add to rhs
					rhs += sum * vCornerValue(d, k);

				//	set stab shape
					stab_shape_vel(ip, d, d, k) = sum / diag;

				//	Pressure part
					sum = -1.0 * (scvf.global_grad(k))[d];

				//	Add to rhs
					rhs += sum * vCornerValue(_P_, k);

				//	set stab shape
					stab_shape_p(ip, d, k) = sum / diag;
				}

			//	Finally, the can invert this row
				stab_vel(ip)[d] = rhs / diag;
			}
		}
	}
	/// need to solve system
	else
	{
	//	First, we compute the current velocity at the ips
		MathVector<dim> vIPVelCurrent[numIp];

	//	Loop integration points
		for(size_t ip = 0; ip < numIp; ++ip)
		{
		//	get SubControlVolumeFace
			const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);

		// 	Compute current (iterated) velocity in the integration points.
			VecSet(vIPVelCurrent[ip], 0.0);
			for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
				for(int d = 0; d < dim; ++d)
					vIPVelCurrent[ip][d] += scvf.shape(sh) * vCornerValue(d, sh);
		}

	// 	For the FIELDS stabilization, there is no connection between the
	//	velocity components. Thus we can solve a system of size=numIP for
	//	each component of the velocity separately. This results in smaller
	//	matrices, that we have to invert.

	//	loop dimensions (i.e. components of the velocity)
		for(int d = 0; d < dim; ++d)
		{
		//	First, we have to assemble the Matrix, that includes all connections
		//	between the ip velocity component. Note, that in this case the
		//	matrix is non-diagonal and we must invert it.

		//	size of the system
			static const size_t N = numIp;

		//	a fixed size matrix
			DenseMatrix< FixedArray2<number, N, N> > mat;

		//	reset all values of the matrix to zero
			mat = 0.0;

		//	Loop integration points
			for(size_t ip = 0; ip < numIp; ++ip)
			{
			//	Time part
				if(pvCornerValueOldTime != NULL)
					mat(ip, ip) += 1./dt;

			//	Diffusion part
				mat(ip, ip) += kinVisco[ip] * diff_length_sq_inv(ip);

			//	cache this value
				const number scale = -1.0 * VecTwoNorm(vIPVelCurrent[ip]) / conv_length(ip);

			//	Convective Term (standard)
				mat(ip, ip) -= scale;

			//	Convective Term by upwind
				for(size_t ip2 = 0; ip2 < numIp; ++ip2)
					mat(ip, ip2) += upwind_shape_ip(ip, ip2) * scale;
			}

		//	we now create a matrix, where we store the inverse matrix
			typename block_traits<DenseMatrix< FixedArray2<number, N, N> > >::inverse_type inv;

		//	get the inverse
			if(!GetInverse(inv, mat))
			{
		       	UG_LOG("ERROR in 'NavierStokesFIELDSStabilization::update':"
		       			" Could not compute inverse.\n");
		       	return false;
		 	}

		//	create two vectors
			DenseVector< FixedArray1<number, N> > contVel[numSh];
			DenseVector< FixedArray1<number, N> > contP[numSh];

		//	Now, we can create several vector that describes the contribution of the
		//	corner velocities and the corner pressure. For each of this contribution
		//	components, we will apply the inverted matrix to get the stab_shapes

		//	Loop integration points
			for(size_t ip = 0; ip < numIp; ++ip)
			{
			//	get SubControlVolumeFace
				const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);

			//	loop shape functions
				for(size_t k = 0; k < numSh; ++k)
				{
				//	Diffusion part
					contVel[k][ip] = kinVisco[ip]
										* diff_length_sq_inv(ip)
										* scvf.shape(k);

				//	Convection part
					contVel[k][ip] += VecTwoNorm(vIPVelCurrent[ip])
										* upwind_shape_sh(ip, k) /
										conv_length(ip);

				//	Pressure part
					contP[k][ip] = -1.0 * (scvf.global_grad(k))[d];
				}
			}

		//	solution vector
			DenseVector< FixedArray1<number, N> > xVel, xP;

		//	compute all stab_shapes
			for(size_t k = 0; k < numSh; ++k)
			{
			//	apply for vel stab_shape
				MatMult(xVel, 1.0, inv, contVel[k]);

			//	apply for pressure stab_shape
				MatMult(xP, 1.0, inv, contP[k]);

			//	write values in data structure
				//\todo: can we optimize this, e.g. without copy?
				for(size_t ip = 0; ip < numIp; ++ip)
				{
				//	write stab_shape for vel
					stab_shape_vel(ip, d, d, k) = xVel[ip];

				//	write stab_shape for pressure
					stab_shape_p(ip, d, k) = xP[ip];
				}
			}

		//	Finally, we can compute the values of the stabilized velocity for each
		//	integration point

		//	vector of all contributions
			DenseVector< FixedArray1<number, N> > f;

		//	Loop integration points
			for(size_t ip = 0; ip < numIp; ++ip)
			{
			//	get SubControlVolumeFace
				const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(ip);

			//	Source
				f[ip] = 0.0;
				if(pSource != NULL)
					f[ip] = (*pSource)[ip][d];

			//	Time
				if(pvCornerValueOldTime != NULL)
				{
				//	interpolate old time step
				//	\todo: Is this ok? Or do we need the old stabilized vel ?
					number oldIPVel = 0.0;
					for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
						oldIPVel += scvf.shape(sh) * (*pvCornerValueOldTime)(d, sh);

				//	add to rhs
					f[ip] += oldIPVel / dt;
				}
			}

		//	sum up all contributions of vel and p for rhs.
			for(size_t k = 0; k < numSh; ++k)
			{
			//	add velocity contribution
				VecScaleAdd(f, 1.0, f, vCornerValue(d, k), contVel[k]);

			//	add pressure contribution
				VecScaleAdd(f, 1.0, f, vCornerValue(_P_, k), contVel[k]);
			}

		//	invert the system for all contributions
			DenseVector< FixedArray1<number, N> > x;
			MatMult(x, 1.0, inv, f);

		//	write values in data structure
			//\todo: can we optimize this, e.g. without copy?
			for(size_t ip = 0; ip < numIp; ++ip)
			{
			//	write stab_shape for vel
				stab_vel(ip)[d] = x[ip];
			}
		} // end dim loop

	} // end switch for non-diag

//	we're done
	return true;
}

} // end namespace ug

#endif /* NEW_STABILIZATION_IMPL_H___H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__NAVIER_STOKES__FV__STABILIZATION_IMPL__ */
