/*
 * upwind_impl.cpp
 *
 *  Created on: 10.03.2011
 *      Author: andreasvogel
 */

#ifndef NEW_STABILIZATION_IMPL_H___H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__NAVIER_STOKES__FV__UPWIND_IMPL__
#define NEW_STABILIZATION_IMPL_H___H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__NAVIER_STOKES__FV__UPWIND_IMPL__

// for minimum
#include <limits>
#include <algorithm>

// function space, reference element
#include "lib_discretization/common/geometry_util.h"
#include "lib_discretization/local_finite_element/local_shape_function_set.h"
#include "lib_discretization/local_finite_element/local_shape_function_set.h"
#include "common/util/provider.h"

#include "upwind.h"
#include "lib_algebra/algebra_selector.h"

namespace ug {

/////////////////////////////////////////////////////////////////////////////
// Interface for Upwinds
/////////////////////////////////////////////////////////////////////////////

//	register a update function for a Geometry
template <int dim>
template <typename TFVGeom, typename TAssFunc>
void
INavierStokesUpwind<dim>::
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
INavierStokesUpwind<dim>::
set_geometry_type()
{
//	get unique geometry id
	size_t id = GetUniqueFVGeomID<TFVGeom>();

//	check that function exists
	if(id >= m_vUpdateFunc.size() || m_vUpdateFunc[id] == NULL)
	{
		UG_LOG("ERROR in 'INavierStokesUpwind::set_geometry_type':"
				" No update function registered for this Geometry.\n");
		return false;
	}

//	set current geometry
	m_id = id;

//	set sizes
	TFVGeom& geo = Provider<TFVGeom>::get();
	set_sizes(geo.num_scvf(), geo.num_scv());

//	we're done
	return true;
}

//	resize the data arrays
template <int dim>
void
INavierStokesUpwind<dim>::
set_sizes(size_t numScvf, size_t numSh)
{
//	remember sizes
	m_numScvf = numScvf;
	m_numSh = numSh;

//	adjust arrays
	m_vConvLength.resize(m_numScvf, 0);
	m_vIPVel.resize(m_numScvf);
	m_vUpShapeSh.resize(m_numScvf);
	m_vUpShapeIp.resize(m_numScvf);

	for(size_t i = 0; i < m_numScvf; ++i)
	{
		m_vUpShapeSh[i].resize(m_numSh, 0);
		m_vUpShapeIp[i].resize(m_numScvf, 0);
	}
}

///	upwind velocity
template <int dim>
MathVector<dim>
INavierStokesUpwind<dim>::
upwind_vel(size_t scvf) const
{
	UG_ASSERT(m_pCornerValue != NULL, "corner vals not set.");

//	reset result
	MathVector<dim> vel; VecSet(vel, 0.0);

//	add corner shapes
	for(size_t sh = 0; sh < num_sh(); ++sh)
		for(size_t d = 0; d < (size_t)dim; ++d)
			vel[d] += upwind_shape_sh(scvf, sh) * (*m_pCornerValue)(d, sh);

//	done if only depending on shapes
	if(!non_zero_shape_ip()) return vel;

//	compute ip vel
	for(size_t scvf2 = 0; scvf2 < num_scvf(); ++scvf)
		VecScaleAppend(vel, upwind_shape_ip(scvf, scvf2), ip_vel(scvf2));

//	return value
	return vel;
}

/////////////////////////////////////////////////////////////////////////////
// No Upwind
/////////////////////////////////////////////////////////////////////////////

template <int TDim>
template <typename TElem>
bool
NavierStokesNoUpwind<TDim>::
update(const FV1Geometry<TElem, dim>* geo, const local_vector_type& vCornerValue)
{
//	set shapes
	for(size_t i = 0; i < geo->num_scvf(); ++i)
	{
	//	get SubControlVolumeFace
		const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(i);

	//  reset IP velocity values to zero
	    VecSet(ip_vel(i), 0.0);

		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
		{
		//	set upwind shape
			upwind_shape_sh(i,sh) = scvf.shape(sh);

		// 	Compute the Velocity in IPs
			for(size_t d = 0; d < (size_t)dim; ++d)
				ip_vel(i)[d] += scvf.shape(sh) * vCornerValue(d, sh);
		}

	//	compute convection length
	//  \todo: (optional) A convection length is not really defined for no upwind.
	//	       but in the computation of a stabilization the term cancels, so
	//   	   we only have to ensure that the conv_lengh is non-zero
        conv_length(i) = 1.0;
	}

//	we're done
	return true;
}

/////////////////////////////////////////////////////////////////////////////
// Full Upwind
/////////////////////////////////////////////////////////////////////////////

template <int TDim>
template <typename TElem>
bool
NavierStokesFullUpwind<TDim>::
update(const FV1Geometry<TElem, dim>* geo, const local_vector_type& vCornerValue)
{
//	two help vectors
	MathVector<dim> dist;

// 	get corners of elem
    const MathVector<dim>* corners = geo->corners();

// 	set shapes
	for(size_t i = 0; i < geo->num_scvf(); ++i)
	{
    //	get SubControlVolumeFace
		const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(i);

    //  reset IP velocity values to zero
        VecSet(ip_vel(i), 0.0);
    //  reset shapes to zero for all IPs and get Velocity in IPs
        for (size_t sh = 0; sh < scvf.num_sh(); ++sh)
        {
        // 	for all components in corners
        	upwind_shape_sh(i,sh)=0.0;

        // 	Compute the Velocity in IPs
        	for(size_t d = 0; d < (size_t)dim; ++d)
        		ip_vel(i)[d] += scvf.shape(sh) * vCornerValue(d, sh);
        }

    // 	switch upwind
        const number flux = VecDot(scvf.normal(), ip_vel(i));
        if(flux > 0.0)
        {
        	upwind_shape_sh(i,scvf.from()) = 1.0;
            VecSubtract(dist, scvf.global_ip(), corners[scvf.from()]);
        }
        else
        {
        	upwind_shape_sh(i,scvf.to()) = 1.0;
            VecSubtract(dist, scvf.global_ip(), corners[scvf.to()]);
        }

     // compute convection length as distance between upwind point and
     // integration point
        conv_length(i) = VecTwoNorm(dist);
    }

//	we're done
	return true;
}

/////////////////////////////////////////////////////////////////////////////
// Skewed Upwind and Linear Profile Skewed Upwind
/////////////////////////////////////////////////////////////////////////////

/// computes the closest node to a elem side ray intersection
template <typename TRefElem, int TWorldDim>
bool GetNodeNextToCut(size_t& coOut,
                      const MathVector<TWorldDim>& IP,
                      const MathVector<TWorldDim>& IPVel,
                      const MathVector<TWorldDim>* vCornerCoords)
{
//	help variables
	size_t side = 0;
	MathVector<TWorldDim> globalIntersection;
	MathVector<TRefElem::dim> localIntersection;

//	compute intersection of ray in direction of ip velocity with elem side
//	we search the ray only in upwind direction
	if(!ElementSideRayIntersection<TRefElem, TWorldDim>
		(	side, globalIntersection, localIntersection,
			IP, IPVel, false /* i.e. search upwind */, vCornerCoords))
	{
		UG_LOG("ERROR in GetNodeNextToCut: Cannot find cut side.\n");
		return false;
	}

//	get reference element
	static const TRefElem& rRefElem = Provider<TRefElem>::get();
	const int dim = TRefElem::dim;

// 	reset minimum
	number min = std::numeric_limits<number>::max();

// 	loop corners of side
	for(size_t i = 0; i < rRefElem.num(dim-1, side, 0); ++i)
	{
	// 	get corner
		const size_t co = rRefElem.id(dim-1, side, 0, i);

	// 	Compute Distance to intersection
		number dist = VecDistanceSq(globalIntersection, vCornerCoords[co]);

	// 	if distance is smaller, choose this node
		if(dist < min)
		{
			min = dist;
			coOut = co;
		}
	}

//	we're done
	return true;
}

template <int TDim>
template <typename TElem>
bool
NavierStokesSkewedUpwind<TDim>::
update(const FV1Geometry<TElem, dim>* geo, const local_vector_type& vCornerValue)
{
// 	corners of geometry
	const MathVector<dim>* vCornerCoords = geo->corners();

//	loop all scvf
	for(size_t i = 0; i < geo->num_scvf(); ++i)
	{
    //	get SubControlVolumeFace
		const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(i);

	// 	reset shapes to zero and extract ip vel
		VecSet(ip_vel(i), 0.0);
 		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
		{
			upwind_shape_sh(i,sh) = 0.0;
	      	for(size_t d = 0; d < (size_t)dim; ++d)
	      		ip_vel(i)[d] += scvf.shape(sh) * vCornerValue(d, sh);
		}

	// 	upwind corner
		size_t sh = 0;

	// 	find upwind node
		if(!GetNodeNextToCut<typename FV1Geometry<TElem, dim>::ref_elem_type, dim>
					(sh, scvf.global_ip(), ip_vel(i), vCornerCoords))
		{
			UG_LOG("ERROR in GetSkewedUpwindShapes: Cannot find upwind node.\n");
			return false;
		}

	// 	set upwind corner
		upwind_shape_sh(i,sh) = 1.0;

	//	compute convection length
		MathVector<dim> dist;
	    VecSubtract(dist, scvf.global_ip(), vCornerCoords[sh]);
        conv_length(i) = VecTwoNorm(dist);
	}

//	we're done
	return true;
}

template <int TDim>
template <typename TElem>
bool
NavierStokesLinearProfileSkewedUpwind<TDim>::
update(const FV1Geometry<TElem, dim>* geo, const local_vector_type& vCornerValue)
{
// 	corners of geometry
	const MathVector<dim>* vCornerCoords = geo->corners();

//	loop all scvf
	for(size_t i = 0; i < geo->num_scvf(); ++i)
	{
    //	get SubControlVolumeFace
		const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(i);

	// 	reset shapes to zero and extract ip vel
		VecSet(ip_vel(i), 0.0);
 		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
		{
			upwind_shape_sh(i,sh) = 0.0;
	      	for(size_t d = 0; d < (size_t)dim; ++d)
	      		ip_vel(i)[d] += scvf.shape(sh) * vCornerValue(d, sh);
		}

 		if(VecTwoNorm(ip_vel(i)) == 0.0)
 		{
 		//	no upwind -> central differences
 			for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
 				upwind_shape_sh(i,sh) = scvf.shape(sh);

 		//	compute convection length
 		//  \todo: (optional) A convection length is not really defined for no upwind.
 		//	       but in the computation of a stabilization the term cancels, so
 		//   	   we only have to ensure that the conv_lengh is non-zero
 	        conv_length(i) = 1.0;

 	    //	next ip
 			continue;
 		}

 	// 	side and intersection vectors
 		static const int refDim = FV1Geometry<TElem, dim>::dim;
 		size_t side = 0;
 		MathVector<dim> globalIntersection;
 		MathVector<refDim> localIntersection;

 	// 	find local intersection and side
 		if(!ElementSideRayIntersection<typename FV1Geometry<TElem, dim>::ref_elem_type, dim>
 			(	side, globalIntersection, localIntersection,
 				scvf.global_ip(), ip_vel(i), false /* search upwind */, vCornerCoords))
 		{
 			UG_LOG("ERROR in GetLinearProfileSkewedUpwindShapes: Cannot find cut side.\n");
 			return false;
 		}

 	// 	get linear trial space
 		const LocalShapeFunctionSet<typename FV1Geometry<TElem, dim>::ref_elem_type>& TrialSpace =
 				LocalShapeFunctionSetProvider::
 					get<typename FV1Geometry<TElem, dim>::ref_elem_type>
 					(LFEID(LFEID::LAGRANGE, 1));

 	// 	get Reference Element
 		typedef typename FV1Geometry<TElem, dim>::ref_elem_type ref_elem_type;
 		static const ref_elem_type& rRefElem
 			= Provider<ref_elem_type>::get();

 	// 	loop corners of side
 		for(size_t j = 0; j < rRefElem.num(dim-1, side, 0); ++j)
 		{
 		// 	get corner
 			const size_t co = rRefElem.id(dim-1, side, 0, j);

 		//	evaluate trial space
 			upwind_shape_sh(j,co) = TrialSpace.shape(co, localIntersection);
 		}

 	//	compute conv length
		MathVector<dim> dist;
	    VecSubtract(dist, scvf.global_ip(), globalIntersection);
        conv_length(i) = VecTwoNorm(dist);
	}

//	we're done
	return true;
}



template <int TDim>
template <typename TElem>
bool
NavierStokesPositiveUpwind<TDim>::
update(const FV1Geometry<TElem, dim>* geo, const local_vector_type& vCornerValue)
{

//	1. Reset values and compute ip velocities and Compute mass fluxes at ip's

//	vector for flux values
	std::vector<number> vMassFlux(geo->num_scvf(), 0.0);

//	loop all scvf
	for(size_t i = 0; i < geo->num_scvf(); ++i)
	{
	//	get SubControlVolumeFace
		const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(i);

	// 	reset shapes w.r.t corner value to zero and extract ip vel
		VecSet(ip_vel(i), 0.0);
		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
		{
			upwind_shape_sh(i,sh) = 0.0;
			for(size_t d = 0; d < (size_t)dim; ++d)
				ip_vel(i)[d] += scvf.shape(sh) * vCornerValue(d, sh);
		}

	// 	reset shapes w.r.t. ip value to zero and extract ip vel
		for(size_t ip = 0; ip < geo->num_scvf(); ++ip)
			upwind_shape_ip(i,ip) = 0.0;

	//	compute flux
		vMassFlux[i] = VecProd(ip_vel(i), scvf.normal());
	}


//	2. Handle each SCV separately

	for(size_t sh = 0; sh < this->num_sh(); ++sh)
	{
	//	reset inflow, outflow
		number m_in = 0, m_out = 0;
		std::vector<size_t> vSCVIP;
		std::vector<number> vFlux;

	//	loop subcontrol volume faces
		for(size_t i = 0; i < geo->num_scvf(); ++i)
		{
		//	get SubControlVolumeFace
			const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(i);

		//	if scvf is part of the scv, add fluxes
			if(scvf.from() == sh)
			{
			//	normal directed outwards
				vSCVIP.push_back(i);
				vFlux.push_back( vMassFlux[i] );
				m_in += -1.0 * std::min(vMassFlux[i], 0.0);
				m_out += std::max(vMassFlux[i], 0.0);
			}
			else if (scvf.to() == sh)
			{
			//	normal directed inwards
				vSCVIP.push_back(i);
				vFlux.push_back( -1.0 * vMassFlux[i] );
				m_in += -1.0 * std::min(-1.0 *  vMassFlux[i], 0.0);
				m_out += std::max(-1.0 * vMassFlux[i], 0.0);
			}
		}

	//	compute F
		number F = std::max(m_in, m_out);

	//	set shapes
		for(size_t i = 0; i < vSCVIP.size(); ++i)
		{
			if(vFlux[i] > 0)
			{
				number sum = 0.0;
				for(size_t j = 0; j < vSCVIP.size(); ++j)
				{
					if(vFlux[j] < 0)
					{
					//	set ip shapes
						upwind_shape_ip(i, j) = -1.0 * vFlux[j] / F;
						sum += upwind_shape_ip(i, j);
					}
				}
			//	set nodal shapes
				upwind_shape_sh(i, sh) = 1.0 - sum;
			}
		}
	}

//	3. compute convection length

// 	corners of geometry
	const MathVector<dim>* vCornerCoords = geo->corners();

//	compute upwind point
	MathVector<dim> upPos; VecSet(upPos, 0.0);
	for(size_t i = 0; i < geo->num_scvf(); ++i)
	{
	//	get SubControlVolumeFace
		const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(i);

	//	sum up contributions
        for (size_t sh = 0; sh < scvf.num_sh(); ++sh)
        	VecScaleAppend(upPos, upwind_shape_sh(i, sh), vCornerCoords[sh]);
        for (size_t j = 0; j < geo->num_scvf(); ++j)
        	VecScaleAppend(upPos, upwind_shape_sh(i, j), geo->scvf(j).global_ip());

    //	save convection length
		MathVector<dim> dist;
	    VecSubtract(dist, scvf.global_ip(), upPos);
        conv_length(i) = VecTwoNorm(dist);
	}

//	we're done
	return true;
}

} // end namespace ug

#endif /* NEW_STABILIZATION_IMPL_H___H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__NAVIER_STOKES__FV__UPWIND_IMPL__ */
