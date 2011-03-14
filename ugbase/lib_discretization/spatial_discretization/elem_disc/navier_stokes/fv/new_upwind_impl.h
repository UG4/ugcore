/*
 * new_upwind_impl.cpp
 *
 *  Created on: 10.03.2011
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__NAVIER_STOKES__NEW_STABILIZATION__IMPL__
#define __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__NAVIER_STOKES__NEW_STABILIZATION__IMPL__

// for minimum
#include <limits>

// function space, reference element
#include "lib_discretization/common/geometry_util.h"
#include "lib_discretization/local_shape_function_set/local_shape_function_set.h"
#include "lib_discretization/local_shape_function_set/local_shape_function_set_provider.h"
#include "lib_discretization/spatial_discretization/disc_helper/geometry_provider.h"

#include "new_upwind.h"
#include "lib_algebra/algebra_chooser.h"

namespace ug {

/////////////////////////////////////////////////////////////////////////////
// Interface for Upwinds
/////////////////////////////////////////////////////////////////////////////

//	register a update function for a Geometry
template <int dim, typename TAlgebra>
template <typename TFVGeom, typename TAssFunc>
void
INavierStokesUpwind<dim, TAlgebra>::
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
template <int dim, typename TAlgebra>
template <typename TFVGeom>
bool
INavierStokesUpwind<dim, TAlgebra>::
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
	TFVGeom& geo = FVGeometryProvider::get_geom<TFVGeom>();
	set_sizes(geo.num_scvf(), geo.num_scv());

//	we're done
	return true;
}

//	resize the data arrays
template <int dim, typename TAlgebra>
void
INavierStokesUpwind<dim, TAlgebra>::
set_sizes(size_t numScvf, size_t numSh)
{
//	remember sizes
	m_numScvf = numScvf;
	m_numSh = numSh;

//	adjust arrays
	m_vConvLength.resize(m_numScvf, 0);
	m_vUpShapeSh.resize(m_numScvf);
	m_vUpShapeIp.resize(m_numScvf);

	for(size_t i = 0; i < m_numScvf; ++i)
	{
		m_vUpShapeSh[i].resize(m_numSh, 0);
		m_vUpShapeIp[i].resize(m_numScvf, 0);
	}
}

/////////////////////////////////////////////////////////////////////////////
// No Upwind
/////////////////////////////////////////////////////////////////////////////

template <int TDim, typename TAlgebra>
template <typename TElem>
bool
NavierStokesNoUpwind<TDim, TAlgebra>::
update(const FV1Geometry<TElem, dim>* geo, const local_vector_type& vCornerVels)
{
//	set shapes
	for(size_t i = 0; i < geo->num_scvf(); ++i)
	{
	//	get SubControlVolumeFace
		const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(i);

		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
			upwind_shape_sh(i,sh) = scvf.shape(sh, 0);

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

template <int TDim, typename TAlgebra>
template <typename TElem>
bool
NavierStokesFullUpwind<TDim, TAlgebra>::
update(const FV1Geometry<TElem, dim>* geo, const local_vector_type& vCornerValue)
{
//	two help vectors
	MathVector<dim> dist;
    MathVector<dim> vIPVelCurrent;

// 	get corners of elem
    const MathVector<dim>* corners = geo->corners();

// 	set shapes
	for(size_t i = 0; i < geo->num_scvf(); ++i)
	{
    //	get SubControlVolumeFace
		const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(i);

    //  reset IP velocity values to zero
        VecSet(vIPVelCurrent, 0.0);
    //  reset shapes to zero for all IPs and get Velocity in IPs
        for (size_t co = 0; co < scvf.num_sh(); ++co)
        {
        // 	for all components in corners
        	upwind_shape_sh(i,co)=0.0;

        // 	Compute the Velocity in IPs
        	for(size_t d = 0; d < (size_t)dim; ++d)
        		vIPVelCurrent[d] += scvf.shape(co, 0) * vCornerValue(d, co);
        }

    // 	switch upwind
        const number flux = VecDot(scvf.normal(), vIPVelCurrent);
        if(flux > 0.0)
        {
        	upwind_shape_sh(i,scvf.from()) = 1.0;
            VecSubtract(dist, scvf.global_ip(0), corners[scvf.from()]);
        }
        else
        {
        	upwind_shape_sh(i,scvf.to()) = 1.0;
            VecSubtract(dist, scvf.global_ip(0), corners[scvf.to()]);
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
	size_t side;
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

	// TODO: Replace by SingletonProvider
	static TRefElem rRefElem;
	const int dim = TRefElem::dim;

// 	reset minimum
	number min = std::numeric_limits<number>::max();

// 	loop corners of side
	for(size_t i = 0; i < rRefElem.num_obj_of_obj(dim-1, side, 0); ++i)
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

template <int TDim, typename TAlgebra>
template <typename TElem>
bool
NavierStokesSkewedUpwind<TDim, TAlgebra>::
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
		MathVector<dim> IPVel; VecSet(IPVel, 0.0);
 		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
		{
			upwind_shape_sh(i,sh) = 0.0;
	      	for(size_t d = 0; d < (size_t)dim; ++d)
	       		IPVel[d] += scvf.shape(sh, 0) * vCornerValue(d, sh);
		}

	// 	upwind corner
		size_t sh = 0;

	// 	find upwind node
		if(!GetNodeNextToCut<typename FV1Geometry<TElem, dim>::ref_elem_type, dim>
					(sh, scvf.global_ip(0), IPVel, vCornerCoords))
		{
			UG_LOG("ERROR in GetSkewedUpwindShapes: Cannot find upwind node.\n");
			return false;
		}

	// 	set upwind corner
		upwind_shape_sh(i,sh) = 1.0;

	//	compute convection length
		MathVector<dim> dist;
	    VecSubtract(dist, scvf.global_ip(0), vCornerCoords[sh]);
        conv_length(i) = VecTwoNorm(dist);
	}

//	we're done
	return true;
}

template <int TDim, typename TAlgebra>
template <typename TElem>
bool
NavierStokesLinearProfileSkewedUpwind<TDim, TAlgebra>::
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
		MathVector<dim> IPVel; VecSet(IPVel, 0.0);
 		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
		{
			upwind_shape_sh(i,sh) = 0.0;
	      	for(size_t d = 0; d < (size_t)dim; ++d)
	       		IPVel[d] += scvf.shape(sh, 0) * vCornerValue(d, sh);
		}

 	// 	side and intersection vectors
 		static const int refDim = FV1Geometry<TElem, dim>::dim;
 		size_t side;
 		MathVector<dim> globalIntersection;
 		MathVector<refDim> localIntersection;

 	// 	find local intersection and side
 		if(!ElementSideRayIntersection<typename FV1Geometry<TElem, dim>::ref_elem_type, dim>
 			(	side, globalIntersection, localIntersection,
 				scvf.global_ip(0), IPVel, false /* search upwind */, vCornerCoords))
 		{
 			UG_LOG("ERROR in GetLinearProfileSkewedUpwindShapes: Cannot find cut side.\n");
 			return false;
 		}

 	// 	get linear trial space
 		const LocalShapeFunctionSet<typename FV1Geometry<TElem, dim>::ref_elem_type>& TrialSpace =
 				LocalShapeFunctionSetProvider::
 					get_local_shape_function_set<typename FV1Geometry<TElem, dim>::ref_elem_type>
 					(LocalShapeFunctionSetID(LocalShapeFunctionSetID::LAGRANGE, 1));

 	// 	TODO: Replace by SingletonProvider
 		static typename FV1Geometry<TElem, dim>::ref_elem_type rRefElem;

 	// 	loop corners of side
 		for(size_t j = 0; j < rRefElem.num_obj_of_obj(dim-1, side, 0); ++j)
 		{
 		// 	get corner
 			const size_t co = rRefElem.id(dim-1, side, 0, j);

 		//	evaluate trial space
 			upwind_shape_sh(j,co) = TrialSpace.shape(co, localIntersection);
 		}

 	//	compute conv length
		MathVector<dim> dist;
	    VecSubtract(dist, scvf.global_ip(0), globalIntersection);
        conv_length(i) = VecTwoNorm(dist);
	}

//	we're done
	return true;
}

} // end namespace ug

#endif /* __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__NAVIER_STOKES__NEW_STABILIZATION__IMPL__ */
