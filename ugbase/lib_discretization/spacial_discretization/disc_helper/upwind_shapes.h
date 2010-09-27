/*
 * upwind_shapes.h
 *
 *  Created on: 23.09.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__DISC_HELPER__UPWIND_SHAPES__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__DISC_HELPER__UPWIND_SHAPES__

#include <limits>

#include "lib_discretization/common/geometry_util.h"

namespace ug{

template <typename TRefElem, int TWorldDim>
bool GetNodeNextToCut(size_t& coOut, const MathVector<TWorldDim>& IP,
								const MathVector<TWorldDim>& IPVel, const MathVector<TWorldDim>* vCornerCoords)
{
	size_t side;
	MathVector<TWorldDim> globalIntersection;
	MathVector<TRefElem::dim> localIntersection;

	if(!ElementSideRayIntersection<TRefElem, TWorldDim>
		(	side, globalIntersection, localIntersection,
			IP, IPVel, false /* search upwind */, vCornerCoords))
	{
		UG_LOG("ERROR in GetNodeNextToCut: Cannot find cut side.\n");
		return false;
	}

	// TODO: Replace by SingletonProvider
	static TRefElem rRefElem;
	const int dim = TRefElem::dim;

	// reset minimum
	number min = std::numeric_limits<number>::max();

	// loop corners of side
	for(size_t i = 0; i < rRefElem.num_obj_of_obj(dim-1, side, 0); ++i)
	{
		// get corner
		const size_t co = rRefElem.num_obj_of_obj(dim-1, side, 0, i);

		// Compute Distance to intersection
		number dist = VecDistanceSq(globalIntersection, vCornerCoords[co]);

		// if smaller
		if(dist < min)
		{
			min = dist;
			coOut = co;
		}
	}

	return true;
}


template <typename TSCVF>
bool GetNoUpwindShapes(const TSCVF& scvf, const MathVector<TSCVF::world_dim>& IPVel, number* shape)
{
	// set shapes
	for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
		shape[sh] = scvf.shape(sh, 0);

	return true;
}

template <typename TSCVF>
bool GetFullUpwindShapes(const TSCVF& scvf, const MathVector<TSCVF::world_dim>& IPVel, number* shape)
{
	// reset shapes to zero
	for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
		shape[sh] = 0.0;

	// switch upwind
	const number sh = VecDot(scvf.normal(), IPVel);
	if(sh > 0.0)
		shape[scvf.from()] = 1.0;
	else
		shape[scvf.to()] = 1.0;

	return true;
}

template <typename TSCVF>
bool GetSkewedUpwindShapes(const TSCVF& scvf, const MathVector<TSCVF::world_dim>& IPVel, number* shape)
{
	// reset shapes to zero
	for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
		shape[sh] = 0.0;

	// get geometry of scvf
	typename TSCVF::geometry_type& geo = scvf.geometry();

	// corners of geometry
	const MathVector<TSCVF::world_dim>* vCornerCoords = geo.corners();

	// upwind corner
	size_t co = 0;

	// find upwind node
	if(!GetNodeNextToCut<typename TSCVF::ref_elem_type, TSCVF::world_dim>(co, scvf.global_ip(0), IPVel, vCornerCoords))
	{
		UG_LOG("ERROR in GetSkewedUpwindShapes: Cannot find upwind node.\n");
		return false;
	}

	// set upwind corner
	shape[co] = 1.0;

	return true;
}

// Linear Profile Skewed Upwind
template <typename TSCVF>
bool GetLinearProfileSkewedUpwindShapes(const TSCVF& scvf, const MathVector<TSCVF::world_dim>& IPVel, number* shape)
{
	// reset shapes to zero
	for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
		shape[sh] = 0.0;

	// get geometry of scvf
	typename TSCVF::geometry_type& geo = scvf.geometry();

	// corners of geometry
	const MathVector<TSCVF::world_dim>* vCornerCoords = geo.corners();

	// reference dimension
	const int dim = TSCVF::dim;

	// side and intersection vectors
	size_t side;
	MathVector<TSCVF::world_dim> globalIntersection;
	MathVector<dim> localIntersection;

	// find local intersection and side
	if(!ElementSideRayIntersection<TSCVF::ref_elem_type, TSCVF::world_dim>
		(	side, globalIntersection, localIntersection,
			scvf.global_ip(0), IPVel, false /* search upwind */, vCornerCoords))
	{
		UG_LOG("ERROR in GetLinearProfileSkewedUpwindShapes: Cannot find cut side.\n");
		return false;
	}

	// get linear trial space
	const LocalShapeFunctionSet<typename TSCVF::ref_elem_type>& TrialSpace =
			LocalShapeFunctionSetFactory::inst().get_local_shape_function_set<typename TSCVF::ref_elem_type>(LSFS_LAGRANGEP1);

	// TODO: Replace by SingletonProvider
	static typename TSCVF::ref_elem_type& rRefElem;

	// loop corners of side
	for(size_t i = 0; i < rRefElem.num_obj_of_obj(dim-1, side, 0); ++i)
	{
		// get corner
		const size_t co = rRefElem.num_obj_of_obj(dim-1, side, 0, i);

		if(!TrialSpace.evaluate(co, localIntersection, shape[co]))
		{
			UG_LOG("Cannot evaluate local shape.\n");
			return false;
		}
	}

	return true;
}



} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__DISC_HELPER__UPWIND_SHAPES__ */
