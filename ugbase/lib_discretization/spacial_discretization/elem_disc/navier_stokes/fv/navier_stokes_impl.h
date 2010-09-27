/*
 * navier_stokes_impl.h
 *
 *  Created on: 20.09.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__NAVIER_STOKES__FV__NAVIER_STOKES_IMPL__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__NAVIER_STOKES__FV__NAVIER_STOKES_IMPL__

#include "navier_stokes.h"

namespace ug{


////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
//  Provide a generic implementation for all elements
//  (since this discretization can be implemented in a generic way)
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////


template<template <class TElem, int TWorldDim> class TFVGeom, typename TDomain, typename TAlgebra>
FVNavierStokesElemDisc<TFVGeom, TDomain, TAlgebra>::
FVNavierStokesElemDisc(TDomain& domain, number upwind_amount,
						KinematicViscosity_fct kinVisc, Rhs_fct rhs)
	: 	m_domain(domain), m_upwindAmount(upwind_amount),
		m_kinematicViscosity(kinVisc), m_Rhs(rhs)
{
	// register all Elements with reference dimension <= world dimension
	register_assemble_functions(Int2Type<dim>());
};



template<template <class TElem, int TWorldDim> class TFVGeom, typename TDomain, typename TAlgebra>
template<typename TElem >
inline
bool
FVNavierStokesElemDisc<TFVGeom, TDomain, TAlgebra>::
prepare_element_loop()
{
	// all this will be performed outside of the loop over the elements.
	// Therefore it is not time critical.

	// resize corner coordinates
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
	m_vCornerCoords.resize(ref_elem_type::num_corners);

	// remember position attachement
	m_aaPos = m_domain.get_position_accessor();

	return true;
}

template<template <class TElem, int TWorldDim> class TFVGeom, typename TDomain, typename TAlgebra>
template<typename TElem >
inline
bool
FVNavierStokesElemDisc<TFVGeom, TDomain, TAlgebra>::
finish_element_loop()
{
	// all this will be performed outside of the loop over the elements.
	// Therefore it is not time critical.

	return true;
}

template<template <class TElem, int TWorldDim> class TFVGeom, typename TDomain, typename TAlgebra>
template<typename TElem >
inline
bool
FVNavierStokesElemDisc<TFVGeom, TDomain, TAlgebra>::
prepare_element(TElem* elem, const local_vector_type& u, const local_index_type& glob_ind)
{
	// this loop will be performed inside the loop over the elements.
	// Therefore, it is TIME CRITICAL

	// load corners of this element
	for(size_t i = 0; i < m_vCornerCoords.size(); ++i)
	{
		VertexBase* vert = elem->vertex(i);
		m_vCornerCoords[i] = m_aaPos[vert];
	}

	// update Geometry for this element
	TFVGeom<TElem, dim>& geo = FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();
	if(!geo.update(elem, m_domain.get_grid(), &m_vCornerCoords[0]))
		{UG_LOG("FVNavierStokesElemDisc::prepare_element: Cannot update Finite Volume Geometry.\n"); return false;}

	return true;
}

template<template <class TElem, int TWorldDim> class TFVGeom, typename TDomain, typename TAlgebra>
template<typename TElem >
inline
bool
FVNavierStokesElemDisc<TFVGeom, TDomain, TAlgebra>::
assemble_JA(local_matrix_type& J, const local_vector_type& u, number time)
{
	// get finite volume geometry
	TFVGeom<TElem, dim>& geo = FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();

	// user data
	number viscosity;

	// loop Sub Control Volume Faces (SCVF)
	for(size_t i = 0; i < geo.num_scvf(); ++i)
	{
		// get current SCVF
		const typename TFVGeom<TElem, dim>::SCVF& scvf = geo.scvf(i);

		// loop integration point of SCVF
		for(size_t ip = 0; ip < scvf.num_ip(); ++ip)
		{
			// Read Data at IP
			m_kinematicViscosity(viscosity, scvf.global_ip(ip), time);

			// loop shape functions
			for(size_t co = 0; co < scvf.num_sh(); ++co)
			{
				////////////////////////////////////////////////////
				// Diffusiv Term
				////////////////////////////////////////////////////

				// Compute flux at IP
				number flux = VecDot(scvf.global_grad(co, ip), scvf.normal());

				// Add flux term to local matrix
				for(size_t vel1 = 0; vel1 < dim; ++vel1)
				{
					J(vel1, scvf.from(), vel1, co) -= viscosity * flux;
					J(vel1, scvf.to()  , vel1, co) += viscosity * flux;

					for(size_t vel2 = 0; vel2 < dim; ++vel2)
					{
						const number flux2 = viscosity * scvf.global_grad(co, ip)[vel1] * scvf.normal()[vel2];
						J(vel1, scvf.from(), vel2, co) -= flux2;
						J(vel1, scvf.to()  , vel2, co) += flux2;
					}
				}

				////////////////////////////////////////////////////
				// Pressure Term
				////////////////////////////////////////////////////
				for(size_t vel = 0; vel < dim; ++vel)
				{
					J(vel, scvf.from(), _P_, co) += flux;
					J(vel, scvf.to()  , _P_, co) -= flux;
				}

				////////////////////////////////////////////////////
				// Convective Term
				// (upwinding_amount == 1.0 -> full upwind;
				//  upwinding_amount == 0.0 -> full central disc)
				////////////////////////////////////////////////////

				// central part convection
				if(m_upwindAmount != 1.0)
				{
					// Compute flux
					flux = (1.- m_upwindAmount) * scvf.shape(j, ip) * VecDot(v, scvf.normal());

					// Add flux to local matrix
					J(_C_, scvf.from(), _C_, j) += flux;
					J(_C_, scvf.to()  , _C_, j) -= flux;

				}
			}

			// upwind part convection (does only depend on one shape function)
			if(m_upwindAmount != 0.0)
			{
				// corner for upwind switch
				int up;

				// Compute flux
				flux = m_upwindAmount * VecDot(v, scvf.normal());

				// switch upwind direction
				if(flux >= 0.0) up = scvf.from(); else up = scvf.to();

				// Add flux to local matrix
				J(_C_, scvf.from(), _C_, up) += flux;
				J(_C_, scvf.to()  , _C_, up) -= flux;
			}
		}
	}


	// we're done
	return true;
}


template<template <class TElem, int TWorldDim> class TFVGeom, typename TDomain, typename TAlgebra>
template<typename TElem >
inline
bool
FVNavierStokesElemDisc<TFVGeom, TDomain, TAlgebra>::
assemble_JM(local_matrix_type& J, const local_vector_type& u, number time)
{
	// get finite volume geometry
	TFVGeom<TElem, dim>& geo = FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();

	// loop Sub Control Volumes (SCV)
	for(size_t i = 0; i < geo.num_scv(); ++i)
	{
		// get current SCV
		const typename TFVGeom<TElem, dim>::SCV& scv = geo.scv(i);

		// get associated node
		const int co = scv.node_id();

		// Add to local matrix
		J(_C_, co, _C_, co) += scv.volume();
	}

	// we're done
	return true;
}


template<template <class TElem, int TWorldDim> class TFVGeom, typename TDomain, typename TAlgebra>
template<typename TElem >
inline
bool
FVNavierStokesElemDisc<TFVGeom, TDomain, TAlgebra>::
assemble_A(local_vector_type& d, const local_vector_type& u, number time)
{
	// get finite volume geometry
	TFVGeom<TElem, dim>& geo = FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();

	number flux;					// flux at ip
	MathVector<dim> grad_u;			// gradient of solution at ip
	number shape_u;					// solution at ip
	MathMatrix<dim,dim> D;			// Diffusion Tensor
	MathVector<dim> v;				// Velocity Field
	MathVector<dim> Dgrad_u;		// Diff.Tensor times gradient of solution

	// loop Sub Control Volume Faces (SCVF)
	for(size_t i = 0; i < geo.num_scvf(); ++i)
	{
		// get current SCVF
		const typename TFVGeom<TElem, dim>::SCVF& scvf = geo.scvf(i);

		// loop integration points
		for(size_t ip = 0; ip < scvf.num_ip(); ++ip)
		{
			// reset values
			VecSet(grad_u, 0.0); shape_u = 0.0;

			// compute gradient and shape at ip
			for(size_t j = 0; j < scvf.num_sh(); ++j)
			{
				VecScaleAppend(grad_u, u(_C_,j), scvf.global_grad(j, ip));
				shape_u += u(_C_,j) * scvf.shape(j, ip);
			}

			// Evaluate user data
			m_Diff_Tensor(D, scvf.global_ip(ip), time);
			m_Conv_Vel(v, scvf.global_ip(ip), time);

			/////////////////////////////////////////////////////
			// Diffusive Term
			// (central discretization)
			/////////////////////////////////////////////////////
			MatVecMult(Dgrad_u, D, grad_u);

			// Compute flux
			flux = VecDot(Dgrad_u, scvf.normal());

			// Add to local matrix
			d(_C_, scvf.from()) -= flux;
			d(_C_, scvf.to()  ) += flux;

			/////////////////////////////////////////////////////
			// Convective Term
			// (upwinding_amount == 1.0 -> full upwind;
			//  upwinding_amount == 0.0 -> full central disc)
			/////////////////////////////////////////////////////

			// central part convection
			if(m_upwindAmount != 1.0)
			{
				// Compute flux at ip
				flux = (1.- m_upwindAmount) * shape_u * VecDot(v, scvf.normal());

				// Add to local matrix
				d(_C_, scvf.from()) += flux;
				d(_C_, scvf.to()  ) -= flux;
			}

			// upwind part convection
			if(m_upwindAmount != 0.0)
			{
				// compute flux at ip
				flux = m_upwindAmount * VecDot(v, scvf.normal());

				// Upwind switch
				if(flux >= 0.0) flux *= u(_C_, scvf.from()); else flux *= u(_C_, scvf.to());

				// Add to local matrix
				d(_C_, scvf.from()) += flux;
				d(_C_, scvf.to()  ) -= flux;
			}

		}
	}

	// loop Sub Control Volumes (SCV)
	for(size_t i = 0; i < geo.num_scv(); ++i)
	{
		// get current SCV
		const typename TFVGeom<TElem, dim>::SCV& scv = geo.scv(i);

		// first ip
		number val;
		m_Reaction(val, scv.global_ip(0), time);

		// loop other ip
		for(size_t ip = 1; ip < scv.num_ip(); ++ip)
		{
			number ip_val;
			m_Reaction(ip_val, scv.global_ip(ip), time);

			// TODO: Add weight for integration
			val += ip_val;
		}

		// scale with volume of SCV
		val *= scv.volume();

		// get associated node
		const int co = scv.node_id();

		// Add to local matrix
		d(_C_, co) += u(_C_, co) * val;
	}

	// we're done
	return true;
}


template<template <class TElem, int TWorldDim> class TFVGeom, typename TDomain, typename TAlgebra>
template<typename TElem >
inline
bool
FVNavierStokesElemDisc<TFVGeom, TDomain, TAlgebra>::
assemble_M(local_vector_type& d, const local_vector_type& u, number time)
{
	// get finite volume geometry
	TFVGeom<TElem, dim>& geo = FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();

	// loop Sub Control Volumes (SCV)
	for(size_t i = 0; i < geo.num_scv(); ++i)
	{
		// get current SCV
		const typename TFVGeom<TElem, dim>::SCV& scv = geo.scv(i);

		// get associated node
		const int co = scv.node_id();

		// Add to local matrix
		d(_C_, co) += u(_C_, co) * scv.volume();
	}

	// we're done
	return true;
}


template<template <class TElem, int TWorldDim> class TFVGeom, typename TDomain, typename TAlgebra>
template<typename TElem >
inline
bool
FVNavierStokesElemDisc<TFVGeom, TDomain, TAlgebra>::
assemble_f(local_vector_type& d, number time)
{
	// get finite volume geometry
	TFVGeom<TElem, dim>& geo = FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();

	// loop Sub Control Volumes (SCV)
	for(size_t i = 0; i < geo.num_scv(); ++i)
	{
		// get current SCV
		const typename TFVGeom<TElem, dim>::SCV& scv = geo.scv(i);

		// first value
		number val = 0.0;
		m_Rhs(val, scv.global_ip(0), time);

		// other values
		for(size_t ip = 1; ip < scv.num_ip(); ++ip)
		{
			number ip_val;
			m_Rhs(ip_val, scv.global_ip(ip), time);

			// TODO: add weights for integration
			val += ip_val;
		}

		// scale with volume of SCV
		val *= scv.volume();

		// get associated node
		const int co = scv.node_id();

		// Add to local matrix
		d(_C_, co) += val;
	}

	// we're done
	return true;
}

} // namespace ug


#endif /*__H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__NAVIER_STOKES__FV__NAVIER_STOKES_IMPL__*/
