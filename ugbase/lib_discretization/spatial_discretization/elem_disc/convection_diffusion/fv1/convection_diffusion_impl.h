/*
 * convection_diffusion_impl.h
 *
 *  Created on: 26.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__CONVECTION_DIFFUSION__FV1__CONVECTION_DIFFUSION_IMPL__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__CONVECTION_DIFFUSION__FV1__CONVECTION_DIFFUSION_IMPL__

#include "convection_diffusion.h"

namespace ug{


////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
//  Provide a generic implementation for all elements
//  (since this discretization can be implemented in a generic way)
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

template<	template <class TElem, int TWorldDim> class TFVGeom,
			typename TDomain,
			typename TAlgebra>
template<typename TElem >
inline
bool
FVConvectionDiffusionElemDisc<TFVGeom, TDomain, TAlgebra>::
prepare_element_loop()
{
	// all this will be performed outside of the loop over the elements.
	// Therefore it is not time critical.

	// resize corner coordinates
	typedef typename reference_element_traits<TElem>::reference_element_type
																ref_elem_type;
	m_vCornerCoords.resize(ref_elem_type::num_corners);

	// remember position attachement
	if(m_pDomain == NULL)
	{
		UG_LOG("ERROR in 'FVConvectionDiffusionElemDisc::prepare_element_loop':"
				" Domain not set.");
		return false;
	}
	m_aaPos = m_pDomain->get_position_accessor();

	return true;
}

template<	template <class TElem, int TWorldDim> class TFVGeom,
			typename TDomain,
			typename TAlgebra>
template<typename TElem >
inline
bool
FVConvectionDiffusionElemDisc<TFVGeom, TDomain, TAlgebra>::
finish_element_loop()
{
	// all this will be performed outside of the loop over the elements.
	// Therefore it is not time critical.

	return true;
}

template<	template <class TElem, int TWorldDim> class TFVGeom,
			typename TDomain,
			typename TAlgebra>
template<typename TElem >
inline
bool
FVConvectionDiffusionElemDisc<TFVGeom, TDomain, TAlgebra>::
prepare_element(TElem* elem, const local_vector_type& u,
								const local_index_type& glob_ind)
{
	// this loop will be performed inside the loop over the elements.
	// Therefore, it is TIME CRITICAL

// 	Load corners of this element
	for(size_t i = 0; i < m_vCornerCoords.size(); ++i)
	{
		VertexBase* vert = elem->vertex(i);
		m_vCornerCoords[i] = m_aaPos[vert];
	}

// 	Update Geometry for this element
	TFVGeom<TElem, dim>& geo = FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();
	if(!geo.update(elem, m_pDomain->get_subset_handler(), &m_vCornerCoords[0]))
	{
		UG_LOG("FVConvectionDiffusionElemDisc::prepare_element:"
				" Cannot update Finite Volume Geometry.\n"); return false;
	}

//	we're done
	return true;
}

template<	template <class TElem, int TWorldDim> class TFVGeom,
			typename TDomain,
			typename TAlgebra>
template<typename TElem >
inline
bool
FVConvectionDiffusionElemDisc<TFVGeom, TDomain, TAlgebra>::
assemble_JA(local_matrix_type& J, const local_vector_type& u, number time)
{
	// get finite volume geometry
	TFVGeom<TElem, dim>& geo = FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();

	// Some variables
	number flux;						// Flux at integration points
	MathMatrix<dim,dim> D;				// Diffusion Tensor
	MathVector<dim> v;					// Velocity Field
	MathVector<dim> Dgrad;				// Diff. Tensor times Gradient

	// loop Sub Control Volume Faces (SCVF)
	for(size_t i = 0; i < geo.num_scvf(); ++i)
	{
		// get current SCVF
		const typename TFVGeom<TElem, dim>::SCVF& scvf = geo.scvf(i);

		// loop integration point of SCVF
		for(size_t ip = 0; ip < scvf.num_ip(); ++ip)
		{
			// Read Data at IP
			m_Diff_Tensor(D, scvf.global_ip(ip), time);
			m_Conv_Vel(v, scvf.global_ip(ip), time);

			// loop shape functions
			for(size_t j = 0; j < scvf.num_sh(); ++j)
			{
				////////////////////////////////////////////////////
				// Diffusive Term
				// (central discretization)
				////////////////////////////////////////////////////

				// Compute Diffusion Tensor times Gradient
				MatVecMult(Dgrad, D, scvf.global_grad(j, ip));

				// Compute flux at IP
				flux = VecDot(Dgrad, scvf.normal());

				// Add flux term to local matrix
				J(_C_, scvf.from(), _C_, j) -= flux;
				J(_C_, scvf.to()  , _C_, j) += flux;

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

	////////////////////////////////////////////////////
	// Reaction Term
	// (using lumping)
	////////////////////////////////////////////////////
	// loop Sub Control Volume (SCV)
	for(size_t i = 0; i < geo.num_scv(); ++i)
	{
		// get current SCV
		const typename TFVGeom<TElem, dim>::SCV& scv = geo.scv(i);

		number val;

		// get first value
		m_Reaction(val, scv.global_ip(0), time);

		// loop other ip positions
		for(size_t ip = 1; ip < scv.num_ip(); ++ip)
		{
			number ip_val;
			m_Reaction(ip_val, scv.global_ip(ip), time);

			// TODO: Add here scaling factor for ip point
			val += ip_val;
		}

		// scale with volume of SCV
		val *= scv.volume();

		// get associated node
		const int co = scv.node_id();

		// Add to local matrix
		J(_C_, co, _C_, co) += val;
	}

	// we're done
	return true;
}


template<	template <class TElem, int TWorldDim> class TFVGeom,
			typename TDomain,
			typename TAlgebra>
template<typename TElem >
inline
bool
FVConvectionDiffusionElemDisc<TFVGeom, TDomain, TAlgebra>::
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


template<	template <class TElem, int TWorldDim> class TFVGeom,
			typename TDomain,
			typename TAlgebra>
template<typename TElem >
inline
bool
FVConvectionDiffusionElemDisc<TFVGeom, TDomain, TAlgebra>::
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


template<	template <class TElem, int TWorldDim> class TFVGeom,
			typename TDomain,
			typename TAlgebra>
template<typename TElem >
inline
bool
FVConvectionDiffusionElemDisc<TFVGeom, TDomain, TAlgebra>::
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


template<	template <class TElem, int TWorldDim> class TFVGeom,
			typename TDomain,
			typename TAlgebra>
template<typename TElem >
inline
bool
FVConvectionDiffusionElemDisc<TFVGeom, TDomain, TAlgebra>::
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


#endif /*__H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__CONVECTION_DIFFUSION__FV1__CONVECTION_DIFFUSION_IMPL__*/
