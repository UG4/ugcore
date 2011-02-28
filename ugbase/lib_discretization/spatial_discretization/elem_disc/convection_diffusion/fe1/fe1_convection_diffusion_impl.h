/*
 * convection_diffusion_impl.h
 *
 *  Created on: 02.08.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__CONVECTION_DIFFUSION__FE1__CONVECTION_DIFFUSION_IMPL__
#define __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__CONVECTION_DIFFUSION__FE1__CONVECTION_DIFFUSION_IMPL__

#include "fe1_convection_diffusion.h"
#include "lib_discretization/spatial_discretization/disc_helper/finite_element_geometry.h"

namespace ug{


////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
//  Provide a generic implementation for all elements
//  (since this discretization can be implemented in a generic way)
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////


template<typename TDomain, typename TAlgebra>
template<typename TElem >
inline
bool
FE1ConvectionDiffusionElemDisc<TDomain, TAlgebra>::
prepare_element_loop()
{
	// all this will be performed outside of the loop over the elements.
	// Therefore it is not time critical.

	typedef typename reference_element_traits<TElem>::reference_element_type
																ref_elem_type;
	static const int refDim = ref_elem_type::dim;

	// resize corner coordinates
	m_vCornerCoords.resize(ref_elem_type::num_corners);

	// remember position attachement
	if(m_pDomain == NULL)
	{
		UG_LOG("ERROR in 'FVConvectionDiffusionElemDisc::prepare_element_loop':"
				" Domain not set.");
		return false;
	}
	m_aaPos = m_pDomain->get_position_accessor();

//	set local positions for rhs
	FEGeometry<TElem, dim>& geo = FEGeometryProvider<TElem, dim>::get_geom(1);
	m_Diff.template 	set_local_ips<refDim>(geo.local_ips(),
											  geo.num_ip());
	m_ConvVel.template 	set_local_ips<refDim>(geo.local_ips(),
											  geo.num_ip());
	m_Rhs.template 		set_local_ips<refDim>(geo.local_ips(),
											  geo.num_ip());
	m_Reaction.template set_local_ips<refDim>(geo.local_ips(),
											  geo.num_ip());
	m_MassScale.template set_local_ips<refDim>(geo.local_ips(),
											   geo.num_ip());
	return true;
}

template<typename TDomain, typename TAlgebra>
template<typename TElem >
inline
bool
FE1ConvectionDiffusionElemDisc<TDomain, TAlgebra>::
finish_element_loop()
{
	// all this will be performed outside of the loop over the elements.
	// Therefore it is not time critical.

	return true;
}

template<typename TDomain, typename TAlgebra>

template<typename TElem >
inline
bool
FE1ConvectionDiffusionElemDisc<TDomain, TAlgebra>::
prepare_element(TElem* elem, const local_vector_type& u, const local_index_type& glob_ind)
{
	// this loop will be performed inside the loop over the elements.
	// Therefore, it is TIME CRITICAL
	typedef typename reference_element_traits<TElem>::reference_element_type
																ref_elem_type;

// 	Load corners of this element
	for(size_t i = 0; i < m_vCornerCoords.size(); ++i)
	{
		VertexBase* vert = elem->vertex(i);
		m_vCornerCoords[i] = m_aaPos[vert];
	}

	// update Geometry for this element
	FEGeometry<TElem, dim>& geo = FEGeometryProvider<TElem, dim>::get_geom(1);
	if(!geo.update(&m_vCornerCoords[0]))
	{
		UG_LOG("FE1ConvectionDiffusionElemDisc::prepare_element:"
				" Cannot update Finite Element Geometry.\n");
		return false;
	}

//	set global positions for rhs
	m_Diff.set_global_ips(geo.global_ips(), geo.num_ip());
	m_ConvVel.set_global_ips(geo.global_ips(), geo.num_ip());
	m_Rhs.set_global_ips(geo.global_ips(), geo.num_ip());
	m_Reaction.set_global_ips(geo.global_ips(), geo.num_ip());
	m_MassScale.set_global_ips(geo.global_ips(), geo.num_ip());

	return true;
}

template<typename TDomain, typename TAlgebra>
template<typename TElem >
inline
bool
FE1ConvectionDiffusionElemDisc<TDomain, TAlgebra>::
assemble_JA(local_matrix_type& J, const local_vector_type& u, number time)
{
	FEGeometry<TElem, dim>& geo = FEGeometryProvider<TElem, dim>::get_geom(1);

	MathVector<dim> v, Dgrad;

	for(size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
		for(size_t j = 0; j < geo.num_sh(); ++j)
		{
			// diffusion
			if(m_Diff.data_given())
				MatVecMult(Dgrad, m_Diff[ip], geo.grad_global(ip, j));
			else
				VecSet(Dgrad, 0.0);

			// convection
			if(m_ConvVel.data_given())
				VecScaleAppend(Dgrad, -1*geo.shape(ip,j), m_ConvVel[ip]);

			for(size_t i = 0; i < geo.num_sh(); ++i)
			{
				number integrand = VecDot(Dgrad, geo.grad_global(ip, i));

				// reaction
				if(m_Reaction.data_given())
					integrand += m_Reaction[ip] * geo.shape(ip, j) * geo.shape(ip, i);

				integrand *= geo.weight(ip);

				J(_C_, i, _C_, j) += integrand;
			}
		}
	}

	return true;
}


template<typename TDomain, typename TAlgebra>
template<typename TElem >
inline
bool
FE1ConvectionDiffusionElemDisc<TDomain, TAlgebra>::
assemble_JM(local_matrix_type& J, const local_vector_type& u, number time)
{
	FEGeometry<TElem, dim>& geo = FEGeometryProvider<TElem, dim>::get_geom(1);

	for(size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
		for(size_t i = 0; i < geo.num_sh(); ++i)
		{
			for(size_t j= 0; j < geo.num_sh(); ++j)
			{
				number val = geo.shape(ip, i) *geo.shape(ip, j) * geo.weight(ip);

				if(m_MassScale.data_given())
					val *= m_MassScale[ip++];

				J(_C_, i, _C_, j) += val;
			}
		}
	}

	return true;
}


template<typename TDomain, typename TAlgebra>
template<typename TElem >
inline
bool
FE1ConvectionDiffusionElemDisc<TDomain, TAlgebra>::
assemble_A(local_vector_type& d, const local_vector_type& u, number time)
{
	static const int dim = TDomain::dim;

	number integrand, shape_u;
	MathMatrix<dim,dim> D;
	MathVector<dim> v, Dgrad_u, grad_u;

	FEGeometry<TElem, dim>& geo = FEGeometryProvider<TElem, dim>::get_geom(1);

	for(size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
		// get current u and grad_u
		VecSet(grad_u, 0.0);
		shape_u = 0.0;
		for(size_t j = 0; j < geo.num_sh(); ++j)
		{
			VecScaleAppend(grad_u, u(_C_,j), geo.grad_global(ip, j));
			shape_u += u(_C_,j) * geo.shape(ip, j);
		}

		// diffusion
		if(m_Diff.data_given())
			MatVecMult(Dgrad_u, m_Diff[ip], grad_u);
		else
			VecSet(Dgrad_u, 0.0);

		// convection
		if(m_Reaction.data_given())
			VecScaleAppend(Dgrad_u, -1*shape_u, m_ConvVel[ip]);

		for(size_t i = 0; i < geo.num_sh(); ++i)
		{
			integrand = VecDot(Dgrad_u, geo.grad_global(ip, i));
			// reaction
			if(m_Reaction.data_given())
				integrand += m_Reaction[ip] * shape_u * geo.shape(ip, i);

			integrand *= geo.weight(ip);

			d(_C_, i) += integrand;
		}
	}

	return true;
}


template<typename TDomain, typename TAlgebra>
template<typename TElem >
inline
bool
FE1ConvectionDiffusionElemDisc<TDomain, TAlgebra>::
assemble_M(local_vector_type& d, const local_vector_type& u, number time)
{
	FEGeometry<TElem, dim>& geo = FEGeometryProvider<TElem, dim>::get_geom(1);

	number shape_u;
	for(size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
		shape_u = 0.0;
		for(size_t j = 0; j < geo.num_sh(); ++j)
		{
			shape_u += u(_C_,j) * geo.shape(ip, j);
		}

		for(size_t i = 0; i < geo.num_sh(); ++i)
		{
			number val = shape_u * geo.shape(ip, i) * geo.weight(ip);
			if(m_MassScale.data_given())
				val *= m_MassScale[ip];

			d(_C_, i) +=  val;
		}
	}

	return true;
}


template<typename TDomain, typename TAlgebra>
template<typename TElem >
inline
bool
FE1ConvectionDiffusionElemDisc<TDomain, TAlgebra>::
assemble_f(local_vector_type& d, number time)
{
	FEGeometry<TElem, dim>& geo = FEGeometryProvider<TElem, dim>::get_geom(1);

	if(!m_Rhs.data_given()) return true;

	for(size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
		for(size_t i = 0; i < geo.num_sh(); ++i)
		{
			d(_C_, i) +=  m_Rhs[ip] * geo.shape(ip, i) * geo.weight(ip);
		}
	}
	return true;
}


} // namespace ug


#endif /*__H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__CONVECTION_DIFFUSION__FE1__CONVECTION_DIFFUSION_IMPL__*/
