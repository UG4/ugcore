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
FE1ConvectionDiffusionElemDisc<TDomain, TAlgebra>::
FE1ConvectionDiffusionElemDisc(TDomain& domain, number upwind_amount,
							Diff_Tensor_fct diff, Conv_Vel_fct vel, Reaction_fct reac, Rhs_fct rhs)
	: 	m_domain(domain), m_upwind_amount(upwind_amount),
		m_Diff_Tensor(diff), m_Conv_Vel(vel), m_Reaction(reac), m_Rhs(rhs)
{
	register_assemble_functions();
};



template<typename TDomain, typename TAlgebra>
template<typename TElem >
inline
bool
FE1ConvectionDiffusionElemDisc<TDomain, TAlgebra>::
prepare_element_loop()
{
	// all this will be performed outside of the loop over the elements.
	// Therefore it is not time critical.

	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
	m_corners = new position_type[ref_elem_type::num_corners];

	// remember position attachement
	m_aaPos = m_domain.get_position_accessor();

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
	delete[] m_corners;

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
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;

	// load corners of this element
	for(int i = 0; i < ref_elem_type::num_corners; ++i)
	{
		VertexBase* vert = elem->vertex(i);
		m_corners[i] = m_aaPos[vert];
	}

	// update Geometry for this element
	FEGeometry<TElem, dim>& geo = FEGeometryProvider<TElem, dim>::get_geom(1);
	if(!geo.update(m_corners))
		{UG_LOG("FE1ConvectionDiffusionElemDisc::prepare_element:"
				" Cannot update Finite Element Geometry.\n"); return false;}

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

	number integrand, reac;
	MathMatrix<dim,dim> D;
	MathVector<dim> v, Dgrad;

	for(size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
		const MathVector<dim>& ipPos = geo.ip_global(ip);
		m_Diff_Tensor(D, ipPos, time);
		m_Conv_Vel(v, ipPos, time);
		m_Reaction(reac, ipPos, time);

		for(size_t j = 0; j < geo.num_sh(); ++j)
		{
			// diffusion
			MatVecMult(Dgrad, D, geo.grad_global(ip, j));

			// convection
			VecScaleAppend(Dgrad, -1*geo.shape(ip,j), v);

			for(size_t i = 0; i < geo.num_sh(); ++i)
			{
				integrand = VecDot(Dgrad, geo.grad_global(ip, i));
				// reaction
				integrand += reac * geo.shape(ip, j) * geo.shape(ip, i);
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
				J(_C_, i, _C_, j) += geo.shape(ip, i) *geo.shape(ip, j) * geo.weight(ip);
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

	number integrand, reac, shape_u;
	MathMatrix<dim,dim> D;
	MathVector<dim> v, Dgrad_u, grad_u;

	FEGeometry<TElem, dim>& geo = FEGeometryProvider<TElem, dim>::get_geom(1);

	for(size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
		const MathVector<dim>& ipPos = geo.ip_global(ip);
		m_Diff_Tensor(D, ipPos, time);
		m_Conv_Vel(v, ipPos, time);
		m_Reaction(reac, ipPos, time);

		// get current u and grad_u
		VecSet(grad_u, 0.0);
		shape_u = 0.0;
		for(size_t j = 0; j < geo.num_sh(); ++j)
		{
			VecScaleAppend(grad_u, u(_C_,j), geo.grad_global(ip, j));
			shape_u += u(_C_,j) * geo.shape(ip, j);
		}

		// diffusion
		MatVecMult(Dgrad_u, D, grad_u);

		// convection
		VecScaleAppend(Dgrad_u, -1*shape_u, v);

		for(size_t i = 0; i < geo.num_sh(); ++i)
		{
			integrand = VecDot(Dgrad_u, geo.grad_global(ip, i));
			// reaction
			integrand += reac * shape_u * geo.shape(ip, i);
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
			d(_C_, i) +=  shape_u * geo.shape(ip, i) * geo.weight(ip);
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

	number fvalue = 0.0;
	for(size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
		m_Rhs(fvalue, geo.ip_global(ip), time);

		for(size_t i = 0; i < geo.num_sh(); ++i)
		{
			d(_C_, i) +=  fvalue * geo.shape(ip, i) * geo.weight(ip);
		}
	}
	return true;
}


} // namespace ug


#endif /*__H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__CONVECTION_DIFFUSION__FE1__CONVECTION_DIFFUSION_IMPL__*/
