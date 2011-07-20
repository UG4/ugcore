/*
 * convection_diffusion_fe1.cpp
 *
 *  Created on: 02.08.2010
 *      Author: andreasvogel
 */

#include "convection_diffusion.h"

#include "lib_discretization/spatial_discretization/disc_util/finite_element_geometry.h"
#include "lib_discretization/spatial_discretization/disc_util/geometry_provider.h"
#include "lib_discretization/local_finite_element/lagrange/lagrange.h"
#include "lib_discretization/local_finite_element/lagrange/lagrangep1.h"
#include "lib_discretization/quadrature/gauss_quad/gauss_quad.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
//  Provide a generic implementation for all elements
//  (since this discretization can be implemented in a generic way)
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////


template<typename TDomain>
template<typename TElem >
inline
bool
ConvectionDiffusionElemDisc<TDomain>::
elem_loop_prepare_fe()
{
	// all this will be performed outside of the loop over the elements.
	// Therefore it is not time critical.

	typedef typename reference_element_traits<TElem>::reference_element_type
																ref_elem_type;
	static const int refDim = ref_elem_type::dim;

//	set local positions for rhs
	FEGeometry<TElem, dim, LagrangeP1, 1, GaussQuadrature, 2>& geo
		= GeomProvider::get<FEGeometry<TElem, dim, LagrangeP1, 1, GaussQuadrature, 2> >();

	m_imDiffusion.template 	set_local_ips<refDim>(geo.local_ips(),
											  geo.num_ip());
	m_imVelocity.template 	set_local_ips<refDim>(geo.local_ips(),
											  geo.num_ip());
	m_imSource.template 		set_local_ips<refDim>(geo.local_ips(),
											  geo.num_ip());
	m_imReaction.template set_local_ips<refDim>(geo.local_ips(),
											  geo.num_ip());
	m_imMassScale.template set_local_ips<refDim>(geo.local_ips(),
											   geo.num_ip());
	return true;
}

template<typename TDomain>
template<typename TElem >
bool
ConvectionDiffusionElemDisc<TDomain>::
elem_loop_finish_fe()
{
	// all this will be performed outside of the loop over the elements.
	// Therefore it is not time critical.

	return true;
}

template<typename TDomain>

template<typename TElem >
bool
ConvectionDiffusionElemDisc<TDomain>::
elem_prepare_fe(TElem* elem, const local_vector_type& u)
{
	// this loop will be performed inside the loop over the elements.
	// Therefore, it is TIME CRITICAL
	typedef typename reference_element_traits<TElem>::reference_element_type
																ref_elem_type;

//	get corners
	m_vCornerCoords = this->template get_element_corners<TElem>(elem);

	// update Geometry for this element
	FEGeometry<TElem, dim, LagrangeP1, 1, GaussQuadrature, 2>& geo
		= GeomProvider::get<FEGeometry<TElem, dim, LagrangeP1, 1, GaussQuadrature, 2> >();

	if(!geo.update(&m_vCornerCoords[0]))
	{
		UG_LOG("ConvectionDiffusionElemDisc::prepare_element:"
				" Cannot update Finite Element Geometry.\n");
		return false;
	}

//	set global positions for rhs
	m_imDiffusion.set_global_ips(geo.global_ips(), geo.num_ip());
	m_imVelocity.set_global_ips(geo.global_ips(), geo.num_ip());
	m_imSource.set_global_ips(geo.global_ips(), geo.num_ip());
	m_imReaction.set_global_ips(geo.global_ips(), geo.num_ip());
	m_imMassScale.set_global_ips(geo.global_ips(), geo.num_ip());

	return true;
}

template<typename TDomain>
template<typename TElem >
inline
bool
ConvectionDiffusionElemDisc<TDomain>::
elem_JA_fe(local_matrix_type& J, const local_vector_type& u)
{
	FEGeometry<TElem, dim, LagrangeP1, 1, GaussQuadrature, 2>& geo
		= GeomProvider::get<FEGeometry<TElem, dim, LagrangeP1, 1, GaussQuadrature, 2> >();

	MathVector<dim> v, Dgrad;

//	loop integration points
	for(size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
	//	loop trial space
		for(size_t j = 0; j < geo.num_sh(); ++j)
		{
		//	Diffusion
			if(m_imDiffusion.data_given())
				MatVecMult(Dgrad, m_imDiffusion[ip], geo.grad_global(ip, j));
			else
				VecSet(Dgrad, 0.0);

		//  Convection
			if(m_imVelocity.data_given())
				VecScaleAppend(Dgrad, -1*geo.shape(ip,j), m_imVelocity[ip]);

		//	loop test space
			for(size_t i = 0; i < geo.num_sh(); ++i)
			{
			//	compute integrand
				number integrand = VecDot(Dgrad, geo.grad_global(ip, i));

			// 	Reaction
				if(m_imReaction.data_given())
					integrand += m_imReaction[ip] * geo.shape(ip, j) * geo.shape(ip, i);

			//	multiply by weight
				integrand *= geo.weight(ip);

			//	add to local matrix
				J(_C_, i, _C_, j) += integrand;
			}
		}
	}

//	done
	return true;
}


template<typename TDomain>
template<typename TElem >
inline
bool
ConvectionDiffusionElemDisc<TDomain>::
elem_JM_fe(local_matrix_type& J, const local_vector_type& u)
{
	FEGeometry<TElem, dim, LagrangeP1, 1, GaussQuadrature, 2>& geo
		= GeomProvider::get<FEGeometry<TElem, dim, LagrangeP1, 1, GaussQuadrature, 2> >();

//	loop integration points
	for(size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
	//	loop test space
		for(size_t i = 0; i < geo.num_sh(); ++i)
		{
		//	loop trial space
			for(size_t j= 0; j < geo.num_sh(); ++j)
			{
			//	compute integrand
				number val = geo.shape(ip, i) *geo.shape(ip, j) * geo.weight(ip);

			//	add MassScale
				if(m_imMassScale.data_given())
					val *= m_imMassScale[ip];

			//	add to local matrix
				J(_C_, i, _C_, j) += val;
			}
		}
	}

//	done
	return true;
}


template<typename TDomain>
template<typename TElem >
inline
bool
ConvectionDiffusionElemDisc<TDomain>::
elem_dA_fe(local_vector_type& d, const local_vector_type& u)
{
	static const int dim = TDomain::dim;

	number integrand, shape_u;
	MathMatrix<dim,dim> D;
	MathVector<dim> v, Dgrad_u, grad_u;

	FEGeometry<TElem, dim, LagrangeP1, 1, GaussQuadrature, 2>& geo
		= GeomProvider::get<FEGeometry<TElem, dim, LagrangeP1, 1, GaussQuadrature, 2> >();

//	loop integration points
	for(size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
	// 	get current u and grad_u
		VecSet(grad_u, 0.0);
		shape_u = 0.0;
		for(size_t j = 0; j < geo.num_sh(); ++j)
		{
			VecScaleAppend(grad_u, u(_C_,j), geo.grad_global(ip, j));
			shape_u += u(_C_,j) * geo.shape(ip, j);
		}

	// 	Diffusion
		if(m_imDiffusion.data_given())
			MatVecMult(Dgrad_u, m_imDiffusion[ip], grad_u);
		else
			VecSet(Dgrad_u, 0.0);

	// 	Convection
		if(m_imReaction.data_given())
			VecScaleAppend(Dgrad_u, -1*shape_u, m_imVelocity[ip]);

	//	loop test spaces
		for(size_t i = 0; i < geo.num_sh(); ++i)
		{
		//	compute integrand
			integrand = VecDot(Dgrad_u, geo.grad_global(ip, i));

		// 	add Reaction
			if(m_imReaction.data_given())
				integrand += m_imReaction[ip] * shape_u * geo.shape(ip, i);

		//	multiply by integration weight
			integrand *= geo.weight(ip);

		//	add to local defect
			d(_C_, i) += integrand;
		}
	}

	return true;
}


template<typename TDomain>
template<typename TElem >
inline
bool
ConvectionDiffusionElemDisc<TDomain>::
elem_dM_fe(local_vector_type& d, const local_vector_type& u)
{
	FEGeometry<TElem, dim, LagrangeP1, 1, GaussQuadrature, 2>& geo
		= GeomProvider::get<FEGeometry<TElem, dim, LagrangeP1, 1, GaussQuadrature, 2> >();

	number shape_u;

//	loop integration points
	for(size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
	//	compute value of current solution at ip
		shape_u = 0.0;
		for(size_t j = 0; j < geo.num_sh(); ++j)
			shape_u += u(_C_,j) * geo.shape(ip, j);

	//	loop test spaces
		for(size_t i = 0; i < geo.num_sh(); ++i)
		{
		//	compute contribution
			number val = shape_u * geo.shape(ip, i) * geo.weight(ip);

		//	add MassScaling
			if(m_imMassScale.data_given())
				val *= m_imMassScale[ip];

		//	add to local defect
			d(_C_, i) +=  val;
		}
	}

//	done
	return true;
}


template<typename TDomain>
template<typename TElem >
inline
bool
ConvectionDiffusionElemDisc<TDomain>::
elem_rhs_fe(local_vector_type& d)
{
	FEGeometry<TElem, dim, LagrangeP1, 1, GaussQuadrature, 2>& geo
		= GeomProvider::get<FEGeometry<TElem, dim, LagrangeP1, 1, GaussQuadrature, 2> >();

	if(!m_imSource.data_given()) return true;

//	loop integration points
	for(size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
	//	loop test spaces
		for(size_t i = 0; i < geo.num_sh(); ++i)
		{
		//	add contribution to local defect
			d(_C_, i) +=  m_imSource[ip] * geo.shape(ip, i) * geo.weight(ip);
		}
	}

//	done
	return true;
}

////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
ConvectionDiffusionElemDisc<TDomain>::RegisterFE1::
RegisterFE1(this_type* pThis) : m_pThis(pThis){}

// register for 1D
template<typename TDomain>
void
ConvectionDiffusionElemDisc<TDomain>::
register_all_fe1_funcs()
{
//	get all grid element types in this dimension and below
	typedef typename GridElemTypes<dim>::AllElemList ElemList;

//	assemble functions
	boost::mpl::for_each<ElemList>( RegisterFE1(this) );
}

template<typename TDomain>
template<typename TElem>
void
ConvectionDiffusionElemDisc<TDomain>::
register_fe1_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef this_type T;

	reg_prepare_elem_loop_fct(id, &T::template elem_loop_prepare_fe<TElem>);
	reg_prepare_elem_fct(	 id, &T::template elem_prepare_fe<TElem>);
	reg_finish_elem_loop_fct( id, &T::template elem_loop_finish_fe<TElem>);
	reg_ass_JA_elem_fct(		 id, &T::template elem_JA_fe<TElem>);
	reg_ass_JM_elem_fct(		 id, &T::template elem_JM_fe<TElem>);
	reg_ass_dA_elem_fct(		 id, &T::template elem_dA_fe<TElem>);
	reg_ass_dM_elem_fct(		 id, &T::template elem_dM_fe<TElem>);
	reg_ass_rhs_elem_fct(	 id, &T::template elem_rhs_fe<TElem>);
}

} // namespace ug

