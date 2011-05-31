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
#include "lib_discretization/spatial_discretization/disc_helper/geometry_provider.h"
#include "lib_discretization/local_shape_function_set/lagrange/lagrange.h"
#include "lib_discretization/quadrature/gauss_quad/gauss_quad.h"

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

//	set local positions for rhs
	FEGeometry<TElem, dim, LagrangeLSFS, 1, GaussQuadrature, 2>& geo
		= GeomProvider::get<FEGeometry<TElem, dim, LagrangeLSFS, 1, GaussQuadrature, 2> >();

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

template<typename TDomain, typename TAlgebra>
template<typename TElem >
bool
FE1ConvectionDiffusionElemDisc<TDomain, TAlgebra>::
finish_element_loop()
{
	// all this will be performed outside of the loop over the elements.
	// Therefore it is not time critical.

	return true;
}

template<typename TDomain, typename TAlgebra>
bool
FE1ConvectionDiffusionElemDisc<TDomain, TAlgebra>::
time_point_changed(number time)
{
//	set new time point at imports
	m_imDiffusion.set_time(time);
	m_imVelocity.set_time(time);
	m_imSource.set_time(time);
	m_imReaction.set_time(time);
	m_imMassScale.set_time(time);

//	this disc does not need the old time solutions, thus, return false
	return false;
}


template<typename TDomain, typename TAlgebra>

template<typename TElem >
bool
FE1ConvectionDiffusionElemDisc<TDomain, TAlgebra>::
prepare_element(TElem* elem, const local_vector_type& u, const local_index_type& glob_ind)
{
	// this loop will be performed inside the loop over the elements.
	// Therefore, it is TIME CRITICAL
	typedef typename reference_element_traits<TElem>::reference_element_type
																ref_elem_type;

//	get corners
	m_vCornerCoords = this->template get_element_corners<TElem>(elem);

	// update Geometry for this element
	FEGeometry<TElem, dim, LagrangeLSFS, 1, GaussQuadrature, 2>& geo
		= GeomProvider::get<FEGeometry<TElem, dim, LagrangeLSFS, 1, GaussQuadrature, 2> >();

	if(!geo.update(&m_vCornerCoords[0]))
	{
		UG_LOG("FE1ConvectionDiffusionElemDisc::prepare_element:"
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

template<typename TDomain, typename TAlgebra>
template<typename TElem >
inline
bool
FE1ConvectionDiffusionElemDisc<TDomain, TAlgebra>::
assemble_JA(local_matrix_type& J, const local_vector_type& u)
{
	FEGeometry<TElem, dim, LagrangeLSFS, 1, GaussQuadrature, 2>& geo
		= GeomProvider::get<FEGeometry<TElem, dim, LagrangeLSFS, 1, GaussQuadrature, 2> >();

	MathVector<dim> v, Dgrad;

	for(size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
		for(size_t j = 0; j < geo.num_sh(); ++j)
		{
			// diffusion
			if(m_imDiffusion.data_given())
				MatVecMult(Dgrad, m_imDiffusion[ip], geo.grad_global(ip, j));
			else
				VecSet(Dgrad, 0.0);

			// convection
			if(m_imVelocity.data_given())
				VecScaleAppend(Dgrad, -1*geo.shape(ip,j), m_imVelocity[ip]);

			for(size_t i = 0; i < geo.num_sh(); ++i)
			{
				number integrand = VecDot(Dgrad, geo.grad_global(ip, i));

				// reaction
				if(m_imReaction.data_given())
					integrand += m_imReaction[ip] * geo.shape(ip, j) * geo.shape(ip, i);

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
assemble_JM(local_matrix_type& J, const local_vector_type& u)
{
	FEGeometry<TElem, dim, LagrangeLSFS, 1, GaussQuadrature, 2>& geo
		= GeomProvider::get<FEGeometry<TElem, dim, LagrangeLSFS, 1, GaussQuadrature, 2> >();

	for(size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
		for(size_t i = 0; i < geo.num_sh(); ++i)
		{
			for(size_t j= 0; j < geo.num_sh(); ++j)
			{
				number val = geo.shape(ip, i) *geo.shape(ip, j) * geo.weight(ip);

				if(m_imMassScale.data_given())
					val *= m_imMassScale[ip++];

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
assemble_A(local_vector_type& d, const local_vector_type& u)
{
	static const int dim = TDomain::dim;

	number integrand, shape_u;
	MathMatrix<dim,dim> D;
	MathVector<dim> v, Dgrad_u, grad_u;

	FEGeometry<TElem, dim, LagrangeLSFS, 1, GaussQuadrature, 2>& geo
		= GeomProvider::get<FEGeometry<TElem, dim, LagrangeLSFS, 1, GaussQuadrature, 2> >();

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
		if(m_imDiffusion.data_given())
			MatVecMult(Dgrad_u, m_imDiffusion[ip], grad_u);
		else
			VecSet(Dgrad_u, 0.0);

		// convection
		if(m_imReaction.data_given())
			VecScaleAppend(Dgrad_u, -1*shape_u, m_imVelocity[ip]);

		for(size_t i = 0; i < geo.num_sh(); ++i)
		{
			integrand = VecDot(Dgrad_u, geo.grad_global(ip, i));
			// reaction
			if(m_imReaction.data_given())
				integrand += m_imReaction[ip] * shape_u * geo.shape(ip, i);

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
assemble_M(local_vector_type& d, const local_vector_type& u)
{
	FEGeometry<TElem, dim, LagrangeLSFS, 1, GaussQuadrature, 2>& geo
		= GeomProvider::get<FEGeometry<TElem, dim, LagrangeLSFS, 1, GaussQuadrature, 2> >();

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
			if(m_imMassScale.data_given())
				val *= m_imMassScale[ip];

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
assemble_f(local_vector_type& d)
{
	FEGeometry<TElem, dim, LagrangeLSFS, 1, GaussQuadrature, 2>& geo
		= GeomProvider::get<FEGeometry<TElem, dim, LagrangeLSFS, 1, GaussQuadrature, 2> >();

	if(!m_imSource.data_given()) return true;

	for(size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
		for(size_t i = 0; i < geo.num_sh(); ++i)
		{
			d(_C_, i) +=  m_imSource[ip] * geo.shape(ip, i) * geo.weight(ip);
		}
	}
	return true;
}


} // namespace ug


#endif /*__H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__CONVECTION_DIFFUSION__FE1__CONVECTION_DIFFUSION_IMPL__*/
