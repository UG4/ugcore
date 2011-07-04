/*
 * convection_diffusion_fe1.cpp
 *
 *  Created on: 02.08.2010
 *      Author: andreasvogel
 */

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


template<typename TDomain>
template<typename TElem >
inline
bool
FE1ConvectionDiffusionElemDisc<TDomain>::
prepare_element_loop()
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
FE1ConvectionDiffusionElemDisc<TDomain>::
finish_element_loop()
{
	// all this will be performed outside of the loop over the elements.
	// Therefore it is not time critical.

	return true;
}

template<typename TDomain>
bool
FE1ConvectionDiffusionElemDisc<TDomain>::
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


template<typename TDomain>

template<typename TElem >
bool
FE1ConvectionDiffusionElemDisc<TDomain>::
prepare_element(TElem* elem, const local_vector_type& u, const local_index_type& glob_ind)
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

template<typename TDomain>
template<typename TElem >
inline
bool
FE1ConvectionDiffusionElemDisc<TDomain>::
assemble_JA(local_matrix_type& J, const local_vector_type& u)
{
	FEGeometry<TElem, dim, LagrangeP1, 1, GaussQuadrature, 2>& geo
		= GeomProvider::get<FEGeometry<TElem, dim, LagrangeP1, 1, GaussQuadrature, 2> >();

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


template<typename TDomain>
template<typename TElem >
inline
bool
FE1ConvectionDiffusionElemDisc<TDomain>::
assemble_JM(local_matrix_type& J, const local_vector_type& u)
{
	FEGeometry<TElem, dim, LagrangeP1, 1, GaussQuadrature, 2>& geo
		= GeomProvider::get<FEGeometry<TElem, dim, LagrangeP1, 1, GaussQuadrature, 2> >();

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


template<typename TDomain>
template<typename TElem >
inline
bool
FE1ConvectionDiffusionElemDisc<TDomain>::
assemble_A(local_vector_type& d, const local_vector_type& u)
{
	static const int dim = TDomain::dim;

	number integrand, shape_u;
	MathMatrix<dim,dim> D;
	MathVector<dim> v, Dgrad_u, grad_u;

	FEGeometry<TElem, dim, LagrangeP1, 1, GaussQuadrature, 2>& geo
		= GeomProvider::get<FEGeometry<TElem, dim, LagrangeP1, 1, GaussQuadrature, 2> >();

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


template<typename TDomain>
template<typename TElem >
inline
bool
FE1ConvectionDiffusionElemDisc<TDomain>::
assemble_M(local_vector_type& d, const local_vector_type& u)
{
	FEGeometry<TElem, dim, LagrangeP1, 1, GaussQuadrature, 2>& geo
		= GeomProvider::get<FEGeometry<TElem, dim, LagrangeP1, 1, GaussQuadrature, 2> >();

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


template<typename TDomain>
template<typename TElem >
inline
bool
FE1ConvectionDiffusionElemDisc<TDomain>::
assemble_f(local_vector_type& d)
{
	FEGeometry<TElem, dim, LagrangeP1, 1, GaussQuadrature, 2>& geo
		= GeomProvider::get<FEGeometry<TElem, dim, LagrangeP1, 1, GaussQuadrature, 2> >();

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

////////////////////////////////////////////////////////////////////////////////
//	Constructor
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
FE1ConvectionDiffusionElemDisc<TDomain>::
FE1ConvectionDiffusionElemDisc()
{
//	register assemling functions
	register_all_fe1_funcs();

//	register imports
	register_import(m_imDiffusion);
	register_import(m_imVelocity);
	register_import(m_imReaction);
	register_import(m_imSource);
	register_import(m_imMassScale);
}

////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////

// register for 1D
template<typename TDomain>
void
FE1ConvectionDiffusionElemDisc<TDomain>::
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
FE1ConvectionDiffusionElemDisc<TDomain>::
register_fe1_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef this_type T;

	reg_prepare_vol_loop_fct(id, &T::template prepare_element_loop<TElem>);
	reg_prepare_vol_fct(	 id, &T::template prepare_element<TElem>);
	reg_finish_vol_loop_fct( id, &T::template finish_element_loop<TElem>);
	reg_ass_JA_vol_fct(		 id, &T::template assemble_JA<TElem>);
	reg_ass_JM_vol_fct(		 id, &T::template assemble_JM<TElem>);
	reg_ass_dA_vol_fct(		 id, &T::template assemble_A<TElem>);
	reg_ass_dM_vol_fct(		 id, &T::template assemble_M<TElem>);
	reg_ass_rhs_vol_fct(	 id, &T::template assemble_f<TElem>);
}

////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

template class FE1ConvectionDiffusionElemDisc<Domain1d>;
template class FE1ConvectionDiffusionElemDisc<Domain2d>;
template class FE1ConvectionDiffusionElemDisc<Domain3d>;


} // namespace ug

