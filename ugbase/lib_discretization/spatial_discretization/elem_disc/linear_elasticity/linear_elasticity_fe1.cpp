/*
 * linear_elasticity_impl.h
 *
 *  Created on: 05.08.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__LINEAR_ELASTICITY__FE1_LINEAR_ELASTICITY_IMPL__
#define __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__LINEAR_ELASTICITY__FE1_LINEAR_ELASTICITY_IMPL__

#include "fe1_linear_elasticity.h"

#include "lib_discretization/spatial_discretization/disc_util/finite_element_geometry.h"
#include "lib_discretization/spatial_discretization/disc_util/geometry_provider.h"
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
FE1LinearElasticityElemDisc<TDomain>::
FE1LinearElasticityElemDisc(Elasticity_Tensor_fct elast)
	: 	m_ElasticityTensorFct(elast)
{
	register_all_fe1_funcs();
};



template<typename TDomain>
template<typename TElem >
bool
FE1LinearElasticityElemDisc<TDomain>::
prepare_element_loop()
{
	// all this will be performed outside of the loop over the elements.
	// Therefore it is not time critical.

// 	evaluate Elasticity Tensor
	m_ElasticityTensorFct(m_ElasticityTensor);

	return true;
}

template<typename TDomain>
template<typename TElem >
inline
bool
FE1LinearElasticityElemDisc<TDomain>::
finish_element_loop()
{
//	nothing to do
	return true;
}

template<typename TDomain>

template<typename TElem >
inline
bool
FE1LinearElasticityElemDisc<TDomain>::
prepare_element(TElem* elem, const local_vector_type& u, const local_index_type& glob_ind)
{
//	get corners
	m_corners = this->template get_element_corners<TElem>(elem);

// 	update Geometry for this element
	FEGeometry<TElem, dim, LagrangeLSFS, 1, GaussQuadrature, 2>& geo
		= GeomProvider::get<FEGeometry<TElem, dim, LagrangeLSFS, 1, GaussQuadrature, 2> >();

	if(!geo.update(m_corners))
		{UG_LOG("FE1LinearElasticityElemDisc::prepare_element:"
				" Cannot update Finite Element Geometry.\n"); return false;}

	return true;
}

template<typename TDomain>
template<typename TElem >
inline
bool
FE1LinearElasticityElemDisc<TDomain>::
assemble_JA(local_matrix_type& J, const local_vector_type& u)
{
	FEGeometry<TElem, dim, LagrangeLSFS, 1, GaussQuadrature, 2>& geo
		= GeomProvider::get<FEGeometry<TElem, dim, LagrangeLSFS, 1, GaussQuadrature, 2> >();

	for(size_t i = 0; i < geo.num_sh(); ++i) // loop corner
	{
		for(size_t c1 = 0; c1 < num_fct(); ++c1) // loop component
		{
			for(size_t j = 0; j < geo.num_sh(); ++j) // loop corner
			{
				for(size_t c2 = 0; c2 < num_fct(); ++c2) // loop component
				{
					// Compute entry A_{i, c1, j, c2}
					number integrand = 0;
					for(size_t ip = 0; ip < geo.num_ip(); ++ip) // loop ip
					{
						for(size_t d1 = 0; d1 < (size_t)dim; ++d1) // loop dimension
						{
							for(size_t d2 = 0; d2 < (size_t)dim; ++d2) // loop dimension
							{
								integrand += geo.grad_global(ip, i)[d1] * m_ElasticityTensor[c1][d1][c2][d2] * geo.grad_global(ip, j)[d2];
							}
						}
						integrand *= geo.weight(ip);
					}

					J(c1, i, c2, j) += integrand;
				}
			}
		}
	}

	return true;
}


template<typename TDomain>
template<typename TElem >
inline
bool
FE1LinearElasticityElemDisc<TDomain>::
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
				// same value for all components
				number value = geo.shape(ip, i) *geo.shape(ip, j) * geo.weight(ip);
				for(size_t c = 0; c < num_fct(); ++c)
				{
					J(c, i, c, j) += value;
				}
			}
		}
	}

	return true;
}


template<typename TDomain>
template<typename TElem >
inline
bool
FE1LinearElasticityElemDisc<TDomain>::
assemble_A(local_vector_type& d, const local_vector_type& u)
{
	// Not implemented
	return false;
}


template<typename TDomain>
template<typename TElem >
inline
bool
FE1LinearElasticityElemDisc<TDomain>::
assemble_M(local_vector_type& d, const local_vector_type& u)
{
	// Not implemented
	return false;
}


template<typename TDomain>
template<typename TElem >
inline
bool
FE1LinearElasticityElemDisc<TDomain>::
assemble_f(local_vector_type& d)
{
	// Not implemented
	return false;
}

////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////

// register for 1D
template<typename TDomain>
void
FE1LinearElasticityElemDisc<TDomain>::
register_all_fe1_funcs()
{
//	get all grid element types in this dimension and below
	typedef typename GridElemTypes<dim>::DimElemList ElemList;

//	switch assemble functions
	 boost::mpl::for_each<ElemList>( RegisterFE1(this) );
}

template<typename TDomain>
template<typename TElem>
void
FE1LinearElasticityElemDisc<TDomain>::
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

template class FE1LinearElasticityElemDisc<Domain1d>;
template class FE1LinearElasticityElemDisc<Domain2d>;
template class FE1LinearElasticityElemDisc<Domain3d>;


} // namespace ug


#endif /*__H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__LINEAR_ELASTICITY__FE1_LINEAR_ELASTICITY_IMPL__*/
