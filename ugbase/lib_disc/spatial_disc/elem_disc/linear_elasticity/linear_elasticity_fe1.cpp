/*
 * linear_elasticity_impl.h
 *
 *  Created on: 05.08.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__LINEAR_ELASTICITY__FE1_LINEAR_ELASTICITY_IMPL__
#define __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__LINEAR_ELASTICITY__FE1_LINEAR_ELASTICITY_IMPL__

#include "fe1_linear_elasticity.h"

#include "lib_disc/spatial_disc/disc_util/finite_element_geometry.h"
#include "common/util/provider.h"
#include "lib_disc/local_finite_element/lagrange/lagrange.h"
#include "lib_disc/quadrature/gauss_quad/gauss_quad.h"

namespace ug{


////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
//  Provide a generic implementation for all elements
//  (since this discretization can be implemented in a generic way)
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////


template<typename TDomain>
FE1LinearElasticityElemDisc<TDomain>::
FE1LinearElasticityElemDisc(const char* functions, const char* subsets)
	: 	IDomainElemDisc<TDomain>(TDomain::dim, functions,subsets), m_ElasticityTensorFct(NULL)
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
prepare_element(TElem* elem, const LocalVector& u)
{
//	get corners
	m_corners = this->template get_element_corners<TElem>(elem);

	typedef typename reference_element_traits<TElem>::reference_element_type
			ref_elem_type;

// 	update Geometry for this element
	FEGeometry<TElem, dim, LagrangeLSFS<ref_elem_type, 1>, GaussQuadrature<ref_elem_type, 2> >& geo
		= Provider<FEGeometry<TElem, dim, LagrangeLSFS<ref_elem_type, 1>, GaussQuadrature<ref_elem_type, 2> > >::get();

	if(!geo.update(elem, m_corners, LFEID(LFEID::LAGRANGE,1), 2))
		{UG_LOG("FE1LinearElasticityElemDisc::prepare_element:"
				" Cannot update Finite Element Geometry.\n"); return false;}

	return true;
}

template<typename TDomain>
template<typename TElem >
inline
bool
FE1LinearElasticityElemDisc<TDomain>::
assemble_JA(LocalMatrix& J, const LocalVector& u)
{
	typedef typename reference_element_traits<TElem>::reference_element_type
			ref_elem_type;

	FEGeometry<TElem, dim, LagrangeLSFS<ref_elem_type, 1>, GaussQuadrature<ref_elem_type, 2> >& geo
		= Provider<FEGeometry<TElem, dim, LagrangeLSFS<ref_elem_type, 1>, GaussQuadrature<ref_elem_type, 2> > >::get();

	for(size_t i = 0; i < geo.num_sh(); ++i) // loop corner
	{
		for(size_t c1 = 0; c1 < (size_t)dim; ++c1) // loop component
		{
			for(size_t j = 0; j < geo.num_sh(); ++j) // loop corner
			{
				for(size_t c2 = 0; c2 < (size_t)dim; ++c2) // loop component
				{
					// Compute entry A_{i, c1, j, c2}
					number integrand = 0;
					for(size_t ip = 0; ip < geo.num_ip(); ++ip) // loop ip
					{
						for(size_t d1 = 0; d1 < (size_t)dim; ++d1) // loop dimension
						{
							for(size_t d2 = 0; d2 < (size_t)dim; ++d2) // loop dimension
							{
								integrand += geo.global_grad(ip, i)[d1] * m_ElasticityTensor[c1][d1][c2][d2] * geo.global_grad(ip, j)[d2];
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
assemble_JM(LocalMatrix& J, const LocalVector& u)
{
	typedef typename reference_element_traits<TElem>::reference_element_type
			ref_elem_type;

	FEGeometry<TElem, dim, LagrangeLSFS<ref_elem_type, 1>, GaussQuadrature<ref_elem_type, 2> >& geo
		= Provider<FEGeometry<TElem, dim, LagrangeLSFS<ref_elem_type, 1>, GaussQuadrature<ref_elem_type, 2> > >::get();

	for(size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
		for(size_t i = 0; i < geo.num_sh(); ++i)
		{
			for(size_t j= 0; j < geo.num_sh(); ++j)
			{
				// same value for all components
				number value = geo.shape(ip, i) *geo.shape(ip, j) * geo.weight(ip);
				for(size_t c = 0; c < (size_t)dim; ++c)
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
assemble_A(LocalVector& d, const LocalVector& u)
{
	// Not implemented
	return false;
}


template<typename TDomain>
template<typename TElem >
inline
bool
FE1LinearElasticityElemDisc<TDomain>::
assemble_M(LocalVector& d, const LocalVector& u)
{
	// Not implemented
	return false;
}


template<typename TDomain>
template<typename TElem >
inline
bool
FE1LinearElasticityElemDisc<TDomain>::
assemble_f(LocalVector& d)
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
	typedef typename domain_traits<dim>::DimElemList ElemList;

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

	set_prep_elem_loop_fct(id, &T::template prepare_element_loop<TElem>);
	set_prep_elem_fct(	 id, &T::template prepare_element<TElem>);
	set_fsh_elem_loop_fct( id, &T::template finish_element_loop<TElem>);
	set_ass_JA_elem_fct(		 id, &T::template assemble_JA<TElem>);
	set_ass_JM_elem_fct(		 id, &T::template assemble_JM<TElem>);
	set_ass_dA_elem_fct(		 id, &T::template assemble_A<TElem>);
	set_ass_dM_elem_fct(		 id, &T::template assemble_M<TElem>);
	set_ass_rhs_elem_fct(	 id, &T::template assemble_f<TElem>);
}

////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

template class FE1LinearElasticityElemDisc<Domain1d>;
template class FE1LinearElasticityElemDisc<Domain2d>;
template class FE1LinearElasticityElemDisc<Domain3d>;


} // namespace ug


#endif /*__H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__LINEAR_ELASTICITY__FE1_LINEAR_ELASTICITY_IMPL__*/
