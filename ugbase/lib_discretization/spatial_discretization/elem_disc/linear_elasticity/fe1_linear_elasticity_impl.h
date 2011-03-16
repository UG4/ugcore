/*
 * linear_elastictiy_impl.h
 *
 *  Created on: 05.08.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__LINEAR_ELASTICITY__FE1_LINEAR_ELASTICITY_IMPL__
#define __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__LINEAR_ELASTICITY__FE1_LINEAR_ELASTICITY_IMPL__

#include "fe1_linear_elasticity.h"
#include "lib_discretization/spatial_discretization/disc_helper/finite_element_geometry.h"

namespace ug{


////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
//  Provide a generic implementation for all elements
//  (since this discretization can be implemented in a generic way)
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////


template<typename TDomain, typename TAlgebra>
FE1LinearElasticityElemDisc<TDomain, TAlgebra>::
FE1LinearElasticityElemDisc(Elasticity_Tensor_fct elast)
	: 	m_ElasticityTensorFct(elast)
{
	register_assemble_functions();
};



template<typename TDomain, typename TAlgebra>
template<typename TElem >
bool
FE1LinearElasticityElemDisc<TDomain, TAlgebra>::
prepare_element_loop()
{
	// all this will be performed outside of the loop over the elements.
	// Therefore it is not time critical.

// 	evaluate Elasticity Tensor
	m_ElasticityTensorFct(m_ElasticityTensor);

	return true;
}

template<typename TDomain, typename TAlgebra>
template<typename TElem >
inline
bool
FE1LinearElasticityElemDisc<TDomain, TAlgebra>::
finish_element_loop()
{
//	nothing to do
	return true;
}

template<typename TDomain, typename TAlgebra>

template<typename TElem >
inline
bool
FE1LinearElasticityElemDisc<TDomain, TAlgebra>::
prepare_element(TElem* elem, const local_vector_type& u, const local_index_type& glob_ind)
{
//	get corners
	m_corners = this->template get_element_corners<TElem>(elem);

// 	update Geometry for this element
	FEGeometry<TElem, dim>& geo = FEGeometryProvider<TElem, dim>::get_geom(1);
	if(!geo.update(m_corners))
		{UG_LOG("FE1LinearElasticityElemDisc::prepare_element:"
				" Cannot update Finite Element Geometry.\n"); return false;}

	return true;
}

template<typename TDomain, typename TAlgebra>
template<typename TElem >
inline
bool
FE1LinearElasticityElemDisc<TDomain, TAlgebra>::
assemble_JA(local_matrix_type& J, const local_vector_type& u)
{
	FEGeometry<TElem, dim>& geo = FEGeometryProvider<TElem, dim>::get_geom(1);

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


template<typename TDomain, typename TAlgebra>
template<typename TElem >
inline
bool
FE1LinearElasticityElemDisc<TDomain, TAlgebra>::
assemble_JM(local_matrix_type& J, const local_vector_type& u)
{
	FEGeometry<TElem, dim>& geo = FEGeometryProvider<TElem, dim>::get_geom(1);

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


template<typename TDomain, typename TAlgebra>
template<typename TElem >
inline
bool
FE1LinearElasticityElemDisc<TDomain, TAlgebra>::
assemble_A(local_vector_type& d, const local_vector_type& u)
{
	// Not implemented
	return false;
}


template<typename TDomain, typename TAlgebra>
template<typename TElem >
inline
bool
FE1LinearElasticityElemDisc<TDomain, TAlgebra>::
assemble_M(local_vector_type& d, const local_vector_type& u)
{
	// Not implemented
	return false;
}


template<typename TDomain, typename TAlgebra>
template<typename TElem >
inline
bool
FE1LinearElasticityElemDisc<TDomain, TAlgebra>::
assemble_f(local_vector_type& d)
{
	// Not implemented
	return false;
}


} // namespace ug


#endif /*__H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__LINEAR_ELASTICITY__FE1_LINEAR_ELASTICITY_IMPL__*/
