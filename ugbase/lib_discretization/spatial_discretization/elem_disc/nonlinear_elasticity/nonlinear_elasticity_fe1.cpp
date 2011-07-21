/*
 * fe1_nonlinear_elastictiy_fe1.cpp
 *
 *  Created on: 18.05.2011
 *      Author: raphaelprohl
 */

#include "fe1_nonlinear_elasticity.h"

#include "lib_discretization/spatial_discretization/disc_util/finite_element_geometry.h"
#include "common/util/provider.h"
#include "lib_discretization/local_finite_element/lagrange/lagrange.h"
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
bool
FE1NonlinearElasticityElemDisc<TDomain>::
prepare_element_loop()
{
	// all this will be performed outside of the loop over the elements.
	// Therefore it is not time critical.

// 	evaluate Elasticity Tensor
	m_ElasticityTensorFct(m_ElasticityTensor);
// 	evaluate Stress Tensor
	m_StressTensorFct(m_StressTensor);
	return true;
}

template<typename TDomain>
template<typename TElem >
inline
bool
FE1NonlinearElasticityElemDisc<TDomain>::
finish_element_loop()
{
//	nothing to do
	return true;
}

template<typename TDomain>

template<typename TElem >
inline
bool
FE1NonlinearElasticityElemDisc<TDomain>::
prepare_element(TElem* elem, const local_vector_type& u)
{
//	get corners
	m_corners = this->template get_element_corners<TElem>(elem);

// 	update Geometry for this element
	FEGeometry<TElem, dim, LagrangeLSFS, 1, GaussQuadrature, 2>& geo
			= Provider::get<FEGeometry<TElem, dim, LagrangeLSFS, 1, GaussQuadrature, 2> >();

	if(!geo.update(m_corners))
		{UG_LOG("FE1NonlinearElasticityElemDisc::prepare_element:"
				" Cannot update Finite Element Geometry.\n"); return false;}

	return true;
}

template<typename TDomain>
template<typename TElem >
inline
bool
FE1NonlinearElasticityElemDisc<TDomain>::
assemble_JA(local_matrix_type& J, const local_vector_type& u)
{
	FEGeometry<TElem, dim, LagrangeLSFS, 1, GaussQuadrature, 2>& geo
			= Provider::get<FEGeometry<TElem, dim, LagrangeLSFS, 1, GaussQuadrature, 2> >();

	//using Lagrange description of the linearized equations (cf. Bonet/Wood 1997 chapter 8/9)

	number CI[dim][dim], FI[dim][dim], FT[dim][dim], F[dim][dim], DEa[dim][dim], DEb[dim][dim];

	//TODO: mean Dilatation Term fehlt noch!
	/*for(size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
		for(size_t i = 0; i < (size_t)dim; ++i)
			{
				for(size_t j = 0; j < (size_t)dim; ++j)
				{
					FT[i][j] = F[j][i];
				}
				FT[i][i] += 1.0;
			}
		//sowas wie M_DIM_INVERT Funktion!
		vol +=  geo.weight(ip) * J1;
		Vol +=  geo.weight(ip);
	}*/

	//computing deformationgradient: Funktionsaufruf fŸr ShapeGradient
	//computing elasticityMatrix: Funktionsaufruf fŸr TangentExactSpatial

	//computing Inverse of deformationgradient
	for(size_t I = 0; I < (size_t)dim; ++I)
	{
		for(size_t j = 0; j < (size_t)dim; ++j)
		{
			FT[I][j] = F[j][I];
		}
		FT[I][I] += 1.0;
	}
	//sowas wie M_DIM_INVERT Funktion fehlt noch!

	//computing Inverse of right Cauchy-Green-Tensor
	for(size_t I = 0; I < (size_t)dim; ++I)
	{
		for(size_t J = 0; J < (size_t)dim; ++J)
		{
			number sum1 = 0.0;
			for(size_t i = 0; i < (size_t)dim; ++i)
			{
				sum1 += FI[i][I] * FI[i][J];
			}
			CI[I][J] = sum1;
		}
	}

	//computing 2.PK-stress tensor: Stress_2PK-Aufruf!

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
					number integrandC = 0;
					number integrandS = 0;
					for(size_t ip = 0; ip < geo.num_ip(); ++ip) // loop ip
					{
						for(size_t d1 = 0; d1 < (size_t)dim; ++d1) // loop dimension
						{
							for(size_t d2 = 0; d2 < (size_t)dim; ++d2) // loop dimension
							{
								number Gsyma = 0.0;
								number Gsymb = 0.0;
								DEa[d1][d2] = 0.0;
								DEb[d1][d2] = 0.0;

								Gsyma += 0.5 * (FT[d1][c1] * geo.grad_global(ip, i)[d2] + FT[d2][c1] * geo.grad_global(ip, i)[d1]);
								DEa[d1][d2] += Gsyma;
								Gsymb += 0.5 * (FT[d1][c2] * geo.grad_global(ip, j)[d2] + FT[d2][c2] * geo.grad_global(ip, j)[d1]);
								DEb[d1][d2] += Gsymb;
							}
						}

						for(size_t d1 = 0; d1 < (size_t)dim; ++d1) // loop dimension
						{
							for(size_t d2 = 0; d2 < (size_t)dim; ++d2) // loop dimension
							{
								for(size_t C1 = 0; C1 < num_fct(); ++C1) // loop component
								{
									for(size_t C2 = 0; C2 < num_fct(); ++C2) // loop component
									{
										integrandC += DEa[C1][d1] * m_ElasticityTensor[C1][d1][C2][d2] * DEb[C2][d2]; //geo.grad_global(ip, i)[d1] * m_ElasticityTensor[c1][d1][c2][d2] * geo.grad_global(ip, j)[d2];
									}
								}
								integrandS += geo.grad_global(ip, j)[d1] * m_StressTensor[d1][d2] * geo.grad_global(ip, i)[d2];
							}
						}
						integrand = geo.weight(ip) * (integrandC + integrandS);
					}

					J(c1, i, c2, j) += integrand; //richtige Reihenfolge der Komponenten?
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
FE1NonlinearElasticityElemDisc<TDomain>::
assemble_JM(local_matrix_type& J, const local_vector_type& u)
{
	FEGeometry<TElem, dim, LagrangeLSFS, 1, GaussQuadrature, 2>& geo
			= Provider::get<FEGeometry<TElem, dim, LagrangeLSFS, 1, GaussQuadrature, 2> >();

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
FE1NonlinearElasticityElemDisc<TDomain>::
assemble_A(local_vector_type& d, const local_vector_type& u)
{
	// to be implemented
	FEGeometry<TElem, dim, LagrangeLSFS, 1, GaussQuadrature, 2>& geo
			= Provider::get<FEGeometry<TElem, dim, LagrangeLSFS, 1, GaussQuadrature, 2> >();
	//TODO: mean Dilatation Term fehlt noch!

	number DE[dim][dim], FT[dim][dim], F[dim][dim];

	//computing deformationgradient: Aufruf von ShapeGradients(n,x,GIP_LOCAL(gip),mm,Gradeps,Grad,sol,F,&w);

	//computing Inverse of deformationgradient
	for(size_t I = 0; I < (size_t)dim; ++I)
	{
		for(size_t j = 0; j < (size_t)dim; ++j)
		{
			FT[I][j] = F[j][I];
		}
		FT[I][I] += 1.0;
	}
	//sowas wie M_DIM_INVERT Funktion fehlt noch!

	//Aufruf von Stress_2PK(ABS(time-time_old),F,data,T);
	for(size_t i = 0; i < geo.num_sh(); ++i) // loop corner
		{
			for(size_t c = 0; c < num_fct(); ++c) // loop component
			{
				number integrand = 0;
				for(size_t ip = 0; ip < geo.num_ip(); ++ip) // loop ip
				{
					for(size_t d1 = 0; d1 < (size_t)dim; ++d1) // loop dimension
					{
						//integrand += geo.grad_global(ip, i)[d1] * m_StressTensor[c][d1];
						for(size_t d2 = 0; d2 < (size_t)dim; ++d2) // loop dimension
						{
							number Gsym = 0.0;

							DE[d1][d2] = 0.0;
							Gsym += 0.5 * (FT[d1][c] * geo.grad_global(ip, i)[d2] + FT[d2][c] * geo.grad_global(ip, i)[d1]);
							DE[d1][d2] += Gsym;

							integrand += DE[d1][d2] * m_StressTensor[d1][d2];
						}

					}
					integrand *= geo.weight(ip);
				}

				d(i,c) -= integrand; //+= oder -=?
			}
		}
	return true; //false;
}


template<typename TDomain>
template<typename TElem >
inline
bool
FE1NonlinearElasticityElemDisc<TDomain>::
assemble_M(local_vector_type& d, const local_vector_type& u)
{
	// Not implemented
	return false;
}


template<typename TDomain>
template<typename TElem >
inline
bool
FE1NonlinearElasticityElemDisc<TDomain>::
assemble_f(local_vector_type& d)
{
	// Not implemented
	return false;
}

////////////////////////////////////////////////////////////////////////////////
//	Constructor
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
FE1NonlinearElasticityElemDisc<TDomain>::
FE1NonlinearElasticityElemDisc()
{
	register_all_fe1_funcs();
};


////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////

// register for 1D
template<typename TDomain>
void
FE1NonlinearElasticityElemDisc<TDomain>::
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
FE1NonlinearElasticityElemDisc<TDomain>::
register_fe1_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef this_type T;

	reg_prepare_elem_loop_fct(id, &T::template prepare_element_loop<TElem>);
	reg_prepare_elem_fct(	 id, &T::template prepare_element<TElem>);
	reg_finish_elem_loop_fct( id, &T::template finish_element_loop<TElem>);
	reg_ass_JA_elem_fct(		 id, &T::template assemble_JA<TElem>);
	reg_ass_JM_elem_fct(		 id, &T::template assemble_JM<TElem>);
	reg_ass_dA_elem_fct(		 id, &T::template assemble_A<TElem>);
	reg_ass_dM_elem_fct(		 id, &T::template assemble_M<TElem>);
	reg_ass_rhs_elem_fct(	 id, &T::template assemble_f<TElem>);
}

////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

template class FE1NonlinearElasticityElemDisc<Domain1d>;
template class FE1NonlinearElasticityElemDisc<Domain2d>;
template class FE1NonlinearElasticityElemDisc<Domain3d>;

} // namespace ug
