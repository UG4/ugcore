/*
 * convection_diffusion_impl.h
 *
 *  Created on: 26.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__CONVECTION_DIFFUSION_EQUATION__CONVECTION_DIFFUSION_IMPL__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__CONVECTION_DIFFUSION_EQUATION__CONVECTION_DIFFUSION_IMPL__

#include "convection_diffusion.h"

namespace ug{


////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
//  Provide a generic implementation for all elements
//  (since this discretization can be implemented in a generic way)
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////


template<typename TDomain, typename TAlgebra>
ConvectionDiffusionElemDisc<TDomain, TAlgebra>::
ConvectionDiffusionElemDisc(TDomain& domain, number upwind_amount,
							Diff_Tensor_fct diff, Conv_Vel_fct vel, Reaction_fct reac, Rhs_fct rhs)
	: 	m_domain(domain), m_upwind_amount(upwind_amount),
		m_Diff_Tensor(diff), m_Conv_Vel(vel), m_Reaction(reac), m_Rhs(rhs)
{
	register_all_assemble_functions<Triangle>(RET_TRIANGLE);
	register_all_assemble_functions<Quadrilateral>(RET_QUADRILATERAL);
};



template<typename TDomain, typename TAlgebra>
template<typename TElem >
inline
bool
ConvectionDiffusionElemDisc<TDomain, TAlgebra>::
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
ConvectionDiffusionElemDisc<TDomain, TAlgebra>::
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
ConvectionDiffusionElemDisc<TDomain, TAlgebra>::
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
	get_fvgeom<TElem>().update(m_corners);

	return true;
}

template<typename TDomain, typename TAlgebra>
template<typename TElem >
inline
bool
ConvectionDiffusionElemDisc<TDomain, TAlgebra>::
assemble_JA(local_matrix_type& J, const local_vector_type& u, number time)
{
	static const int dim = TDomain::dim;

	number flux;
	MathMatrix<dim,dim> D;
	MathVector<dim> v, lin_Defect;
	MathVector<dim> Dgrad;
	for(size_t i = 0; i < get_fvgeom<TElem>().num_scvf(); ++i)
	{
		const SubControlVolumeFace<TElem, TDomain::dim>& scvf = get_fvgeom<TElem>().scvf(i);

		for(size_t ip = 0; ip < scvf.num_ip(); ++ip)
		{
			const SD_Values<TElem, TDomain::dim>& sdv = scvf.sdv();

			m_Diff_Tensor(D, scvf.global_ip(), time);
			m_Conv_Vel(v, scvf.global_ip(), time);

			for(size_t j = 0; j < sdv.num_sh(); ++j)
			{
				////////////////////////////////////
				// diffusiv term (central discretization)
				MatVecMult(Dgrad, D, sdv.grad_global(j));
				flux = VecDot(Dgrad, scvf.normal());

				J(scvf.from(), j) -= flux;
				J(scvf.to()  , j) += flux;

				////////////////////////////////////
				// convective term
				// (upwinding_amount == 1.0 -> full upwind;
				//  upwinding_amount == 0.0 -> full central disc)

				// central part convection
				if(m_upwind_amount != 1.0)
				{
					flux = (1.- m_upwind_amount) * sdv.shape(j) * VecDot(v, scvf.normal());

					// coupling 'from' with j  (i.e. A[from][j]) and 'to' with j (i.e. A[to][j])
					J(scvf.from(), j) += flux;
					J(scvf.to()  , j) -= flux;

				}
			}
			// upwind part convection
			if(m_upwind_amount != 0.0)
			{
				int up;
				flux = m_upwind_amount * VecDot(v, scvf.normal());
				if(flux >= 0.0) up = scvf.from(); else up = scvf.to();
				J(scvf.from(), up) += flux;
				J(scvf.to()  , up) -= flux;
			}

		}
	}
	int co;
	number reac_val;
	for(size_t i = 0; i < get_fvgeom<TElem>().num_scv(); ++i)
	{
		const SubControlVolume<TElem, TDomain::dim>& scv = get_fvgeom<TElem>().scv(i);

		co = scv.local_corner_id();

		m_Reaction(reac_val, scv.global_corner_pos(), time);

		J(co , co) += reac_val * scv.volume();
	}

	return true;
}


template<typename TDomain, typename TAlgebra>
template<typename TElem >
inline
bool
ConvectionDiffusionElemDisc<TDomain, TAlgebra>::
assemble_JM(local_matrix_type& J, const local_vector_type& u, number time)
{
	int co;
	for(size_t i = 0; i < get_fvgeom<TElem>().num_scv(); ++i)
	{
		const SubControlVolume<TElem, TDomain::dim>& scv = get_fvgeom<TElem>().scv(i);

		co = scv.local_corner_id();

		J(co , co) += 1.0 * scv.volume();
	}

	return true;
}


template<typename TDomain, typename TAlgebra>
template<typename TElem >
inline
bool
ConvectionDiffusionElemDisc<TDomain, TAlgebra>::
assemble_A(local_vector_type& d, const local_vector_type& u, number time)
{
	static const int dim = TDomain::dim;

	number flux;
	MathVector<dim> grad_u;
	number shape_u;
	MathMatrix<dim,dim> D;
	MathVector<dim> v;
	MathVector<dim> Dgrad_u;
	for(size_t i = 0; i < get_fvgeom<TElem>().num_scvf(); ++i)
	{
		const SubControlVolumeFace<TElem, TDomain::dim>& scvf = get_fvgeom<TElem>().scvf(i);

		for(size_t ip = 0; ip < scvf.num_ip(); ++ip)
		{
			const SD_Values<TElem, TDomain::dim>& sdv = scvf.sdv();

			VecSet(grad_u, 0.0);
			shape_u = 0.0;
			for(size_t j = 0; j < sdv.num_sh(); ++j)
			{
				VecScaleAppend(grad_u, u[j], sdv.grad_global(j));
				shape_u += u[j] * sdv.shape(j);
			}

			m_Diff_Tensor(D, scvf.global_ip(), time);
			m_Conv_Vel(v, scvf.global_ip(), time);

			////////////////////////////////////
			// diffusiv term (central discretization)
			MatVecMult(Dgrad_u, D, grad_u);
			flux = VecDot(Dgrad_u, scvf.normal());

			d[scvf.from()] -= flux;
			d[scvf.to()] += flux;

			////////////////////////////////////
			// convective term
			// (upwinding_amount == 1.0 -> full upwind;
			//  upwinding_amount == 0.0 -> full central disc)

			// central part convection
			if(m_upwind_amount != 1.0)
			{
				flux = (1.- m_upwind_amount) * shape_u * VecDot(v, scvf.normal());

				d[scvf.from()] += flux;
				d[scvf.to()] -= flux;
			}

			// upwind part convection
			if(m_upwind_amount != 0.0)
			{
				flux = m_upwind_amount * VecDot(v, scvf.normal());
				if(flux >= 0.0) flux *= u[scvf.from()]; else flux *= u[scvf.to()];
				d[scvf.from()] += flux;
				d[scvf.to()] -= flux;
			}

		}
	}
	int co;
	number reac_val;
	for(size_t i = 0; i < get_fvgeom<TElem>().num_scv(); ++i)
	{
		const SubControlVolume<TElem, TDomain::dim>& scv = get_fvgeom<TElem>().scv(i);

		co = scv.local_corner_id();

		m_Reaction(reac_val, scv.global_corner_pos(), time);

		d[co] += reac_val * u[co] * scv.volume();
	}

	return true;
}


template<typename TDomain, typename TAlgebra>
template<typename TElem >
inline
bool
ConvectionDiffusionElemDisc<TDomain, TAlgebra>::
assemble_M(local_vector_type& d, const local_vector_type& u, number time)
{
	int co;
	for(size_t i = 0; i < get_fvgeom<TElem>().num_scv(); ++i)
	{
		const SubControlVolume<TElem, TDomain::dim>& scv = get_fvgeom<TElem>().scv(i);

		co = scv.local_corner_id();

		d[co] += u[co] * scv.volume();
	}

	return true;
}


template<typename TDomain, typename TAlgebra>
template<typename TElem >
inline
bool
ConvectionDiffusionElemDisc<TDomain, TAlgebra>::
assemble_f(local_vector_type& d, number time)
{
	number fvalue = 0.0;
	for(size_t i = 0; i < get_fvgeom<TElem>().num_scv(); ++i)
	{
		const SubControlVolume<TElem, TDomain::dim>& scv = get_fvgeom<TElem>().scv(i);

		m_Rhs(fvalue, scv.global_corner_pos(), time);
		d[scv.local_corner_id()] += fvalue * scv.volume();
	}

	return true;
}


} // namespace ug


#endif /*__H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__CONVECTION_DIFFUSION_EQUATION__CONVECTION_DIFFUSION_IMPL__*/
