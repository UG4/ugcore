/*
 * convection_diffusion_assemble.cpp
 *
 *  Created on: 26.02.2010
 *      Author: andreasvogel
 */

#ifndef CONV_DIFF_IMPL_H
#define CONV_DIFF_IMPL_H

namespace ug{


////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
//  Provide a generic implementation for all elements
//  (since this discretization can be implemented in a generic way)
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

template<typename TDomain, typename TElem >
inline
void
ConvectionDiffusionEquation<TDomain, TElem>::
prepare_element_loop()
{
	// all this will be performed outside of the loop over the elements.
	// Therefore it is not time critical.

	// create new Geometry
	m_geo = new FVElementGeometry<TElem>();
	assert(m_geo != NULL);

	// remember position attachement
	m_aaPos = m_domain.get_position_accessor();

}

template<typename TDomain, typename TElem >
inline
void
ConvectionDiffusionEquation<TDomain, TElem>::
finish_element_loop()
{
	// all this will be performed outside of the loop over the elements.
	// Therefore it is not time critical.

	// delete Geometry
	if(m_geo != NULL)
		delete m_geo;
}

template<typename TDomain, typename TElem >
inline
void
ConvectionDiffusionEquation<TDomain, TElem>::
prepare_element(TElem* elem)
{
	// this loop will be performed inside the loop over the elements.
	// Therefore, it is TIME CRITICAL

	// load corners of this element
	for(int i = 0; i < reference_element_traits<TElem>::num_corners; ++i)
	{
		VertexBase* vert = elem->vertex(i);
		m_corners[i] = m_aaPos[vert];
	}

	// update Geometry for this element
	m_geo->update(m_corners);
}

template<typename TDomain, typename TElem >
inline
void
ConvectionDiffusionEquation<TDomain, TElem>::
assemble_element_JA(local_matrix_type& J, const local_vector_type& u, number time)
{
	static const int dim = TDomain::dim;

	number flux;
	MathMatrix<dim,dim> D;
	MathVector<dim> v;
	MathVector<dim> Dgrad;
	for(uint i = 0; i < m_geo->num_scvf(); ++i)
	{
		const SubControlVolumeFace<TElem>& scvf = m_geo->scvf(i);

		for(uint ip = 0; ip < scvf.num_ip(); ++ip)
		{
			const SD_Values<TElem>& sdv = scvf.sdv();

			m_Diff_Tensor(D, scvf.global_ip(), time);
			m_Conv_Vel(v, scvf.global_ip(), time);

			for(uint j = 0; j < sdv.num_sh(); ++j)
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
	for(uint i = 0; i < m_geo->num_scv(); ++i)
	{
		const SubControlVolume<TElem>& scv = m_geo->scv(i);

		co = scv.local_corner_id();

		m_Reaction(reac_val, scv.global_corner_pos(), time);

		J(co , co) += reac_val * scv.volume();
	}
}


template<typename TDomain, typename TElem >
inline
void
ConvectionDiffusionEquation<TDomain, TElem>::
assemble_element_JM(local_matrix_type& J, const local_vector_type& u, number time)
{
	int co;
	for(uint i = 0; i < m_geo->num_scv(); ++i)
	{
		const SubControlVolume<TElem>& scv = m_geo->scv(i);

		co = scv.local_corner_id();

		J(co , co) += 1.0 * scv.volume();
	}
}


template<typename TDomain, typename TElem >
inline
void
ConvectionDiffusionEquation<TDomain, TElem>::
assemble_element_A(local_vector_type& d, const local_vector_type& u, number time)
{
	static const int dim = TDomain::dim;

	number flux;
	MathVector<dim> grad_u;
	number shape_u;
	MathMatrix<dim,dim> D;
	MathVector<dim> v;
	MathVector<dim> Dgrad_u;
	for(uint i = 0; i < m_geo->num_scvf(); ++i)
	{
		const SubControlVolumeFace<TElem>& scvf = m_geo->scvf(i);

		for(uint ip = 0; ip < scvf.num_ip(); ++ip)
		{
			const SD_Values<TElem>& sdv = scvf.sdv();

			VecSet(grad_u, 0.0);
			shape_u = 0.0;
			for(uint j = 0; j < sdv.num_sh(); ++j)
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
	for(uint i = 0; i < m_geo->num_scv(); ++i)
	{
		const SubControlVolume<TElem>& scv = m_geo->scv(i);

		co = scv.local_corner_id();

		m_Reaction(reac_val, scv.global_corner_pos(), time);

		d[co] += reac_val * u[co] * scv.volume();
	}
}


template<typename TDomain, typename TElem >
inline
void
ConvectionDiffusionEquation<TDomain, TElem>::
assemble_element_M(local_vector_type& d, const local_vector_type& u, number time)
{
	int co;
	for(uint i = 0; i < m_geo->num_scv(); ++i)
	{
		const SubControlVolume<TElem>& scv = m_geo->scv(i);

		co = scv.local_corner_id();

		d[co] += u[co] * scv.volume();
	}
}


template<typename TDomain, typename TElem >
inline
void
ConvectionDiffusionEquation<TDomain, TElem>::
assemble_element_f(local_vector_type& d, number time)
{
	number fvalue = 0.0;
	for(uint i = 0; i < m_geo->num_scv(); ++i)
	{
		const SubControlVolume<TElem>& scv = m_geo->scv(i);

		m_Rhs(fvalue, scv.global_corner_pos(), time);
		d[scv.local_corner_id()] += fvalue * scv.volume();
	}
}


} // namespace ug


#endif
