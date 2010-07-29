/*
 * convection_diffusion_assemble.cpp
 *
 *  Created on: 26.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__CPL_CONVECTION_DIFFUSION_EQUATION__CPL_CONVECTION_DIFFUSION_ASSEMBLE_IMPL__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__CPL_CONVECTION_DIFFUSION_EQUATION__CPL_CONVECTION_DIFFUSION_ASSEMBLE_IMPL__


namespace ug{


////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
//  Provide a generic implementation for all elements
//  (since this discretization can be implemented in a generic way)
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

template<typename TDomain, typename TAlgebra>
CplConvectionDiffusionElemDisc<TDomain, TAlgebra>::
CplConvectionDiffusionElemDisc(TDomain& domain, number upwind_amount,
								Diff_Tensor_fct diff, Conv_Scale_fct conv_scale,
								Mass_Scale_fct mass_scale, Mass_Const_fct mass_const,
								Reaction_fct reac, Rhs_fct rhs)
: 	m_Velocity("Velocity"),
	m_domain(domain), m_upwind_amount(upwind_amount),
	m_Diff_Tensor(diff), m_Conv_Scale(conv_scale),
	m_Mass_Scale(mass_scale), m_Mass_Const(mass_const),
	m_Reaction(reac), m_Rhs(rhs)
{
	register_assemble_functions();
};


template<typename TDomain, typename TAlgebra>
template <typename TElem>
inline
bool
CplConvectionDiffusionElemDisc<TDomain, TAlgebra>::
prepare_element_loop()
{
	// all this will be performed outside of the loop over the elements.
	// Therefore it is not time critical.

	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
	m_corners = new position_type[ref_elem_type::num_corners];

	// remember position attachement
	m_aaPos = m_domain.get_position_accessor();

	typename std::vector<MathVector<ref_elem_type::dim> > pos;
	for(size_t i = 0; i < get_fvgeom<TElem>().num_scvf(); ++i)
	{
		const SubControlVolumeFace<TElem, TDomain::dim>& scvf = get_fvgeom<TElem>().scvf(i);

		for(size_t ip = 0; ip < scvf.num_ip(); ++ip)
		{
			pos.push_back(scvf.local_ip());
		}
	}

	// overwrite old positions (maybe form other element type)
	m_Velocity.template set_positions<ref_elem_type::dim>(pos, true);
	m_Velocity.set_num_eq(ref_elem_type::num_corners);

	return true;
}

template<typename TDomain, typename TAlgebra>
template <typename TElem>
inline
bool
CplConvectionDiffusionElemDisc<TDomain, TAlgebra>::
finish_element_loop()
{
	// all this will be performed outside of the loop over the elements.
	// Therefore it is not time critical.
	delete[] m_corners;

	return true;
}

template<typename TDomain, typename TAlgebra>
template <typename TElem>
inline
bool
CplConvectionDiffusionElemDisc<TDomain, TAlgebra>::
prepare_element(TElem* elem, const local_vector_type& u, const local_index_type& glob_ind)
{
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
	// this loop will be performed inside the loop over the elements.
	// Therefore, it is TIME CRITICAL

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
template <typename TElem>
inline
bool
CplConvectionDiffusionElemDisc<TDomain, TAlgebra>::
assemble_JA(local_matrix_type& J, const local_vector_type& u, number time)
{

	static const int dim = TDomain::dim;
	size_t ip_pos = 0;

	number flux, conv_scale;
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
			v = m_Velocity[ip_pos];
			m_Conv_Scale(conv_scale, scvf.global_ip(), time);
			v *= conv_scale;

			for(size_t j = 0; j < sdv.num_sh(); ++j)
			{
				////////////////////////////////////
				// diffusiv term (central discretization)
				MatVecMult(Dgrad, D, sdv.grad_global(j));
				flux = VecDot(Dgrad, scvf.normal());

				J(_C_, scvf.from(), _C_, j) -= flux;
				J(_C_, scvf.to()  , _C_, j) += flux;

				////////////////////////////////////
				// convective term
				// (upwinding_amount == 1.0 -> full upwind;
				//  upwinding_amount == 0.0 -> full central disc)

				// central part convection
				if(m_upwind_amount != 1.0)
				{
					flux = (1.- m_upwind_amount) * sdv.shape(j) * VecDot(v, scvf.normal());

					// coupling 'from' with j  (i.e. A[from][j]) and 'to' with j (i.e. A[to][j])
					J(_C_, scvf.from(), _C_, j) += flux;
					J(_C_, scvf.to()  , _C_, j) -= flux;

					// linearization of defect

					if(m_Velocity.num_sh() > 0)
					{
						number shape_u = 0.0;
						for(size_t j = 0; j < sdv.num_sh(); ++j)
						{
							shape_u += u(_C_, j) * sdv.shape(j);
						}
						shape_u *= (1.- m_upwind_amount);
						lin_Defect = scvf.normal();
						VecScale(lin_Defect, lin_Defect, shape_u);
						m_Velocity.lin_defect(scvf.from(), ip_pos) = lin_Defect;
						VecScale(m_Velocity.lin_defect(scvf.to(), ip_pos), lin_Defect, -1.);
					}
				}
			}
			// upwind part convection
			if(m_upwind_amount != 0.0)
			{
				int up;
				flux = m_upwind_amount * VecDot(v, scvf.normal());
				if(flux >= 0.0) up = scvf.from(); else up = scvf.to();
				J(_C_, scvf.from(), _C_, up) += flux;
				J(_C_, scvf.to()  , _C_, up) -= flux;

				// linearization of defect
				if(m_Velocity.num_sh() > 0)
				{
					number shape_u = u(_C_, up) * m_upwind_amount;
					lin_Defect = scvf.normal();
					VecScale(lin_Defect, lin_Defect, shape_u);
					m_Velocity.lin_defect(scvf.from(), ip_pos) = lin_Defect;
					VecScale(m_Velocity.lin_defect(scvf.to(), ip_pos), lin_Defect, -1.);
				}
			}

			// next ip position
			++ip_pos;
		}
	}
	int co;
	number reac_val;
	for(size_t i = 0; i < get_fvgeom<TElem>().num_scv(); ++i)
	{
		const SubControlVolume<TElem, TDomain::dim>& scv = get_fvgeom<TElem>().scv(i);

		co = scv.local_corner_id();

		m_Reaction(reac_val, scv.global_corner_pos(), time);

		J(_C_, co , _C_, co) += reac_val * scv.volume();
	}

	return true;
}


template<typename TDomain, typename TAlgebra> template <typename TElem>
inline
bool
CplConvectionDiffusionElemDisc<TDomain, TAlgebra>::
assemble_JM(local_matrix_type& J, const local_vector_type& u, number time)
{
	int co;
	number mass_scale;
	for(size_t i = 0; i < get_fvgeom<TElem>().num_scv(); ++i)
	{
		const SubControlVolume<TElem, TDomain::dim>& scv = get_fvgeom<TElem>().scv(i);

		co = scv.local_corner_id();
		m_Mass_Scale(mass_scale, scv.global_corner_pos(), time);

		J(_C_, co , _C_, co) += mass_scale * scv.volume();
	}

	return true;
}


template<typename TDomain, typename TAlgebra>
template <typename TElem>
inline
bool
CplConvectionDiffusionElemDisc<TDomain, TAlgebra>::
assemble_A(local_vector_type& d, const local_vector_type& u, number time)
{
	static const int dim = TDomain::dim;
	size_t ip_pos = 0;

	number flux, conv_scale;
	MathVector<dim> grad_u;
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
			number shape_u = 0.0;
			for(size_t j = 0; j < sdv.num_sh(); ++j)
			{
				VecScaleAppend(grad_u, u(_C_, j), sdv.grad_global(j));
				shape_u += u(_C_, j) * sdv.shape(j);
			}

			m_Diff_Tensor(D, scvf.global_ip(), time);
			v = m_Velocity[ip_pos++];
			m_Conv_Scale(conv_scale, scvf.global_ip(), time);
			v *= conv_scale;

			////////////////////////////////////
			// diffusiv term (central discretization)
			MatVecMult(Dgrad_u, D, grad_u);
			flux = VecDot(Dgrad_u, scvf.normal());

			d(_C_, scvf.from()) -= flux;
			d(_C_, scvf.to()) += flux;

			////////////////////////////////////
			// convective term
			// (upwinding_amount == 1.0 -> full upwind;
			//  upwinding_amount == 0.0 -> full central disc)

			// central part convection
			if(m_upwind_amount != 1.0)
			{
				flux = (1.- m_upwind_amount) * shape_u * VecDot(v, scvf.normal());

				d(_C_, scvf.from()) += flux;
				d(_C_, scvf.to()) -= flux;
			}

			// upwind part convection
			if(m_upwind_amount != 0.0)
			{
				flux = m_upwind_amount * VecDot(v, scvf.normal());
				if(flux >= 0.0) flux *= u(_C_, scvf.from()); else flux *= u(_C_, scvf.to());
				d(_C_, scvf.from()) += flux;
				d(_C_, scvf.to()) -= flux;
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

		d(_C_, co) += reac_val * u(_C_, co) * scv.volume();
	}

	return true;
}


template<typename TDomain, typename TAlgebra>
template <typename TElem>
inline
bool
CplConvectionDiffusionElemDisc<TDomain, TAlgebra>::
assemble_M(local_vector_type& d, const local_vector_type& u, number time)
{
	int co;
	number mass_scale , mass_const;
	for(size_t i = 0; i < get_fvgeom<TElem>().num_scv(); ++i)
	{
		const SubControlVolume<TElem, TDomain::dim>& scv = get_fvgeom<TElem>().scv(i);

		co = scv.local_corner_id();

		m_Mass_Scale(mass_scale, scv.global_corner_pos(), time);
		m_Mass_Const(mass_const, scv.global_corner_pos(), time);

		d(_C_, co) += (mass_const * mass_scale * u(_C_, co)) * scv.volume();
	}

	return true;
}


template<typename TDomain, typename TAlgebra>
template <typename TElem>
inline
bool
CplConvectionDiffusionElemDisc<TDomain, TAlgebra>::
assemble_f(local_vector_type& d, number time)
{
	number fvalue = 0.0;
	for(size_t i = 0; i < get_fvgeom<TElem>().num_scv(); ++i)
	{
		const SubControlVolume<TElem, TDomain::dim>& scv = get_fvgeom<TElem>().scv(i);

		m_Rhs(fvalue, scv.global_corner_pos(), time);
		d(_C_, scv.local_corner_id()) += fvalue * scv.volume();
	}

	return true;
}


} // namespace ug


#endif /*__H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__CPL_CONVECTION_DIFFUSION_EQUATION__CPL_CONVECTION_DIFFUSION_ASSEMBLE_IMPL__*/
