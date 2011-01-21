/*
 * convection_diffusion_assemble.cpp
 *
 *  Created on: 26.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__CPL_CONVECTION_DIFFUSION_EQUATION__CPL_CONVECTION_DIFFUSION_ASSEMBLE_IMPL__
#define __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__CPL_CONVECTION_DIFFUSION_EQUATION__CPL_CONVECTION_DIFFUSION_ASSEMBLE_IMPL__


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
								Mass_Scale_fct mass_scale,
								Reaction_fct reac, Rhs_fct rhs)
: 	m_Velocity("Velocity"), m_SolutionGrad("SolGrad", this, 0),
	m_domain(domain), m_upwind_amount(upwind_amount),
	m_Diff_Tensor(diff), m_Conv_Scale(conv_scale),
	m_Mass_Scale(mass_scale),
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
	if(!m_Velocity.template set_positions<ref_elem_type::dim>(pos, true))
		{UG_LOG("CplConvectionDiffusionElemDisc::prepare_element_loop: Cannot set positions for Velocity.\n"); return false;}

	if(!m_Velocity.set_num_eq_fct(1))
		{UG_LOG("CplConvectionDiffusionElemDisc::prepare_element_loop: Cannot set number of fct.\n"); return false;}
	if(!m_Velocity.set_num_eq_dofs(_C_, ref_elem_type::num_corners))
		{UG_LOG("CplConvectionDiffusionElemDisc::prepare_element_loop: Cannot set number of dofs.\n"); return false;}

	if(!m_SolutionGrad.set_num_fct(1))
		{UG_LOG("CplDensityDrivenFlowElemDisc::prepare_element_loop: Cannot set num_fct for Velocity.\n"); return false;}
	if(!m_SolutionGrad.set_num_dofs(_C_, ref_elem_type::num_corners))
		{UG_LOG("CplDensityDrivenFlowElemDisc::prepare_element_loop: Cannot set num_dofs for Velocity.\n"); return false;}

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

	//prepare Export for computation
	m_SolutionGrad.set_local_solution(u);

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

			// reset lin_defect
			if(m_Velocity.num_fct() > 0)
			{
				for(size_t k = 0; k < sdv.num_sh(); ++k)
				{
					m_Velocity.lin_defect(ip_pos, _C_, k) = 0.0;
				}
			}

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
			}

			// linearization of defect
			if(m_Velocity.num_fct() > 0)
			{
				// central part convection
				if(m_upwind_amount != 1.0)
				{
					lin_Defect = scvf.normal();
					lin_Defect *= conv_scale * (1.- m_upwind_amount);
					number shape_u = 0.0; for(size_t k = 0; k < sdv.num_sh(); ++k) shape_u += u(_C_, k) * sdv.shape(k);
					lin_Defect *= shape_u;

					m_Velocity.lin_defect(ip_pos, _C_, scvf.from()) += lin_Defect;
					m_Velocity.lin_defect(ip_pos, _C_, scvf.to()) -= lin_Defect;
				}

				// upwind part convection
				if(m_upwind_amount != 0.0)
				{
					lin_Defect = scvf.normal();
					lin_Defect *= conv_scale * m_upwind_amount;
					flux = m_upwind_amount * VecDot(v, scvf.normal());
					if(flux >= 0.0) lin_Defect *= u(_C_, scvf.from()); else lin_Defect *= u(_C_, scvf.to());

					m_Velocity.lin_defect(ip_pos, _C_, scvf.from()) += lin_Defect;
					m_Velocity.lin_defect(ip_pos, _C_, scvf.to()) -= lin_Defect;
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
	number mass_scale;
	for(size_t i = 0; i < get_fvgeom<TElem>().num_scv(); ++i)
	{
		const SubControlVolume<TElem, TDomain::dim>& scv = get_fvgeom<TElem>().scv(i);

		co = scv.local_corner_id();

		m_Mass_Scale(mass_scale, scv.global_corner_pos(), time);

		d(_C_, co) += (mass_scale * u(_C_, co)) * scv.volume();
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


#endif /*__H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__CPL_CONVECTION_DIFFUSION_EQUATION__CPL_CONVECTION_DIFFUSION_ASSEMBLE_IMPL__*/
