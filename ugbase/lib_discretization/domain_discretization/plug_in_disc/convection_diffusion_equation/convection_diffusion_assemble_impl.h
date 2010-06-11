/*
 * convection_diffusion_assemble.cpp
 *
 *  Created on: 26.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__DOMAIN_DISCRETIZATION__PLUG_IN_DISC__CONVECTION_DIFFUSION_EQUATION__CONVECTION_DIFFUSION_ASSEMBLE_IMPL__
#define __H__LIB_DISCRETIZATION__DOMAIN_DISCRETIZATION__PLUG_IN_DISC__CONVECTION_DIFFUSION_EQUATION__CONVECTION_DIFFUSION_ASSEMBLE_IMPL__


namespace ug{


////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
//  Provide a generic implementation for all elements
//  (since this discretization can be implemented in a generic way)
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

template<typename TDomain, typename TAlgebra, typename TElem >
inline
IPlugInReturn
ConvectionDiffusionEquation<TDomain, TAlgebra, TElem>::
prepare_element_loop()
{
	// all this will be performed outside of the loop over the elements.
	// Therefore it is not time critical.

	// create new Geometry
	m_geo = new FVElementGeometry<TElem>();
	assert(m_geo != NULL);

	// remember position attachement
	m_aaPos = m_domain.get_position_accessor();

	typename std::vector<MathVector<TDomain::dim> > pos;
	for(size_t i = 0; i < m_geo->num_scvf(); ++i)
	{
		const SubControlVolumeFace<TElem>& scvf = m_geo->scvf(i);

		for(size_t ip = 0; ip < scvf.num_ip(); ++ip)
		{
			pos.push_back(scvf.local_ip());
		}
	}

	// overwrite old positions (maybe form other element type)
	//m_Velocity.set_positions(pos, true);
	//m_Velocity.set_num_eq(reference_element_traits<TElem>::num_corners);

	return IPlugInReturn_OK;
}

template<typename TDomain, typename TAlgebra, typename TElem >
inline
IPlugInReturn
ConvectionDiffusionEquation<TDomain, TAlgebra, TElem>::
finish_element_loop()
{
	// all this will be performed outside of the loop over the elements.
	// Therefore it is not time critical.

	// delete Geometry
	if(m_geo != NULL)
		delete m_geo;

	return IPlugInReturn_OK;
}

template<typename TDomain, typename TAlgebra, typename TElem >
inline
IPlugInReturn
ConvectionDiffusionEquation<TDomain, TAlgebra, TElem>::
prepare_element(TElem* elem, const local_vector_type& u, const local_index_type& glob_ind)
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

	return IPlugInReturn_OK;
}

template<typename TDomain, typename TAlgebra, typename TElem >
inline
IPlugInReturn
ConvectionDiffusionEquation<TDomain, TAlgebra, TElem>::
assemble_element_JA(local_matrix_type& J, const local_vector_type& u, number time)
{

	static const int dim = TDomain::dim;
	size_t ip_pos = 0;

	number flux;
	//number shape_u;
	MathMatrix<dim,dim> D;
	MathVector<dim> v, lin_Defect;
	MathVector<dim> Dgrad;
	for(size_t i = 0; i < m_geo->num_scvf(); ++i)
	{
		const SubControlVolumeFace<TElem>& scvf = m_geo->scvf(i);

		for(size_t ip = 0; ip < scvf.num_ip(); ++ip)
		{
			const SD_Values<TElem>& sdv = scvf.sdv();

			m_Diff_Tensor(D, scvf.global_ip(), time);
			m_Conv_Vel(v, scvf.global_ip(), time);

			//v = m_Velocity[ip_pos];

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

					// linearization of defect
					/*
					if(m_Velocity.num_sh() > 0)
					{
						shape_u = 0.0;
						for(size_t j = 0; j < sdv.num_sh(); ++j)
						{
							shape_u += u[j] * sdv.shape(j);
						}
						shape_u *= (1.- m_upwind_amount);
						lin_Defect = scvf.normal();
						VecScale(lin_Defect, lin_Defect, shape_u);
						m_Velocity.lin_defect(scvf.from(), ip_pos) = lin_Defect;
						VecScale(m_Velocity.lin_defect(scvf.to(), ip_pos), lin_Defect, -1.);
					}*/
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

				// linearization of defect
				/*if(m_Velocity.num_sh() > 0)
				{
					shape_u = u[up] * m_upwind_amount;
					lin_Defect = scvf.normal();
					VecScale(lin_Defect, lin_Defect, shape_u);
					m_Velocity.lin_defect(scvf.from(), ip_pos) = lin_Defect;
					VecScale(m_Velocity.lin_defect(scvf.to(), ip_pos), lin_Defect, -1.);
				}*/
			}

			// next ip position
			++ip_pos;
		}
	}
	int co;
	number reac_val;
	for(size_t i = 0; i < m_geo->num_scv(); ++i)
	{
		const SubControlVolume<TElem>& scv = m_geo->scv(i);

		co = scv.local_corner_id();

		m_Reaction(reac_val, scv.global_corner_pos(), time);

		J(co , co) += reac_val * scv.volume();
	}

	return IPlugInReturn_OK;
}


template<typename TDomain, typename TAlgebra, typename TElem >
inline
IPlugInReturn
ConvectionDiffusionEquation<TDomain, TAlgebra, TElem>::
assemble_element_JM(local_matrix_type& J, const local_vector_type& u, number time)
{
	int co;
	for(size_t i = 0; i < m_geo->num_scv(); ++i)
	{
		const SubControlVolume<TElem>& scv = m_geo->scv(i);

		co = scv.local_corner_id();

		J(co , co) += 1.0 * scv.volume();
	}

	return IPlugInReturn_OK;
}


template<typename TDomain, typename TAlgebra, typename TElem >
inline
IPlugInReturn
ConvectionDiffusionEquation<TDomain, TAlgebra, TElem>::
assemble_element_A(local_vector_type& d, const local_vector_type& u, number time)
{
	static const int dim = TDomain::dim;
//	size_t ip_pos = 0;

	number flux;
	MathVector<dim> grad_u;
	number shape_u;
	MathMatrix<dim,dim> D;
	MathVector<dim> v;
	MathVector<dim> Dgrad_u;
	for(size_t i = 0; i < m_geo->num_scvf(); ++i)
	{
		const SubControlVolumeFace<TElem>& scvf = m_geo->scvf(i);

		for(size_t ip = 0; ip < scvf.num_ip(); ++ip)
		{
			const SD_Values<TElem>& sdv = scvf.sdv();

			VecSet(grad_u, 0.0);
			shape_u = 0.0;
			for(size_t j = 0; j < sdv.num_sh(); ++j)
			{
				VecScaleAppend(grad_u, u[j], sdv.grad_global(j));
				shape_u += u[j] * sdv.shape(j);
			}

			m_Diff_Tensor(D, scvf.global_ip(), time);
			m_Conv_Vel(v, scvf.global_ip(), time);

			//v = m_Velocity[ip_pos++];

			////////////////////////////////////
			// diffusiv term (central discretization)
			MatVecMult(Dgrad_u, D, grad_u);
			flux = VecDot(Dgrad_u, scvf.normal());

			d[scvf.from()] -= flux;
			d[scvf.to()] += flux;

			//UG_LOG("diff flux: " << flux << "\n");

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
			//	UG_LOG("conv flux: " << flux << "\n");
			}

		}
	}
	int co;
	number reac_val;
	for(size_t i = 0; i < m_geo->num_scv(); ++i)
	{
		const SubControlVolume<TElem>& scv = m_geo->scv(i);

		co = scv.local_corner_id();

		m_Reaction(reac_val, scv.global_corner_pos(), time);

		d[co] += reac_val * u[co] * scv.volume();
	}

	return IPlugInReturn_OK;
}


template<typename TDomain, typename TAlgebra, typename TElem >
inline
IPlugInReturn
ConvectionDiffusionEquation<TDomain, TAlgebra, TElem>::
assemble_element_M(local_vector_type& d, const local_vector_type& u, number time)
{
	int co;
	for(size_t i = 0; i < m_geo->num_scv(); ++i)
	{
		const SubControlVolume<TElem>& scv = m_geo->scv(i);

		co = scv.local_corner_id();

		d[co] += u[co] * scv.volume();
	}

	return IPlugInReturn_OK;
}


template<typename TDomain, typename TAlgebra, typename TElem >
inline
IPlugInReturn
ConvectionDiffusionEquation<TDomain, TAlgebra, TElem>::
assemble_element_f(local_vector_type& d, number time)
{
	number fvalue = 0.0;
	for(size_t i = 0; i < m_geo->num_scv(); ++i)
	{
		const SubControlVolume<TElem>& scv = m_geo->scv(i);

		m_Rhs(fvalue, scv.global_corner_pos(), time);
		d[scv.local_corner_id()] += fvalue * scv.volume();
	}

	return IPlugInReturn_OK;
}


} // namespace ug


#endif /*__H__LIB_DISCRETIZATION__DOMAIN_DISCRETIZATION__PLUG_IN_DISC__CONVECTION_DIFFUSION_EQUATION__CONVECTION_DIFFUSION_ASSEMBLE_IMPL__*/
