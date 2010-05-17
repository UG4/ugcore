/*
 * convection_diffusion_assemble.cpp
 *
 *  Created on: 26.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__DOMAIN_DISCRETIZATION__PLUG_IN_DISC__DENSITY_DRIVEN_FLOW__DENSITY_DRIVEN_FLOW_ASSEMBLE_IMPL__
#define __H__LIB_DISCRETIZATION__DOMAIN_DISCRETIZATION__PLUG_IN_DISC__DENSITY_DRIVEN_FLOW__DENSITY_DRIVEN_FLOW_ASSEMBLE_IMPL__


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
DensityDrivenFlow<TDomain, TAlgebra, TElem>::
prepare_element_loop()
{
	// all this will be performed outside of the loop over the elements.
	// Therefore it is not time critical.

	// create new Geometry
	m_geo = new FVElementGeometry<TElem>();
	assert(m_geo != NULL);

	// remember position attachement
	m_aaPos = m_domain.get_position_accessor();

	m_Darcy_Velocity_export.set_eval_function(&DataExportingClass<MathVector<TDomain::dim>, MathVector<TDomain::dim>,TAlgebra>::export1, this);

	return IPlugInReturn_OK;
}

template<typename TDomain, typename TAlgebra, typename TElem >
inline
IPlugInReturn
DensityDrivenFlow<TDomain, TAlgebra, TElem>::
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
DensityDrivenFlow<TDomain, TAlgebra, TElem>::
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

	// user function
	m_Porosity(m_porosity);

	//prepare Export for computation
	m_Darcy_Velocity_export.set_local_solution(u);

	return IPlugInReturn_OK;
}

#define _C_ 0
#define _P_ 1
#define J(fct1, fct2, i, j) ( J( (reference_element_traits<TElem>::num_corners)*(fct1) + i, (reference_element_traits<TElem>::num_corners)*(fct2) + j) )
#define d(fct, i)    ( d[reference_element_traits<TElem>::num_corners*(fct) + (i)])
#define u(fct, i)    ( u[reference_element_traits<TElem>::num_corners*(fct) + (i)])


template<typename TDomain, typename TAlgebra, typename TElem >
inline
IPlugInReturn
DensityDrivenFlow<TDomain, TAlgebra, TElem>::
assemble_element_JA(local_matrix_type& J, const local_vector_type& u, number time)
{
	const int num_co = reference_element_traits<TElem>::num_corners;
	number flux, flux_c, flux_p;
	MathMatrix<dim,dim> D;
	MathVector<dim> grad_p_ip, grad_c_ip;
	MathVector<dim> Darcy_vel, D_Darcy_vel_c[num_co], D_Darcy_vel_p[num_co];
	number c_ip;
	MathVector<dim> Dgrad;
	for(uint i = 0; i < m_geo->num_scvf(); ++i)
	{
		const SubControlVolumeFace<TElem>& scvf = m_geo->scvf(i);

		for(uint ip = 0; ip < scvf.num_ip(); ++ip)
		{
			const SD_Values<TElem>& sdv = scvf.sdv();

			VecSet(grad_p_ip, 0.0); VecSet(grad_c_ip, 0.0);
			c_ip = 0.0;
			for(uint j = 0; j < sdv.num_sh(); ++j)
			{
				VecScaleAppend(grad_p_ip, u(_P_,j), sdv.grad_global(j));
				VecScaleAppend(grad_c_ip, u(_C_,j), sdv.grad_global(j));
				c_ip += u(_C_, j) * sdv.shape(j);
			}

			compute_D_ip_Darcy_velocity(scvf, Darcy_vel, D_Darcy_vel_c, D_Darcy_vel_p, c_ip, grad_p_ip);
			m_Mol_Diff_Tensor(D);

			for(uint j = 0; j < sdv.num_sh(); ++j)
			{
				////////////////////////////////////
				// diffusiv term (central discretization)
				MatVecMult(Dgrad, D, sdv.grad_global(j));
				flux = VecDot(Dgrad, scvf.normal());

				J(_C_, _C_, scvf.from(), j) -= flux;
				J(_C_, _C_, scvf.to(), j) += flux;

				//J(_C_, _P_, scvf.from(), j) -= 0.0;
				//J(_C_, _P_, scvf.to(), j) += 0.0;;
				////////////////////////////////////
				// convective term
				// (upwinding_amount == 1.0 -> full upwind;
				//  upwinding_amount == 0.0 -> full central disc)

				// central part convection
				if(m_upwind_amount != 1.0)
				{
					flux_c = (1.- m_upwind_amount) * (sdv.shape(j) * VecDot(Darcy_vel, scvf.normal()) + c_ip * VecDot(D_Darcy_vel_c[j], scvf.normal()));
					flux_p = (1.- m_upwind_amount) * (											0.0	 + c_ip * VecDot(D_Darcy_vel_p[j], scvf.normal()));

					// coupling 'from' with j  (i.e. A[from][j]) and 'to' with j (i.e. A[to][j])
					J(_C_, _C_, scvf.from(), j) += flux_c;
					J(_C_, _C_, scvf.to(),j) -= flux_c;

					J(_C_, _P_, scvf.from(), j) += flux_p;
					J(_C_, _P_, scvf.to(),j) -= flux_p;
				}
			}
			// upwind part convection
			// TODO: Something wrong here or in Newton solver, since Newton solver does not converge by quadratic rate if upwind is used
			if(m_upwind_amount != 0.0)
			{
				uint up;
				flux_c = m_upwind_amount * VecDot(Darcy_vel, scvf.normal());
				if(flux_c >= 0.0) up = scvf.from(); else up = scvf.to();

				for(uint j = 0; j < sdv.num_sh(); ++j)
				{
					if(j == up) flux_c = m_upwind_amount * ( sdv.shape(up) * VecDot(Darcy_vel, scvf.normal()) );
					else flux_c = 0.0;
					flux_c += m_upwind_amount * u(_C_, up) * VecDot(D_Darcy_vel_c[j], scvf.normal());

					J(_C_, _C_,scvf.from(), j) += flux_c;
					J(_C_, _C_,scvf.to(), j) -= flux_c;

					flux_p =  m_upwind_amount * ( u(_C_, up) * VecDot(D_Darcy_vel_p[j], scvf.normal()));
					J(_C_, _P_,scvf.from(), j) += flux_p;
					J(_C_, _P_,scvf.to(), j) -= flux_p;
				}
			}


			///////////////////
			///// flow equation
			for(uint j = 0; j < sdv.num_sh(); ++j)
			{
				flux_c = VecDot(D_Darcy_vel_c[j], scvf.normal());
				flux_p = VecDot(D_Darcy_vel_p[j], scvf.normal());

				UG_DLOG(LIB_DISC_D3F, 3, "flux_c = " << flux_c << " (scvf.from = " << scvf.from() << ", scvf.to = " <<  scvf.to() << ", j = " << j <<")\n");
				UG_DLOG(LIB_DISC_D3F, 3, "flux_p = " << flux_p << " (scvf.from = " << scvf.from() << ", scvf.to = " <<  scvf.to() << ", j = " << j <<")\n");

				J(_P_, _C_, scvf.from(), j) += flux_c;
				J(_P_, _C_, scvf.to(), j) -= flux_c;

				J(_P_, _P_, scvf.from(), j) += flux_p;
				J(_P_, _P_, scvf.to(), j) -= flux_p;
			}
		}
	}

	return IPlugInReturn_OK;
}


template<typename TDomain, typename TAlgebra, typename TElem >
inline
IPlugInReturn
DensityDrivenFlow<TDomain, TAlgebra, TElem>::
assemble_element_JM(local_matrix_type& J, const local_vector_type& u, number time)
{
	int co;
	for(uint i = 0; i < m_geo->num_scv(); ++i)
	{
		const SubControlVolume<TElem>& scv = m_geo->scv(i);

		co = scv.local_corner_id();

		J(_C_, _C_, co, co) += m_porosity * scv.volume();
		//J(_C_, _P_, co, co) += 0;
		//J(_P_, _C_, co, co) += 0;
		//J(_P_, _P_, co, co) += 0;
	}

	return IPlugInReturn_OK;
}


template<typename TDomain, typename TAlgebra, typename TElem >
inline
IPlugInReturn
DensityDrivenFlow<TDomain, TAlgebra, TElem>::
assemble_element_A(local_vector_type& d, const local_vector_type& u, number time)
{
	number flux;
	MathVector<dim> grad_p_ip, grad_c_ip;
	number c_ip;
	MathMatrix<dim,dim> D;
	MathVector<dim> Dgrad_c_ip;
	MathVector<dim> Darcy_vel;
	for(uint i = 0; i < m_geo->num_scvf(); ++i)
	{
		const SubControlVolumeFace<TElem>& scvf = m_geo->scvf(i);

		for(uint ip = 0; ip < scvf.num_ip(); ++ip)
		{
			const SD_Values<TElem>& sdv = scvf.sdv();

			VecSet(grad_p_ip, 0.0); VecSet(grad_c_ip, 0.0);
			c_ip = 0.0;
			for(uint j = 0; j < sdv.num_sh(); ++j)
			{
				VecScaleAppend(grad_p_ip, u(_P_,j), sdv.grad_global(j));
				VecScaleAppend(grad_c_ip, u(_C_,j), sdv.grad_global(j));
				c_ip += u(_C_, j) * sdv.shape(j);
			}

			compute_ip_Darcy_velocity(Darcy_vel, c_ip, grad_p_ip);

			m_Mol_Diff_Tensor(D);


			////////////////////////////////////
			// diffusiv term (central discretization)
			MatVecMult(Dgrad_c_ip, D, grad_c_ip);
			flux = VecDot(Dgrad_c_ip, scvf.normal());

			d(_C_, scvf.from()) -= flux;
			d(_C_,scvf.to()) += flux;

			////////////////////////////////////
			// convective term
			// (upwinding_amount == 1.0 -> full upwind;
			//  upwinding_amount == 0.0 -> full central disc)

			// central part convection
			if(m_upwind_amount != 1.0)
			{
				flux = (1.- m_upwind_amount) * c_ip * VecDot(Darcy_vel, scvf.normal());

				d(_C_,scvf.from()) += flux;
				d(_C_,scvf.to()) -= flux;

			}

			// upwind part convection
			if(m_upwind_amount != 0.0)
			{
				flux = m_upwind_amount * VecDot(Darcy_vel, scvf.normal());
				if(flux >= 0.0) flux *= u(_C_,scvf.from()); else flux *= u(_C_,scvf.to());
				d(_C_,scvf.from()) += flux;
				d(_C_,scvf.to()) -= flux;
			}

			flux = VecDot(Darcy_vel, scvf.normal());
			d(_P_,scvf.from()) += flux;
			d(_P_,scvf.to()) -= flux;
		}
	}

	return IPlugInReturn_OK;
}


template<typename TDomain, typename TAlgebra, typename TElem >
inline
IPlugInReturn
DensityDrivenFlow<TDomain, TAlgebra, TElem>::
assemble_element_M(local_vector_type& d, const local_vector_type& u, number time)
{
	int co;
	for(uint i = 0; i < m_geo->num_scv(); ++i)
	{
		const SubControlVolume<TElem>& scv = m_geo->scv(i);

		co = scv.local_corner_id();

		d(_C_,co) += m_porosity * u(_C_,co) * scv.volume();
		d(_P_,co) += m_porosity * scv.volume();
	}
	return IPlugInReturn_OK;
}


template<typename TDomain, typename TAlgebra, typename TElem >
inline
IPlugInReturn
DensityDrivenFlow<TDomain, TAlgebra, TElem>::
assemble_element_f(local_vector_type& d, number time)
{
	return IPlugInReturn_OK;
}


} // namespace ug


#endif /*__H__LIB_DISCRETIZATION__DOMAIN_DISCRETIZATION__PLUG_IN_DISC__DENSITY_DRIVEN_FLOW__DENSITY_DRIVEN_FLOW_ASSEMBLE_IMPL__*/
