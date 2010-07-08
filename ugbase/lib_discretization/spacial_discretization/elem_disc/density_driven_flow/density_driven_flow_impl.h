/*
 * density_driven_flow_impl.h
 *
 *  Created on: 26.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__DENSITY_DRIVEN_FLOW__DENSITY_DRIVEN_FLOW_IMPL__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__DENSITY_DRIVEN_FLOW__DENSITY_DRIVEN_FLOW_IMPL__


namespace ug{


////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
//  Provide a generic implementation for all elements
//  (since this discretization can be implemented in a generic way)
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

template<typename TDomain, typename TAlgebra>
DensityDrivenFlowElemDisc<TDomain, TAlgebra>::
DensityDrivenFlowElemDisc(	TDomain& domain, number upwind_amount,
							Pososity_fct Porosity, Viscosity_fct Viscosity, Density_fct Density, D_Density_fct D_Density,
							Mol_Diff_Tensor_fct Mol_Diff, Permeability_Tensor_fct Permeability_Tensor, Gravity_fct Gravity)
: 	m_domain(domain), m_upwind_amount(upwind_amount),
	m_Porosity(Porosity), m_Viscosity(Viscosity), m_Density(Density), m_D_Density(D_Density),
	m_Mol_Diff_Tensor(Mol_Diff), m_Permeability_Tensor(Permeability_Tensor), m_Gravity(Gravity)
{
	IElemDisc<TAlgebra>:: template register_all_assemble_functions<Triangle, 		DensityDrivenFlowElemDisc>(RET_TRIANGLE);
	IElemDisc<TAlgebra>:: template register_all_assemble_functions<Quadrilateral, 	DensityDrivenFlowElemDisc>(RET_QUADRILATERAL);
};


template<typename TDomain, typename TAlgebra>
inline
bool
DensityDrivenFlowElemDisc<TDomain, TAlgebra>::
compute_ip_Darcy_velocity(MathVector<dim>& Darcy_vel, number c_ip, const MathVector<dim>& grad_p_ip)
{
	number s, viscosity_ip;
	MathVector<dim> vel;
	MathMatrix<dim, dim> K;

	m_Density(s, c_ip);
	m_Viscosity(viscosity_ip, c_ip);
	m_Gravity(vel);
	m_Permeability_Tensor(K);
	VecScale(vel, vel, s);
	VecSubtract(vel, vel, grad_p_ip);
	MatVecMult(Darcy_vel, K, vel);
	VecScale(Darcy_vel, Darcy_vel, 1./viscosity_ip);
	return true;
};

template<typename TDomain, typename TAlgebra>
template<typename TElem >
inline
bool
DensityDrivenFlowElemDisc<TDomain, TAlgebra>::
compute_D_ip_Darcy_velocity(	const SubControlVolumeFace<TElem, dim>& scvf,
								MathVector<dim>& Darcy_vel, MathVector<dim> D_Darcy_vel_c[], MathVector<dim> D_Darcy_vel_p[],
								number c_ip, const MathVector<dim>& grad_p_ip)
{
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
	const int num_co = ref_elem_type::num_corners;
	number s, mu_ip;
	MathVector<dim> vel, gravity;
	MathVector<dim> D_vel_c[num_co], D_vel_p[num_co];
	MathMatrix<dim, dim> K;
	const SD_Values<TElem, dim>& sdv = scvf.sdv();

	m_Density(s, c_ip);
	m_Gravity(gravity);
	m_Viscosity(mu_ip, c_ip);
	m_Permeability_Tensor(K);
	VecScale(vel, gravity, s);
	VecSubtract(vel, vel, grad_p_ip);
	MatVecMult(Darcy_vel, K, vel);
	VecScale(Darcy_vel, Darcy_vel, 1./mu_ip);

	m_D_Density(s, c_ip);
	for(int co = 0; co < num_co; ++co)
	{
		VecScale(D_vel_c[co], gravity, s*sdv.shape(co));
		VecScale(D_vel_p[co], sdv.grad_global(co), -1);
		MatVecMult(D_Darcy_vel_c[co], K, D_vel_c[co]);
		MatVecMult(D_Darcy_vel_p[co], K, D_vel_p[co]);

		VecScale(D_Darcy_vel_c[co],D_Darcy_vel_c[co],1./mu_ip);
		VecScale(D_Darcy_vel_p[co],D_Darcy_vel_p[co],1./mu_ip);
	}

	// D_Viscosity == 0 !!!!
	return true;
};


template<typename TDomain, typename TAlgebra>
template<typename TElem >
inline
bool
DensityDrivenFlowElemDisc<TDomain, TAlgebra>::
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
DensityDrivenFlowElemDisc<TDomain, TAlgebra>::
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
DensityDrivenFlowElemDisc<TDomain, TAlgebra>::
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

	// user function
	m_Porosity(m_porosity);

	return true;
}

#define J(fct1, fct2, i, j) ( J( (ref_elem_type::num_corners)*(fct1) + i, (ref_elem_type::num_corners)*(fct2) + j) )
#define d(fct, i)    ( d[ref_elem_type::num_corners*(fct) + (i)])
#define u(fct, i)    ( u[ref_elem_type::num_corners*(fct) + (i)])


template<typename TDomain, typename TAlgebra>
template<typename TElem >
inline
bool
DensityDrivenFlowElemDisc<TDomain, TAlgebra>::
assemble_JA(local_matrix_type& J, const local_vector_type& u, number time)
{
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
	const int num_co = ref_elem_type::num_corners;
	number flux, flux_c, flux_p;
	MathMatrix<dim,dim> D;
	MathVector<dim> grad_p_ip, grad_c_ip;
	MathVector<dim> Darcy_vel, D_Darcy_vel_c[num_co], D_Darcy_vel_p[num_co];
	number c_ip;
	MathVector<dim> Dgrad;
	for(size_t i = 0; i < get_fvgeom<TElem>().num_scvf(); ++i)
	{
		const SubControlVolumeFace<TElem, dim>& scvf = get_fvgeom<TElem>().scvf(i);

		for(size_t ip = 0; ip < scvf.num_ip(); ++ip)
		{
			const SD_Values<TElem, dim>& sdv = scvf.sdv();

			VecSet(grad_p_ip, 0.0); VecSet(grad_c_ip, 0.0);
			c_ip = 0.0;
			for(size_t j = 0; j < sdv.num_sh(); ++j)
			{
				VecScaleAppend(grad_p_ip, u(_P_,j), sdv.grad_global(j));
				VecScaleAppend(grad_c_ip, u(_C_,j), sdv.grad_global(j));
				c_ip += u(_C_, j) * sdv.shape(j);
			}

			compute_D_ip_Darcy_velocity(scvf, Darcy_vel, D_Darcy_vel_c, D_Darcy_vel_p, c_ip, grad_p_ip);
			m_Mol_Diff_Tensor(D);
			MatScale(D, m_porosity, D);

			for(size_t j = 0; j < sdv.num_sh(); ++j)
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
			if(m_upwind_amount != 0.0)
			{
				size_t up;
				flux_c = m_upwind_amount * VecDot(Darcy_vel, scvf.normal());
				if(flux_c >= 0.0) up = scvf.from(); else up = scvf.to();

				for(size_t j = 0; j < sdv.num_sh(); ++j)
				{
					if(j == up) flux_c = m_upwind_amount * ( VecDot(Darcy_vel, scvf.normal()) );
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
			for(size_t j = 0; j < sdv.num_sh(); ++j)
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

	return true;
}

template<typename TDomain, typename TAlgebra>
template<typename TElem >
inline
bool
DensityDrivenFlowElemDisc<TDomain, TAlgebra>::
assemble_JM(local_matrix_type& J, const local_vector_type& u, number time)
{
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
	int co;
	for(size_t i = 0; i < get_fvgeom<TElem>().num_scv(); ++i)
	{
		const SubControlVolume<TElem, dim>& scv = get_fvgeom<TElem>().scv(i);

		co = scv.local_corner_id();

		J(_C_, _C_, co, co) += m_porosity * scv.volume();
		//J(_C_, _P_, co, co) += 0;
		//J(_P_, _C_, co, co) += 0;
		//J(_P_, _P_, co, co) += 0;
	}

	return true;
}


template<typename TDomain, typename TAlgebra>
template<typename TElem >
inline
bool
DensityDrivenFlowElemDisc<TDomain, TAlgebra>::
assemble_A(local_vector_type& d, const local_vector_type& u, number time)
{
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
	number flux;
	MathVector<dim> grad_p_ip, grad_c_ip;
	number c_ip;
	MathMatrix<dim,dim> D;
	MathVector<dim> Dgrad_c_ip;
	MathVector<dim> Darcy_vel;
	for(size_t i = 0; i < get_fvgeom<TElem>().num_scvf(); ++i)
	{
		const SubControlVolumeFace<TElem, dim>& scvf = get_fvgeom<TElem>().scvf(i);

		for(size_t ip = 0; ip < scvf.num_ip(); ++ip)
		{
			const SD_Values<TElem, dim>& sdv = scvf.sdv();

			VecSet(grad_p_ip, 0.0); VecSet(grad_c_ip, 0.0);
			c_ip = 0.0;
			for(size_t j = 0; j < sdv.num_sh(); ++j)
			{
				VecScaleAppend(grad_p_ip, u(_P_,j), sdv.grad_global(j));
				VecScaleAppend(grad_c_ip, u(_C_,j), sdv.grad_global(j));
				c_ip += u(_C_, j) * sdv.shape(j);
			}

			compute_ip_Darcy_velocity(Darcy_vel, c_ip, grad_p_ip);

			m_Mol_Diff_Tensor(D);
			MatScale(D, m_porosity, D);

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

	return true;
}

template<typename TDomain, typename TAlgebra>
template<typename TElem >
inline
bool
DensityDrivenFlowElemDisc<TDomain, TAlgebra>::
assemble_M(local_vector_type& d, const local_vector_type& u, number time)
{
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;

	int co;
	for(size_t i = 0; i < get_fvgeom<TElem>().num_scv(); ++i)
	{
		const SubControlVolume<TElem, dim>& scv = get_fvgeom<TElem>().scv(i);

		co = scv.local_corner_id();

		d(_C_,co) += m_porosity * u(_C_,co) * scv.volume();
		d(_P_,co) += m_porosity * scv.volume();
	}
	return true;
}


template<typename TDomain, typename TAlgebra>
template<typename TElem >
inline
bool
DensityDrivenFlowElemDisc<TDomain, TAlgebra>::
assemble_f(local_vector_type& d, number time)
{
	// Here we implement the right hand side and the boundary conditions, that do not depend on the solution

	return true;
}

#undef J
#undef d
#undef u
} // namespace ug


#endif /*__H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__DENSITY_DRIVEN_FLOW__DENSITY_DRIVEN_FLOW_IMPL__*/
