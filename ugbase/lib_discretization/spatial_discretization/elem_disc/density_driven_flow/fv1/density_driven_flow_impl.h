/*
 * density_driven_flow_impl.h
 *
 *  Created on: 26.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__DENSITY_DRIVEN_FLOW__FV1__DENSITY_DRIVEN_FLOW_IMPL__
#define __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__DENSITY_DRIVEN_FLOW__FV1__DENSITY_DRIVEN_FLOW_IMPL__

//#define PROFILE_D3F
#ifdef PROFILE_D3F
	#define D3F_PROFILE_FUNC()		PROFILE_FUNC()
	#define D3F_PROFILE_BEGIN(name)	PROFILE_BEGIN(name)
	#define D3F_PROFILE_END()		PROFILE_END()
#else
	#define D3F_PROFILE_FUNC()
	#define D3F_PROFILE_BEGIN(name)
	#define D3F_PROFILE_END()
#endif

namespace ug{


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  Provide a generic implementation for all elements
//  (since this discretization can be implemented in a generic way)
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain, typename TAlgebra>
void
DensityDrivenFlowElemDisc<TDomain, TAlgebra>::
compute_darcy_velocity_ip_std(MathVector<dim>& DarcyVel,
                              const MathMatrix<dim, dim>& Permeability,
                              number Viscosity,
                              const MathVector<dim>& DensityTimesGravity,
                              const MathVector<dim>& PressureGrad,
                              bool compDeriv,
                              MathVector<dim>* DarcyVel_c,
                              MathVector<dim>* DarcyVel_p,
                              const MathVector<dim>* DensityTimesGravity_c,
                              const number* Viscosity_c,
                              const MathVector<dim>* PressureGrad_p,
                              size_t numSh)
{
//	Variables
	MathVector<dim> Vel;

//	compute inverse viscocity
	const number InvVisco = 1./Viscosity;

// 	compute rho*g - \nabla p
	VecSubtract(Vel, DensityTimesGravity, PressureGrad);

//	compute Darcy velocity q := K / mu * (rho*g - \nabla p)
	MatVecMult(DarcyVel, Permeability, Vel);
	VecScale(DarcyVel, DarcyVel, InvVisco);

//	if no derivative needed, done for this point
	if(!compDeriv) return;

//	check correct data
	UG_ASSERT(DarcyVel_c != NULL, "zero pointer");
	UG_ASSERT(DarcyVel_p != NULL, "zero pointer");

//	Derivative of Darcy Vel w.r.t. pressure
	for(size_t sh = 0; sh < numSh; ++sh)
	{
	//	scale by permeability
		MatVecMult(DarcyVel_p[sh], Permeability, PressureGrad_p[sh]);

	//	scale by viscosity
		VecScale(DarcyVel_p[sh], DarcyVel_p[sh], -InvVisco);
	}

//	Derivative of Vel w.r.t. brine mass fraction
	for(size_t sh = 0; sh < numSh; ++sh)
	{
	//	scale by permeability
		MatVecMult(DarcyVel_c[sh], Permeability, DensityTimesGravity_c[sh]);

	//	scale by viscosity
		VecScale(DarcyVel_c[sh], DarcyVel_c[sh], InvVisco);
	}

//	Viscosity derivative
	if(Viscosity_c != NULL)
		for(size_t sh = 0; sh < numSh; ++sh)
		{
		//  DarcyVel_c[sh] -= mu_c_sh(c) / mu * q
			VecSubtract(DarcyVel_c[sh], DarcyVel, - Viscosity_c[sh] * InvVisco);
		}
}


template<typename TDomain, typename TAlgebra>
template <typename TElem>
bool
DensityDrivenFlowElemDisc<TDomain, TAlgebra>::
compute_darcy_export_std(const local_vector_type& u, bool compDeriv)
{
// 	Get finite volume geometry
	static const FV1Geometry<TElem, dim>& geo =
				FVGeometryProvider::get_geom<FV1Geometry, TElem,dim>();

//	Constants
	static const size_t numSh = FV1Geometry<TElem, dim>::numSCV;
	static const int refDim = FV1Geometry<TElem, dim>::dim;

//	Loop all series
	for(size_t s = 0; s < m_exDarcyVel.num_series(); ++s)
	{
	//	currently only for FV1 elem dofs implemented
		if(m_exDarcyVel.template local_ips<refDim>(s) != geo.scvf_local_ips())
		{
			UG_LOG("ERROR in 'compute_darcy_export': Currently export"
					" of Darcy Velocity only implemented for FV1 and"
					"SCVF integration points.\n");
			return false;
		}

	//	Loop ips
		for(size_t ip = 0; ip < m_exDarcyVel.num_ip(s); ++ip)
		{
		//	Darcy Velocity to be filled
			MathVector<dim>& DarcyVel = m_exDarcyVel.value(s, ip);

		//	Derivatives to be filled
			MathVector<dim>* DarcyVel_c = NULL;
			MathVector<dim>* DarcyVel_p = NULL;

			const number* Viscosity_c = NULL;
			const number* Density_c = NULL;
			const MathVector<dim>* PressureGrad_p = NULL;

			MathVector<dim> DensityTimesGravity;
			MathVector<dim> DensityTimesGravity_c[numSh];

		//	get data fields to fill
			if(compDeriv)
			{
				DarcyVel_c = m_exDarcyVel.deriv(s, ip, _C_);
				DarcyVel_p = m_exDarcyVel.deriv(s, ip, _P_);

			//	get derivative of viscosity w.r.t mass fraction
				if(!m_imViscosityScvf.constant_data())
					Viscosity_c = m_imViscosityScvf.deriv(ip, _C_);

			//	compute rho_c_sh * g
				if(!m_imDensityScvf.constant_data())
				{
					Density_c = m_imDensityScvf.deriv(ip, _C_);
					for(size_t sh = 0; sh < numSh; ++sh)
						VecScale(DensityTimesGravity_c[sh], m_Gravity, Density_c[sh]);
				}

			//	get derivative of gradient derivative
				PressureGrad_p = m_imPressureGradScvf.deriv(ip, _P_);
			}

		//	Compute DensityTimesGravity = rho * g
			VecScale(DensityTimesGravity, m_Gravity, m_imDensityScvf[ip]);

		//	compute Darcy velocity (and derivative if required)
			compute_darcy_velocity_ip_std(DarcyVel, m_imPermeabilityScvf[ip],
			                              m_imViscosityScvf[ip], DensityTimesGravity,
			                              m_imPressureGradScvf[ip],
			                              compDeriv,
			                              DarcyVel_c, DarcyVel_p, DensityTimesGravity_c,
			                              Viscosity_c, PressureGrad_p, numSh);
		}
	}

//	we're done
	return true;
}


template<typename TDomain, typename TAlgebra>
template <typename TElem>
bool
DensityDrivenFlowElemDisc<TDomain, TAlgebra>::
compute_darcy_export_cons_grav(const local_vector_type& u, bool compDeriv)
{
// 	Get finite volume geometry
	static const FV1Geometry<TElem, dim>& geo =
				FVGeometryProvider::get_geom<FV1Geometry, TElem,dim>();

//	Constants
	static const size_t numSh = FV1Geometry<TElem, dim>::numSCV;
	static const int refDim = FV1Geometry<TElem, dim>::dim;

/// Consistent gravity and its derivative at corners
	MathVector<dim> vConsGravity[dim*dim];
	MathVector<dim> vvDConsGravity[dim*dim][dim*dim];

// 	Prepare Density in Corners
	if(!PrepareConsistentGravity<dim>(
			&vConsGravity[0], numSh, &m_vCornerCoords[0], m_imDensityScv.values(),m_Gravity))
	{
		UG_LOG("ERROR in assemble_JA: Cannot prepare Consistent Gravity.\n");
		return false;
	}

// 	Prepare DensityDerivative in Corners
	if(compDeriv) {
		number DCoVal[numSh]; memset(DCoVal, 0, sizeof(number)*numSh);
		for(size_t sh = 0; sh < numSh; sh++) {
			DCoVal[sh] = m_imDensityScv.deriv(sh, _C_, sh);
			if(!PrepareConsistentGravity<dim>(
					&vvDConsGravity[sh][0], numSh, &m_vCornerCoords[0],
					&DCoVal[0], m_Gravity))
			{
				UG_LOG("ERROR in assemble_JA: Cannot prepare Consistent Gravity.\n");
				return false;
			}
			DCoVal[sh] = 0.0;
		}
	}

//	Loop all series
	for(size_t s = 0; s < m_exDarcyVel.num_series(); ++s)
	{
	//	currently only for FV1 elem dofs implemented
		if(m_exDarcyVel.template local_ips<refDim>(s) != geo.scvf_local_ips())
		{
			UG_LOG("ERROR in 'compute_darcy_export': Currently export"
					" of Darcy Velocity only implemented for FV1 and"
					"SCVF integration points.\n");
			return false;
		}

	//	only one ip present
		const static size_t j = 0;

	//	Loop Sub Control Volume Faces (SCVF)
		for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
		{
		// 	Get current SCVF
			const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo.scvf(ip);

		//	Darcy Velocity to be filled
			MathVector<dim>& DarcyVel = m_exDarcyVel.value(s, ip);

		//	Derivatives to be filled
			MathVector<dim>* DarcyVel_c = NULL;
			MathVector<dim>* DarcyVel_p = NULL;

			const number* Viscosity_c = NULL;
			const MathVector<dim>* PressureGrad_p = NULL;

			MathVector<dim> DensityTimesGravity;
			MathVector<dim> DensityTimesGravity_c[numSh];

		//	get data fields to fill
			if(compDeriv)
			{
				DarcyVel_c = m_exDarcyVel.deriv(s, ip, _C_);
				DarcyVel_p = m_exDarcyVel.deriv(s, ip, _P_);

			//	get derivative of viscosity w.r.t mass fraction
				if(!m_imViscosityScvf.constant_data())
					Viscosity_c = m_imViscosityScvf.deriv(ip, _C_);

			//	compute rho_c_sh * g
				for(size_t sh = 0; sh < numSh; ++sh)
					if(!ComputeConsistentGravity<dim>(
							DensityTimesGravity_c[sh], numSh, scvf.JTInv(j),
							&(scvf.local_grad_vector(j))[0], &vvDConsGravity[sh][0]))
					{
						UG_LOG("ERROR in compute_ip_Darcy_velocity: Cannot "
								"Compute Consistent Gravity.\n");
						return false;
					}

			//	get derivative of gradient derivative
				PressureGrad_p = m_imPressureGradScvf.deriv(ip, _P_);
			}

		//	Compute DensityTimesGravity = rho * g
			if(!ComputeConsistentGravity<dim>(
					DensityTimesGravity, numSh, scvf.JTInv(j),
					&(scvf.local_grad_vector(j))[0], vConsGravity))
			{
				UG_LOG("ERROR in compute_ip_Darcy_velocity: Cannot "
						"Compute Consistent Gravity.\n");
				return false;
			}

		//	compute Darcy velocity (and derivative if required)
			compute_darcy_velocity_ip_std(DarcyVel, m_imPermeabilityScvf[ip],
										  m_imViscosityScvf[ip], DensityTimesGravity,
										  m_imPressureGradScvf[ip],
										  compDeriv,
										  DarcyVel_c, DarcyVel_p, DensityTimesGravity_c,
										  Viscosity_c, PressureGrad_p, numSh);
		}
	}

//	we're done
	return true;
}

template<typename TDomain, typename TAlgebra>
template<typename TElem >
bool
DensityDrivenFlowElemDisc<TDomain, TAlgebra>::
compute_brine_export(const local_vector_type& u, bool compDeriv)
{
// 	Get finite volume geometry
	static const FV1Geometry<TElem, dim>& geo =
				FVGeometryProvider::get_geom<FV1Geometry, TElem,dim>();

//	Constants
	static const size_t numSh = FV1Geometry<TElem, dim>::numSCV;
	static const int refDim = FV1Geometry<TElem, dim>::dim;

	for(size_t s = 0; s < m_exBrine.num_series(); ++s)
	{
	//	FV1 SCVF ip
		if(m_exBrine.template local_ips<refDim>(s)
				== geo.scvf_local_ips())
		{
		//	Loop Sub Control Volume Faces (SCVF)
			for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
			{
			// 	Get current SCVF
				const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo.scvf(ip);

			//	number of ips is fixed
				const static size_t j = 0;

			//	Compute Gradients and concentration at ip
				number& cIP = m_exBrine.value(s, ip);
				cIP = 0.0;
				for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
					cIP += u(_C_, sh) * scvf.shape(sh, j);

				if(compDeriv)
				{
					number* cIP_c = m_exBrine.deriv(s, ip, _C_);
					number* cIP_p = m_exBrine.deriv(s, ip, _P_);

					for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
					{
						cIP_c[sh] = scvf.shape(sh, j);
						cIP_p[sh] = 0.0;
					}
				}
			}
			continue;
		}

	//	FV1 SCV ip
		if(m_exBrine.template local_ips<refDim>(s)
				== geo.scv_local_ips())
		{
		//	Loop Corners
			for(size_t ip = 0; ip < numSh; ip++)
			{
			//	Compute Gradients and concentration at ip
				m_exBrine.value(s, ip) = u(_C_, ip);

				if(compDeriv)
				{
					number* cIP_c = m_exBrine.deriv(s, ip, _C_);
					number* cIP_p = m_exBrine.deriv(s, ip, _P_);

					for(size_t sh = 0; sh < numSh; ++sh)
					{
						cIP_c[sh] = (ip==sh) ? 1.0 : 0.0;
						cIP_p[sh] = 0.0;
					}
				}
			}
			continue;
		}

		// others not implemented
		UG_LOG("Evaluation not implemented.");
		return false;
	}

//	we're done
	return true;
}

template<typename TDomain, typename TAlgebra>
template<typename TElem >
bool
DensityDrivenFlowElemDisc<TDomain, TAlgebra>::
compute_brine_grad_export(const local_vector_type& u, bool compDeriv)
{
// 	Get finite volume geometry
	static const FV1Geometry<TElem, dim>& geo =
				FVGeometryProvider::get_geom<FV1Geometry, TElem,dim>();

//	Constants
	static const int refDim = FV1Geometry<TElem, dim>::dim;

	for(size_t s = 0; s < m_exBrineGrad.num_series(); ++s)
	{
	//	FV1 SCVF ip
		if(m_exBrineGrad.template local_ips<refDim>(s)
				== geo.scvf_local_ips())
		{
		//	Loop Sub Control Volume Faces (SCVF)
			for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
			{
			// 	Get current SCVF
				const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo.scvf(ip);

			//	number of ips is fixed
				const static size_t j = 0;

			//	Compute Gradients and concentration at ip
				MathVector<dim>& cIP = m_exBrineGrad.value(s, ip);

				VecSet(cIP, 0.0);
				for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
					VecScaleAppend(cIP, u(_C_, sh), scvf.global_grad(sh, j));

				if(compDeriv)
				{
					MathVector<dim>* cIP_c = m_exBrineGrad.deriv(s, ip, _C_);
					MathVector<dim>* cIP_p = m_exBrineGrad.deriv(s, ip, _P_);

					for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
					{
						cIP_c[sh] = scvf.global_grad(sh, j);
						VecSet(cIP_p[sh], 0.0);
					}
				}
			}
			continue;
		}

		// others not implemented
		UG_LOG("Evaluation not implemented.");
		return false;
	}

//	we're done
	return true;
}

template<typename TDomain, typename TAlgebra>
template<typename TElem >
bool
DensityDrivenFlowElemDisc<TDomain, TAlgebra>::
compute_pressure_grad_export(const local_vector_type& u, bool compDeriv)
{
// 	Get finite volume geometry
	static const FV1Geometry<TElem, dim>& geo =
				FVGeometryProvider::get_geom<FV1Geometry, TElem,dim>();

//	Constants
	static const int refDim = FV1Geometry<TElem, dim>::dim;

	for(size_t s = 0; s < m_exPressureGrad.num_series(); ++s)
	{
	//	FV1 SCVF ip
		if(m_exPressureGrad.template local_ips<refDim>(s)
				== geo.scvf_local_ips())
		{
		//	number of ips is fixed
			const static size_t j = 0;

			for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
			{
			// 	Get current SCVF
				const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo.scvf(ip);

			//	Compute Gradients and concentration at ip
				MathVector<dim>& pressGrad = m_exPressureGrad.value(s, ip);

				VecSet(pressGrad, 0.0);
				for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
					VecScaleAppend(pressGrad, u(_P_, sh), scvf.global_grad(sh, j));

				if(compDeriv)
				{
					MathVector<dim>* pressGrad_c = m_exPressureGrad.deriv(s, ip, _C_);
					MathVector<dim>* pressGrad_p = m_exPressureGrad.deriv(s, ip, _P_);

					for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
					{
						pressGrad_p[sh] = scvf.global_grad(sh, j);
						VecSet(pressGrad_c[sh], 0.0);
					}
				}
			}
			continue;
		}

		// others not implemented
		UG_LOG("Evaluation not implemented.");
		return false;
	}

//	we're done
	return true;
}

template<typename TDomain, typename TAlgebra>
template<typename TElem >
bool
DensityDrivenFlowElemDisc<TDomain, TAlgebra>::
prepare_element_loop()
{
// 	Resizes
	typedef typename reference_element_traits<TElem>::reference_element_type
					ref_elem_type;

//	reference dimension
	static const int refDim = ref_elem_type::dim;

//	check, that upwind has been set
	if(m_pUpwind == NULL)
	{
		UG_LOG("ERROR in 'DensityDrivenFlowElemDisc::prepare_element_loop':"
				" Upwind has not been set.\n");
		return false;
	}

//	init upwind for element type
	if(!m_pUpwind->template set_geometry_type<FV1Geometry<TElem, dim> >())
	{
		UG_LOG("ERROR in 'DensityDrivenFlowElemDisc::prepare_element_loop':"
				" Cannot init upwind for element type.\n");
		return false;
	}

//	check necessary imports
	if(!m_imBrineScvf.data_given())
	{
		UG_LOG("DensityDrivenFlowElemDisc::prepare_element_loop:"
				" Missing Import: Brine mass fraction.\n");
		return false;
	}
	if(!m_imBrineGradScvf.data_given())
	{
		UG_LOG("DensityDrivenFlowElemDisc::prepare_element_loop:"
				" Missing Import: Gradient of Brine mass fraction.\n");
		return false;
	}
	if(!m_imPressureGradScvf.data_given())
	{
		UG_LOG("DensityDrivenFlowElemDisc::prepare_element_loop:"
				" Missing Import: Gradient of Pressure.\n");
		return false;
	}
	if(!m_imPorosityScvf.data_given())
	{
		UG_LOG("DensityDrivenFlowElemDisc::prepare_element_loop:"
				" Missing user function: Porosity.\n");
		return false;
	}
	if(!m_imPermeabilityScvf.data_given())
	{
		UG_LOG("DensityDrivenFlowElemDisc::prepare_element_loop:"
				" Missing user function: Permeability.\n");
		return false;
	}
	if(!m_imMolDiffusionScvf.data_given())
	{
		UG_LOG("DensityDrivenFlowElemDisc::prepare_element_loop:"
				" Missing user function: Molecular Diffusion.\n");
		return false;
	}
	if(!m_imViscosityScvf.data_given())
	{
		UG_LOG("DensityDrivenFlowElemDisc::prepare_element_loop:"
				" Missing user function: Viscosity.\n");
		return false;
	}
	if(!m_imDensityScvf.data_given())
	{
		UG_LOG("DensityDrivenFlowElemDisc::prepare_element_loop:"
				" Missing user function: Density.\n");
		return false;
	}
	if(!m_imDarcyVelScvf.data_given())
	{
		UG_LOG("DensityDrivenFlowElemDisc::prepare_element_loop:"
				" Missing Import: Darcy Velocity.\n");
		return false;
	}
	if(!m_imPorosityScv.data_given())
	{
		UG_LOG("DensityDrivenFlowElemDisc::prepare_element_loop:"
				" Missing user function: Porosity.\n");
		return false;
	}
	if(!m_imDensityScv.data_given())
	{
		UG_LOG("DensityDrivenFlowElemDisc::prepare_element_loop:"
				" Missing user function: Density.\n");
		return false;
	}
	if(m_imConstGravity.constant_data())
	{
	//	cast gravity export
		ConstUserVector<dim>* pGrav
			= dynamic_cast<ConstUserVector<dim>*>(m_imConstGravity.get_data());

	//	check success
		if(pGrav == NULL)
		{
			UG_LOG("ERROR in prepare_element_loop: Cannot cast constant gravity.\n");
			return false;
		}

	//	evaluate constant data
		MathVector<dim> dummy;
		pGrav->operator ()(m_Gravity, dummy, 0.0);
	}
	else
	{
		UG_LOG("ERROR in prepare_element_loop: Gravity must be constant.\n");
		return false;
	}


//	set local positions for user data
	FV1Geometry<TElem, dim>& geo = FVGeometryProvider::get_geom<FV1Geometry, TElem,dim>();

	const MathVector<refDim>* vSCVFip = geo.scvf_local_ips();
	size_t numSCVFip = geo.num_scvf_local_ips();
	m_imBrineScvf.			template set_local_ips<refDim>(vSCVFip, numSCVFip);
	m_imBrineGradScvf.		template set_local_ips<refDim>(vSCVFip, numSCVFip);
	m_imPressureGradScvf.	template set_local_ips<refDim>(vSCVFip, numSCVFip);
	m_imPorosityScvf.		template set_local_ips<refDim>(vSCVFip, numSCVFip);
	m_imPermeabilityScvf.	template set_local_ips<refDim>(vSCVFip, numSCVFip);
	m_imMolDiffusionScvf.	template set_local_ips<refDim>(vSCVFip, numSCVFip);
	m_imViscosityScvf.		template set_local_ips<refDim>(vSCVFip, numSCVFip);
	m_imDensityScvf.		template set_local_ips<refDim>(vSCVFip, numSCVFip);
	m_imDarcyVelScvf.		template set_local_ips<refDim>(vSCVFip, numSCVFip);

	const MathVector<refDim>* vSCVip = geo.scv_local_ips();
	size_t numSCVip = geo.num_scv_local_ips();
	m_imPorosityScv.	template set_local_ips<refDim>(vSCVip, numSCVip);
	m_imDensityScv.		template set_local_ips<refDim>(vSCVip, numSCVip);

//	we're done
	return true;
}


template<typename TDomain, typename TAlgebra>
bool
DensityDrivenFlowElemDisc<TDomain, TAlgebra>::
time_point_changed(number time)
{
//	set new time point at imports
	m_imBrineScvf.set_time(time);
	m_imBrineGradScvf.set_time(time);
	m_imPressureGradScvf.set_time(time);
	m_imPressureGradScvf.set_time(time);
	m_imPorosityScvf.set_time(time);
	m_imPermeabilityScvf.set_time(time);
	m_imMolDiffusionScvf.set_time(time);
	m_imViscosityScvf.set_time(time);
	m_imDensityScvf.set_time(time);
	m_imDarcyVelScvf.set_time(time);
	m_imPorosityScv.set_time(time);
	m_imDensityScv.set_time(time);

//	this disc does not need the old time solutions, thus, return false
	return false;
}

template<typename TDomain, typename TAlgebra>
template<typename TElem >
bool
DensityDrivenFlowElemDisc<TDomain, TAlgebra>::
finish_element_loop()
{
	return true;
}

template<typename TDomain, typename TAlgebra>
template<typename TElem >
bool
DensityDrivenFlowElemDisc<TDomain, TAlgebra>::
prepare_element(TElem* elem, const local_vector_type& u, const local_index_type& glob_ind)
{
//	reference element
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;

//	get corners
	m_vCornerCoords = this->template get_element_corners<TElem>(elem);

// 	Update Geometry for this element
	static FV1Geometry<TElem, dim>& geo =
				FVGeometryProvider::get_geom<FV1Geometry, TElem,dim>();

	if(!geo.update(elem, this->get_subset_handler(), &m_vCornerCoords[0]))
	{
		UG_LOG("FVConvectionDiffusionElemDisc::prepare_element:"
				" Cannot update Finite Volume Geometry.\n"); return false;
	}

//	set global positions for user data
	const MathVector<dim>* vSCVFip = geo.scvf_global_ips();
	size_t numSCVFip = geo.num_scvf_global_ips();
	m_imBrineScvf.			set_global_ips(vSCVFip, numSCVFip);
	m_imBrineGradScvf.		set_global_ips(vSCVFip, numSCVFip);
	m_imPressureGradScvf.	set_global_ips(vSCVFip, numSCVFip);
	m_imPorosityScvf.		set_global_ips(vSCVFip, numSCVFip);
	m_imPermeabilityScvf.	set_global_ips(vSCVFip, numSCVFip);
	m_imMolDiffusionScvf.	set_global_ips(vSCVFip, numSCVFip);
	m_imViscosityScvf.		set_global_ips(vSCVFip, numSCVFip);
	m_imDensityScvf.		set_global_ips(vSCVFip, numSCVFip);
	m_imDarcyVelScvf.		set_global_ips(vSCVFip, numSCVFip);

	const MathVector<dim>* vSCVip = geo.scv_global_ips();
	size_t numSCVip = geo.num_scv_global_ips();
	m_imPorosityScv.		set_global_ips(vSCVip, numSCVip);
	m_imDensityScv.			set_global_ips(vSCVip, numSCVip);

//	we're done
	return true;
}

template<typename TDomain, typename TAlgebra>
template<typename TElem >
bool
DensityDrivenFlowElemDisc<TDomain, TAlgebra>::
assemble_JA(local_matrix_type& J, const local_vector_type& u)
{
// 	Get finite volume geometry
	static const FV1Geometry<TElem, dim>& geo =
				FVGeometryProvider::get_geom<FV1Geometry, TElem,dim>();

//	static numbers, known at compile-time
	static const size_t numSh = FV1Geometry<TElem, dim>::numSCV;
	static const size_t numScvf = FV1Geometry<TElem, dim>::numSCVF;

// 	Some variables
	MathMatrix<dim,dim> D, Diffusion[numScvf];
	MathVector<dim> Dgrad;
	number vDFlux_c[numSh], vDFlux_p[numSh];

//	compute Diffusion part of the Diffusion-Dispersion Tensor
	for(size_t s = 0; s < numScvf; ++s)
		MatScale(Diffusion[s], m_imPorosityScvf[s], m_imMolDiffusionScvf[s]);

//	compute Dispersion part
	// todo: Compute DarcyVelocity for Dispersion
	// todo: Compute Dispersion

//	compute upwind shapes for transport equation
	if(!m_pUpwind->update(&geo, m_imDarcyVelScvf.values(), Diffusion, true))
	{
		UG_LOG("ERROR in 'DensityDrivenFlowElemDisc::assemble_JA': "
				"Cannot compute convection shapes.\n");
		return false;
	}

//	get a const (!!) reference to the upwind
	const IConvectionShapes<dim>& convShape
		= *const_cast<const IConvectionShapes<dim>*>(m_pUpwind);

//	there is only on integration point for first order
	const size_t j = 0;

//	Loop Sub Control Volume Faces (SCVF)
	for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
	{
	// 	Get current SCVF
		const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo.scvf(ip);

		///////////////////////////////////////////
		// Darcy Velocity at ip
		///////////////////////////////////////////

		const MathVector<dim>* vDDarcyVel_c = m_imDarcyVelScvf.deriv(ip, _C_);
		const MathVector<dim>* vDDarcyVel_p = m_imDarcyVelScvf.deriv(ip, _P_);

		////////////////////////
		// Transport Equation
		////////////////////////

	//	Loop Shape Functions
		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
		{
		//	Compute Derivative of Convective Flux
			vDFlux_c[sh] = convShape.conv_shape(ip, sh);
			vDFlux_p[sh] = 0.0;

		//	Derivative w.r.t. Velocity
			for(size_t sh1 = 0; sh1 < scvf.num_sh(); ++sh1)
			{
				vDFlux_c[sh] += u(_C_, sh1) * VecDot(convShape.conv_shape_vel(ip, sh1), vDDarcyVel_c[sh]);
				vDFlux_p[sh] += u(_C_, sh1) * VecDot(convShape.conv_shape_vel(ip, sh1), vDDarcyVel_p[sh]);
			}

			//todo: Derivative of Dispersion

		//	Add Derivative of Diffusive Flux
			MatVecMult(Dgrad, Diffusion[ip], scvf.global_grad(sh, j));
			vDFlux_c[sh] -= VecDot(Dgrad, scvf.normal());

			// todo: Derivative of Dispersion
		}

	//	Handle density in case of full equation
		if(!m_BoussinesqTransport)
		{
		//	Convective Flux
			number flux = 0;
			for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
				flux += convShape.conv_shape(ip, sh) * u(_C_, sh);

		//	Diffusive Flux
			MatVecMult(Dgrad, Diffusion[ip], m_imBrineGradScvf[ip]);
			flux -= VecDot(Dgrad, scvf.normal());

		//	Derivative of product
			for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
			{
				vDFlux_c[sh] = m_imDensityScvf[ip] * vDFlux_c[sh] +
							   m_imDensityScvf.deriv(ip, _C_, sh) * flux;
				vDFlux_p[sh] *= m_imDensityScvf[ip];
			}
		}

	//	Add Flux contribution
		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
		{
			J(_C_, scvf.from(), _C_, sh) += vDFlux_c[sh];
			J(_C_, scvf.to(),   _C_, sh) -= vDFlux_c[sh];
			J(_C_, scvf.from(), _P_, sh) += vDFlux_p[sh];
			J(_C_, scvf.to(),   _P_, sh) -= vDFlux_p[sh];
		}

		///////////////////
		// Flow equation
		///////////////////

		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
		{
			vDFlux_c[sh] = VecDot(vDDarcyVel_c[sh], scvf.normal());
			vDFlux_p[sh] = VecDot(vDDarcyVel_p[sh], scvf.normal());
		}

		if(!m_BoussinesqFlow)
		{
			number flux = VecDot(m_imDarcyVelScvf[ip], scvf.normal());
			for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
			{
				vDFlux_c[sh] = m_imDensityScvf[ip] * vDFlux_c[sh] +
							   m_imDensityScvf.deriv(ip, _C_, sh) * flux;
				vDFlux_p[sh] *= m_imDensityScvf[ip];
			}
		}

	//	Add Flux contribution
		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
		{
			J(_P_, scvf.from(), _C_, sh) += vDFlux_c[sh];
			J(_P_, scvf.to(),   _C_, sh) -= vDFlux_c[sh];
			J(_P_, scvf.from(), _P_, sh) += vDFlux_p[sh];
			J(_P_, scvf.to(),   _P_, sh) -= vDFlux_p[sh];
		}
	}

//	we're done
	return true;
}

template<typename TDomain, typename TAlgebra>
template<typename TElem >
bool
DensityDrivenFlowElemDisc<TDomain, TAlgebra>::
assemble_A(local_vector_type& d, const local_vector_type& u)
{
//	Get finite volume geometry
	static const FV1Geometry<TElem, dim>& geo =
			FVGeometryProvider::get_geom<FV1Geometry, TElem,dim>();

//	number of subcontrol volume faces
	static const size_t numScvf = FV1Geometry<TElem, dim>::numSCVF;

//	Some variables
	MathMatrix<dim,dim> D, Diffusion[numScvf];
	MathVector<dim> Dgrad_c_ip;

	for(size_t s = 0; s < numScvf; ++s)
		MatScale(Diffusion[s], m_imPorosityScvf[s], m_imMolDiffusionScvf[s]);

	// todo: Compute Darcy velocity for Dispersion
	// todo: Compute Dispersion
	// todo: Use DiffDisp = Diffusion + Dispersion

//	compute upwind shapes for transport equation
	if(!m_pUpwind->update(&geo, m_imDarcyVelScvf.values(), Diffusion, false))
	{
		UG_LOG("ERROR in 'DensityDrivenFlowElemDisc::assemble_JA': "
				"Cannot compute convection shapes.\n");
		return false;
	}

//	get a const (!!) reference to the upwind
	const IConvectionShapes<dim>& convShape
		= *const_cast<const IConvectionShapes<dim>*>(m_pUpwind);

// 	Loop Sub Control Volume Faces (SCVF)
	for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
	{
	// 	Get current SCVF
		const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo.scvf(ip);

		////////////////////////
		// Transport Equation
		////////////////////////

	//	Compute Convective Flux
		number flux = 0;
		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
			flux += convShape.conv_shape(ip, sh) * u(_C_, sh);

	//	Compute Diffusive Flux
		MatVecMult(Dgrad_c_ip, Diffusion[ip], m_imBrineGradScvf[ip]);
		const number diffFlux = VecDot(Dgrad_c_ip, scvf.normal());

	//	Sum total flux
		flux -= diffFlux;
		if(!m_BoussinesqTransport) flux *= m_imDensityScvf[ip];

	//	Add contribution to transport equation
		d(_C_,scvf.from()) += flux;
		d(_C_,scvf.to()) -= flux;

		////////////////////
		// Flow Equation
		////////////////////

	//	Compute flux
		flux = VecDot(m_imDarcyVelScvf[ip], scvf.normal());
		if(!m_BoussinesqFlow) flux *= m_imDensityScvf[ip];

	//	Add contribution to flow equation
		d(_P_,scvf.from()) += flux;
		d(_P_,scvf.to()) -= flux;
	}

//	we're done
	return true;
}


template<typename TDomain, typename TAlgebra>
template<typename TElem >
bool
DensityDrivenFlowElemDisc<TDomain, TAlgebra>::
assemble_JM(local_matrix_type& J, const local_vector_type& u)
{
// 	get finite volume geometry
	static const FV1Geometry<TElem, dim>& geo =
				FVGeometryProvider::get_geom<FV1Geometry, TElem,dim>();

// 	loop Sub Control Volumes (SCV)
	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
	// 	get current SCV
		const typename FV1Geometry<TElem, dim>::SCV& scv = geo.scv(ip);

	// 	get associated node
		const int co = scv.node_id();

	// 	Add to local matrix
		if(m_BoussinesqTransport)
			J(_C_, co, _C_, co) += m_imPorosityScv[ip] * scv.volume();
		else
			J(_C_, co, _C_, co) +=
				m_imPorosityScv[ip] * scv.volume() *
				(m_imDensityScv[ip]+m_imDensityScv.deriv(ip, _C_, co)*u(_C_, co));

		if(!m_BoussinesqFlow)
			J(_P_, co, _C_, co) += m_imPorosityScv[ip]
			                       * m_imDensityScv.deriv(ip, _C_, co) * scv.volume();
		//else
			//J(_P_, co, _C_, co) += 0;

		// Remark: Other addings are zero
		//J(_C_, co, _P_, co) += 0;
		//J(_P_, co, _P_, co) += 0;
	}

// 	we're done
	return true;
}


template<typename TDomain, typename TAlgebra>
template<typename TElem >
bool
DensityDrivenFlowElemDisc<TDomain, TAlgebra>::
assemble_M(local_vector_type& d, const local_vector_type& u)
{
// 	Get finite volume geometry
	static const FV1Geometry<TElem, dim>& geo =
				FVGeometryProvider::get_geom<FV1Geometry, TElem,dim>();

// 	Loop Sub Control Volumes (SCV)
	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
	// 	Get current SCV
		const typename FV1Geometry<TElem, dim>::SCV& scv = geo.scv(ip);

	// 	Get associated node
		const int co = scv.node_id();

	// 	Add to local matrix
		if(m_BoussinesqTransport)
			d(_C_,co) += m_imPorosityScv[ip] * u(_C_,co) * scv.volume();
		else
			d(_C_,co) += m_imPorosityScv[ip] * m_imDensityScv[ip]
			                                 * u(_C_,co) * scv.volume();

		if(m_BoussinesqFlow)
			d(_P_,co) += m_imPorosityScv[ip] * scv.volume();
		else
			d(_P_,co) += m_imPorosityScv[ip] * m_imDensityScv[ip] * scv.volume();
	}

// 	we're done
	return true;
}


template<typename TDomain, typename TAlgebra>
template<typename TElem >
bool
DensityDrivenFlowElemDisc<TDomain, TAlgebra>::
assemble_f(local_vector_type& d)
{
//	currently there is no contribution that does not depend on the solution

	return true;
}

} // namespace ug


#endif /*__H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__DENSITY_DRIVEN_FLOW__FV1__DENSITY_DRIVEN_FLOW_IMPL__*/
