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


////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
//  Provide a generic implementation for all elements
//  (since this discretization can be implemented in a generic way)
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////


template<	template <class TElem, int TWorldDim> class TFVGeom,
			typename TDomain,
			typename TAlgebra>
template <typename TElem>
bool
DensityDrivenFlowElemDisc<TFVGeom, TDomain, TAlgebra>::
compute_darcy_export_std(const local_vector_type& u, bool compDeriv)
{
// 	Get finite volume geometry
	static const TFVGeom<TElem, dim>& geo =
				FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();

//	Constants
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
	static const size_t numCo = ref_elem_type::num_corners;
	static const int refDim = ref_elem_type::dim;

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
		//	get data fields to fill
			MathVector<dim>& DarcyVel = m_exDarcyVel.value(s, ip);
			MathVector<dim>* DarcyVel_c = m_exDarcyVel.deriv(s, ip, _C_);
			MathVector<dim>* DarcyVel_p = m_exDarcyVel.deriv(s, ip, _P_);

		//	Variables
			number Viscosity = m_imViscosityScvf[ip];
			MathVector<dim> Vel;		// Velocity

		//	Vel = Rho*Gravity
			VecScale(Vel, m_Gravity, m_imDensityScvf[ip]);

		// 	Compute Darcy velocity
			VecSubtract(Vel, Vel, m_imPressureGradScvf[ip]);
			MatVecMult(DarcyVel, m_imPermeabilityScvf[ip], Vel);
			VecScale(DarcyVel, DarcyVel, 1./Viscosity);

		//	if no derivative needed, done for this point
			if(!compDeriv) continue;

			MathVector<dim> Vel_c[numCo], Vel_p[numCo];

		//	Derivative of Vel w.r.t. brine mass fraction
			for(size_t sh = 0; sh < numCo; ++sh)
			{
			// 	D_vel_c[sh] = Rho'(c) * phi_sh * Gravity
				VecScale(Vel_c[sh], m_Gravity,
									  m_imDensityScvf.deriv(ip, _C_, sh));
			}

		// 	Out of the parameters, only the density depends on c
			for(size_t sh = 0; sh < numCo; ++sh)
			{
				VecScale(Vel_p[sh], m_imPressureGradScvf.deriv(ip, _P_, sh), -1.);
				MatVecMult(DarcyVel_c[sh], m_imPermeabilityScvf[ip], Vel_c[sh]);
				MatVecMult(DarcyVel_p[sh], m_imPermeabilityScvf[ip], Vel_p[sh]);

				VecScale(DarcyVel_c[sh],DarcyVel_c[sh],1./Viscosity);
				VecScale(DarcyVel_p[sh],DarcyVel_p[sh],1./Viscosity);
			}

			// D_Viscosity == 0 !!!! since mu is constant
		}
	}

//	we're done
	return true;
}


template<	template <class TElem, int TWorldDim> class TFVGeom,
			typename TDomain,
			typename TAlgebra>
template <typename TElem>
bool
DensityDrivenFlowElemDisc<TFVGeom, TDomain, TAlgebra>::
compute_darcy_export_cons_grav(const local_vector_type& u, bool compDeriv)
{
// 	Get finite volume geometry
	static const TFVGeom<TElem, dim>& geo =
				FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();

//	Constants
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
	static const size_t numCo = ref_elem_type::num_corners;
	static const int refDim = ref_elem_type::dim;

/// Consistent gravity and its derivative at corners
	MathVector<dim> vConsGravity[dim*dim];
	MathVector<dim> vvDConsGravity[dim*dim][dim*dim];

// 	Prepare Density in Corners
	if(!PrepareConsistentGravity<dim>(	&vConsGravity[0],
										numCo,
										&m_vCornerCoords[0],
										m_imDensityScv.values(),
										m_Gravity))
	{
		UG_LOG("ERROR in assemble_JA: Cannot "
				"prepare Consistent Gravity.\n");
		return false;
	}

	if(compDeriv)
	{
	// 	Prepare DensityDerivative in Corners
		number DCoVal[numCo];
		memset(DCoVal, 0, sizeof(number)*numCo);

		for(size_t sh = 0; sh < numCo; sh++)
		{
			DCoVal[sh] = m_imDensityScv.deriv(sh, _C_, sh);
			if(!PrepareConsistentGravity<dim>(	&vvDConsGravity[sh][0],
												numCo,
												&m_vCornerCoords[0],
												&DCoVal[0],
												m_Gravity))
			{
				UG_LOG("ERROR in assemble_JA: Cannot "
						"prepare Consistent Gravity.\n");
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

	//	Loop Sub Control Volume Faces (SCVF)
		size_t ip = 0;
		for(size_t i = 0; i < geo.num_scvf(); ++i)
		{
		// 	Get current SCVF
			const typename TFVGeom<TElem, dim>::SCVF& scvf = geo.scvf(i);

		// 	Loop integration point of SCVF
			for(size_t j = 0; j < scvf.num_ip(); ++j, ++ip)
			{
			//	get data fields to fill
				MathVector<dim>& DarcyVel = m_exDarcyVel.value(s, ip);
				MathVector<dim>* DarcyVel_c = m_exDarcyVel.deriv(s, ip, _C_);
				MathVector<dim>* DarcyVel_p = m_exDarcyVel.deriv(s, ip, _P_);

			//	Variables
				number Viscosity = m_imViscosityScvf[ip];
				MathVector<dim> Vel;		// Velocity
				MathVector<dim> Vel_c[dim*dim], Vel_p[dim*dim];

			//	Vel
				if(!ComputeConsistentGravity<dim>(	Vel,
													numCo,
													scvf.JTInv(j),
													&(scvf.local_grad_vector(j))[0],
													vConsGravity))
				{
					UG_LOG("ERROR in compute_ip_Darcy_velocity: Cannot "
							"Compute Consistent Gravity.\n");
					return false;
				}

			// 	Compute Darcy velocity
				VecSubtract(Vel, Vel, m_imPressureGradScvf[ip]);
				MatVecMult(DarcyVel, m_imPermeabilityScvf[ip], Vel);
				VecScale(DarcyVel, DarcyVel, 1./Viscosity);

			//	if no derivatives needed, proceed with next point
				if(!compDeriv) continue;

			//	Derivative of Vel
				for(size_t sh = 0; sh < numCo; ++sh)
				{
					if(!ComputeConsistentGravity<dim>(	Vel_c[sh],
														numCo,
														scvf.JTInv(j),
														&(scvf.local_grad_vector(j))[0],
														&vvDConsGravity[sh][0]))
					{
						UG_LOG("ERROR in compute_ip_Darcy_velocity: Cannot "
								"Compute Consistent Gravity.\n");
						return false;
					}
				}

			// 	Out of the parameters, only the density depends on c
				for(size_t sh = 0; sh < numCo; ++sh)
				{
					VecScale(Vel_p[sh], scvf.global_grad(sh, j), -1.);
					MatVecMult(DarcyVel_c[sh], m_imPermeabilityScvf[ip], Vel_c[sh]);
					MatVecMult(DarcyVel_p[sh], m_imPermeabilityScvf[ip], Vel_p[sh]);

					VecScale(DarcyVel_c[sh],DarcyVel_c[sh],1./Viscosity);
					VecScale(DarcyVel_p[sh],DarcyVel_p[sh],1./Viscosity);
				}

				// D_Viscosity == 0 !!!! since mu is constant

			}
		}
	}

//	we're done
	return true;
}

template<	template <class TElem, int TWorldDim> class TFVGeom,
			typename TDomain,
			typename TAlgebra>
template<typename TElem >
bool
DensityDrivenFlowElemDisc<TFVGeom, TDomain, TAlgebra>::
compute_brine_export(const local_vector_type& u, bool compDeriv)
{
// 	Get finite volume geometry
	static const TFVGeom<TElem, dim>& geo =
				FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();

	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
	static const size_t refDim = ref_elem_type::dim;
	static const size_t num_co = ref_elem_type::num_corners;

	for(size_t s = 0; s < m_exBrine.num_series(); ++s)
	{
	//	FV1 SCVF ip
		if(m_exBrine.template local_ips<refDim>(s)
				== geo.scvf_local_ips())
		{
		//	Loop Sub Control Volume Faces (SCVF)
			size_t ip = 0;
			for(size_t i = 0; i < geo.num_scvf(); ++i)
			{
			// 	Get current SCVF
				const typename TFVGeom<TElem, dim>::SCVF& scvf = geo.scvf(i);

			// 	Loop integration point of SCVF
				for(size_t j = 0; j < scvf.num_ip(); ++j, ++ip)
				{
				//	Compute Gradients and concentration at ip
					number& cIP = m_exBrine.value(s, ip);
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
			}
			continue;
		}

	//	FV1 SCV ip
		if(m_exBrine.template local_ips<refDim>(s)
				== geo.scv_local_ips())
		{
		//	Loop Corners
			for(size_t ip = 0; ip < num_co; ip++)
			{
			//	Compute Gradients and concentration at ip
				m_exBrine.value(s, ip) = u(_C_, ip);

				if(compDeriv)
				{
					number* cIP_c = m_exBrine.deriv(s, ip, _C_);
					number* cIP_p = m_exBrine.deriv(s, ip, _P_);

					for(size_t sh = 0; sh < num_co; ++sh)
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

template<	template <class TElem, int TWorldDim> class TFVGeom,
			typename TDomain,
			typename TAlgebra>
template<typename TElem >
bool
DensityDrivenFlowElemDisc<TFVGeom, TDomain, TAlgebra>::
compute_brine_grad_export(const local_vector_type& u, bool compDeriv)
{
// 	Get finite volume geometry
	static const TFVGeom<TElem, dim>& geo =
				FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();

	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
	static const size_t refDim = ref_elem_type::dim;

	for(size_t s = 0; s < m_exBrineGrad.num_series(); ++s)
	{
	//	FV1 SCVF ip
		if(m_exBrineGrad.template local_ips<refDim>(s)
				== geo.scvf_local_ips())
		{
		//	Loop Sub Control Volume Faces (SCVF)
			size_t ip = 0;
			for(size_t i = 0; i < geo.num_scvf(); ++i)
			{
			// 	Get current SCVF
				const typename TFVGeom<TElem, dim>::SCVF& scvf = geo.scvf(i);

			// 	Loop integration point of SCVF
				for(size_t j = 0; j < scvf.num_ip(); ++j, ++ip)
				{
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

template<	template <class TElem, int TWorldDim> class TFVGeom,
			typename TDomain,
			typename TAlgebra>
template<typename TElem >
bool
DensityDrivenFlowElemDisc<TFVGeom, TDomain, TAlgebra>::
compute_pressure_grad_export(const local_vector_type& u, bool compDeriv)
{
// 	Get finite volume geometry
	static const TFVGeom<TElem, dim>& geo =
				FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();

	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
	static const size_t refDim = ref_elem_type::dim;

	for(size_t s = 0; s < m_exPressureGrad.num_series(); ++s)
	{
	//	FV1 SCVF ip
		if(m_exPressureGrad.template local_ips<refDim>(s)
				== geo.scvf_local_ips())
		{
		//	Loop Sub Control Volume Faces (SCVF)
			size_t ip = 0;
			for(size_t i = 0; i < geo.num_scvf(); ++i)
			{
			// 	Get current SCVF
				const typename TFVGeom<TElem, dim>::SCVF& scvf = geo.scvf(i);

			// 	Loop integration point of SCVF
				for(size_t j = 0; j < scvf.num_ip(); ++j, ++ip)
				{
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

template<	template <class TElem, int TWorldDim> class TFVGeom,
			typename TDomain,
			typename TAlgebra>
template<typename TElem >
inline
bool
DensityDrivenFlowElemDisc<TFVGeom, TDomain, TAlgebra>::
prepare_element_loop()
{
// 	Resizes
	typedef typename reference_element_traits<TElem>::reference_element_type
					ref_elem_type;
	static const int refDim = ref_elem_type::dim;
	m_vCornerCoords.resize(ref_elem_type::num_corners);

// 	Remember position attachement
	if(m_pDomain == NULL)
	{
		UG_LOG("ERROR in 'DensityDrivenFlowElemDisc::prepare_element_loop':"
				" Domain not set.");
		return false;
	}
	m_aaPos = m_pDomain->get_position_accessor();

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
	TFVGeom<TElem, dim>& geo = FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();

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

template<	template <class TElem, int TWorldDim> class TFVGeom,
			typename TDomain,
			typename TAlgebra>
template<typename TElem >
inline
bool
DensityDrivenFlowElemDisc<TFVGeom, TDomain, TAlgebra>::
finish_element_loop()
{
	return true;
}

template<	template <class TElem, int TWorldDim> class TFVGeom,
			typename TDomain,
			typename TAlgebra>
template<typename TElem >
inline
bool
DensityDrivenFlowElemDisc<TFVGeom, TDomain, TAlgebra>::
prepare_element(TElem* elem, const local_vector_type& u, const local_index_type& glob_ind)
{
//	reference element
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;

// 	load corners of this element
	for(size_t i = 0; i < m_vCornerCoords.size(); ++i)
	{
		VertexBase* vert = elem->vertex(i);
		m_vCornerCoords[i] = m_aaPos[vert];
	}

// 	Update Geometry for this element
	static TFVGeom<TElem, dim>& geo =
				FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();

	if(!geo.update(elem, m_pDomain->get_subset_handler(), &m_vCornerCoords[0]))
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

template<	template <class TElem, int TWorldDim> class TFVGeom,
			typename TDomain,
			typename TAlgebra>
template<typename TElem >
inline
bool
DensityDrivenFlowElemDisc<TFVGeom, TDomain, TAlgebra>::
assemble_JA(local_matrix_type& J, const local_vector_type& u, number time)
{
// 	Get finite volume geometry
	static const TFVGeom<TElem, dim>& geo =
				FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();

//	Number of Corners
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
	static const size_t numCo = ref_elem_type::num_corners;

// 	Some variables
	MathMatrix<dim,dim> D, Diffusion;					// Diffusion Tensor
	MathVector<dim> Dgrad;
	number vConvShape[numCo];
	MathVector<dim> vDConvShape_Vel[numCo];
	MathMatrix<dim,dim> vDConvShape_D[numCo];
	bool bNonZeroDerivD;
	number vDFlux_c[numCo];
	number vDFlux_p[numCo];

//	Loop Sub Control Volume Faces (SCVF)
	size_t ip = 0;
	for(size_t i = 0; i < geo.num_scvf(); ++i)
	{
	// 	Get current SCVF
		const typename TFVGeom<TElem, dim>::SCVF& scvf = geo.scvf(i);

	// 	Loop integration point of SCVF
		for(size_t j = 0; j < scvf.num_ip(); ++j, ip++)
		{
			///////////////////////////////////////////
			// Darcy Velocity at ip
			///////////////////////////////////////////

			const MathVector<dim>& DarcyVel = m_imDarcyVelScvf[ip];
			const MathVector<dim>* vDDarcyVel_c = m_imDarcyVelScvf.deriv(ip, _C_);
			const MathVector<dim>* vDDarcyVel_p = m_imDarcyVelScvf.deriv(ip, _P_);

			// todo: Compute DarcyVelocity for Dispersion
			// todo: Compute Dispersion

			MatScale(Diffusion, m_imPorosityScvf[ip],
			         	 	 	  m_imMolDiffusionScvf[ip]);

			////////////////////////
			// Transport Equation
			////////////////////////

		//	Compute Convection Shapes
			switch(m_Upwind)
			{
			case NO_UPWIND:
				ConvectionShapes_NoUpwind(	geo, scvf, DarcyVel, Diffusion, vConvShape,
											vDConvShape_Vel, vDConvShape_D, &bNonZeroDerivD);
				break;
			case FULL_UPWIND:
				ConvectionShapes_FullUpwind(geo, scvf, DarcyVel, Diffusion, vConvShape,
											vDConvShape_Vel, vDConvShape_D, &bNonZeroDerivD);
				break;
			case PART_UPWIND:
				ConvectionShapes_PartialUpwind(geo, scvf, DarcyVel, Diffusion, vConvShape,
											vDConvShape_Vel, vDConvShape_D, &bNonZeroDerivD);
				break;
			default: UG_LOG("Upwind not found.\n"); return false;
			}

		//	Loop Shape Functions
			for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
			{
			//	Compute Derivative of Convective Flux
				vDFlux_c[sh] = vConvShape[sh];
				vDFlux_p[sh] = 0.0;

			//	Derivative w.r.t. Velocity
				for(size_t sh1 = 0; sh1 < scvf.num_sh(); ++sh1)
				{
					vDFlux_c[sh] += u(_C_, sh1) * VecDot(vDConvShape_Vel[sh1], vDDarcyVel_c[sh]);
					vDFlux_p[sh] += u(_C_, sh1) * VecDot(vDConvShape_Vel[sh1], vDDarcyVel_p[sh]);
				}

				//todo: Derivative of Dispersion

			//	Add Derivative of Diffusive Flux
				MatVecMult(Dgrad, Diffusion, scvf.global_grad(sh, j));
				vDFlux_c[sh] -= VecDot(Dgrad, scvf.normal());

				// todo: Derivative of Dispersion
			}

		//	Handle density in case of full equation
			if(!m_BoussinesqTransport)
			{
			//	Convective Flux
				number flux = 0;
				for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
					flux += vConvShape[sh] * u(_C_, sh);

			//	Diffusive Flux
				MatVecMult(Dgrad, Diffusion, m_imBrineGradScvf[ip]);
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
				number flux = VecDot(DarcyVel, scvf.normal());
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
	}

//	we're done
	return true;
}

template<	template <class TElem, int TWorldDim> class TFVGeom,
			typename TDomain,
			typename TAlgebra>
template<typename TElem >
inline
bool
DensityDrivenFlowElemDisc<TFVGeom, TDomain, TAlgebra>::
assemble_A(local_vector_type& d, const local_vector_type& u, number time)
{
//	Get finite volume geometry
	static const TFVGeom<TElem, dim>& geo =
			FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();

//	Number of Corners
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
	static const size_t numCo = ref_elem_type::num_corners;

//	Some variables
	MathMatrix<dim,dim> D, Diffusion;
	MathVector<dim> Dgrad_c_ip;
	MathVector<dim> DarcyVel;
	number vConvShape[numCo];

// 	Loop Sub Control Volume Faces (SCVF)
	size_t ip = 0;
	for(size_t i = 0; i < geo.num_scvf(); ++i)
	{
	// 	Get current SCVF
		const typename TFVGeom<TElem, dim>::SCVF& scvf = geo.scvf(i);

	// 	Loop integration points
		for(size_t j = 0; j < scvf.num_ip(); ++j, ++ip)
		{
			///////////////////
			// Darcy Velocity
			///////////////////

			const MathVector<dim>& DarcyVel = m_imDarcyVelScvf[ip];

			// todo: Compute Darcy velocity for Dispersion
			// todo: Compute Dispersion
			// todo: Use DiffDisp = Diffusion + Dispersion

			MatScale(Diffusion, m_imPorosityScvf[ip],
			         	 	 	  m_imMolDiffusionScvf[ip]);

			////////////////////////
			// Transport Equation
			////////////////////////

		//	Compute Convection Shapes
			switch(m_Upwind)
			{
			case NO_UPWIND:
				ConvectionShapes_NoUpwind<TFVGeom<TElem, dim> >
					(geo, scvf, DarcyVel, Diffusion, vConvShape, NULL, NULL, NULL); break;
			case FULL_UPWIND:
				ConvectionShapes_FullUpwind<TFVGeom<TElem, dim> >
					(geo, scvf, DarcyVel, Diffusion, vConvShape, NULL, NULL, NULL); break;
			case PART_UPWIND:
				ConvectionShapes_PartialUpwind<TFVGeom<TElem, dim> >
					(geo, scvf, DarcyVel, Diffusion, vConvShape, NULL, NULL, NULL); break;
			default: UG_LOG("Upwind not found.\n"); return false;
			}

		//	Compute Convective Flux
			number flux = 0;
			for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
			{
				flux += vConvShape[sh] * u(_C_, sh);
			}

		//	Compute Diffusive Flux
			MatVecMult(Dgrad_c_ip, Diffusion, m_imBrineGradScvf[ip]);
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
			flux = VecDot(DarcyVel, scvf.normal());
			if(!m_BoussinesqFlow) flux *= m_imDensityScvf[ip];

		//	Add contribution to flow equation
			d(_P_,scvf.from()) += flux;
			d(_P_,scvf.to()) -= flux;
		}
	}

//	we're done
	return true;
}


template<	template <class TElem, int TWorldDim> class TFVGeom,
			typename TDomain,
			typename TAlgebra>
template<typename TElem >
inline
bool
DensityDrivenFlowElemDisc<TFVGeom, TDomain, TAlgebra>::
assemble_JM(local_matrix_type& J, const local_vector_type& u, number time)
{
// 	get finite volume geometry
	static const TFVGeom<TElem, dim>& geo =
				FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();

// 	loop Sub Control Volumes (SCV)
	size_t ip = 0;
	for(size_t i = 0; i < geo.num_scv(); ++i, ++ip)
	{
	// 	get current SCV
		const typename TFVGeom<TElem, dim>::SCV& scv = geo.scv(i);

	// 	get associated node
		const int co = scv.node_id();

	// 	Add to local matrix
		if(m_BoussinesqTransport)
			J(_C_, co, _C_, co) += m_imPorosityScv[ip] * scv.volume();
		else
			J(_C_, co, _C_, co) += m_imPorosityScv[ip] * scv.volume() *
									(m_imDensityScv[ip] + m_imDensityScv.deriv(ip, _C_, co) * u(_C_, co));

		if(!m_BoussinesqFlow)
			J(_P_, co, _C_, co) += m_imPorosityScv[ip] * m_imDensityScv.deriv(ip, _C_, co) * scv.volume();
		//else
			//J(_P_, co, _C_, co) += 0;

		// Remark: Other addings are zero
		//J(_C_, co, _P_, co) += 0;
		//J(_P_, co, _P_, co) += 0;
	}

// 	we're done
	return true;
}


template<	template <class TElem, int TWorldDim> class TFVGeom,
			typename TDomain,
			typename TAlgebra>
template<typename TElem >
inline
bool
DensityDrivenFlowElemDisc<TFVGeom, TDomain, TAlgebra>::
assemble_M(local_vector_type& d, const local_vector_type& u, number time)
{
// 	Get finite volume geometry
	static const TFVGeom<TElem, dim>& geo =
				FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();

// 	Loop Sub Control Volumes (SCV)
	size_t ip = 0;
	for(size_t i = 0; i < geo.num_scv(); ++i, ++ip)
	{
	// 	Get current SCV
		const typename TFVGeom<TElem, dim>::SCV& scv = geo.scv(i);

	// 	Get associated node
		const int co = scv.node_id();

	// 	Add to local matrix
		if(m_BoussinesqTransport)
			d(_C_,co) += m_imPorosityScv[ip] * u(_C_,co) * scv.volume();
		else
			d(_C_,co) += m_imPorosityScv[ip] * m_imDensityScv[ip] * u(_C_,co) * scv.volume();

		if(m_BoussinesqFlow)
			d(_P_,co) += m_imPorosityScv[ip] * scv.volume();
		else
			d(_P_,co) += m_imPorosityScv[ip] * m_imDensityScv[ip] * scv.volume();
	}

// 	we're done
	return true;
}


template<	template <class TElem, int TWorldDim> class TFVGeom,
			typename TDomain,
			typename TAlgebra>
template<typename TElem >
inline
bool
DensityDrivenFlowElemDisc<TFVGeom, TDomain, TAlgebra>::
assemble_f(local_vector_type& d, number time)
{
	// Here we implement the right hand side and the boundary conditions,
	// that do not depend on the solution

	return true;
}

} // namespace ug


#endif /*__H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__DENSITY_DRIVEN_FLOW__FV1__DENSITY_DRIVEN_FLOW_IMPL__*/
