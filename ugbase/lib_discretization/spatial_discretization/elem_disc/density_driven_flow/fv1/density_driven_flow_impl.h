/*
 * density_driven_flow_impl.h
 *
 *  Created on: 26.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__DENSITY_DRIVEN_FLOW__FV1__DENSITY_DRIVEN_FLOW_IMPL__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__DENSITY_DRIVEN_FLOW__FV1__DENSITY_DRIVEN_FLOW_IMPL__

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
inline
bool
DensityDrivenFlowElemDisc<TFVGeom, TDomain, TAlgebra>::
compute_ip_Darcy_velocity(	MathVector<dim>& DarcyVel,
							number c,
							const MathVector<dim>& grad_p,
							const MathMatrix<dim, dim>& JTInv,
							const std::vector<MathVector<dim> >& vLocalGrad,
							size_t scvfIP)
{
//	Number of shapes
	const size_t num_sh = vLocalGrad.size();

//	Variables
	number Rho; 			// Density
	number Viscosity;		// Viscosity
	MathVector<dim> Vel;	// Velocity

// 	Read in user data
	m_ViscosityFct(Viscosity, c);	// Viscosity = Viscosity(c)

//	Compute Gravity part
	if(m_UseConsistentGravity)
	{
		if(!ComputeConsistentGravity<dim>(	Vel,
											num_sh,
											JTInv,
											&vLocalGrad[0],
											&m_vConsGravity[0]))
		{
			UG_LOG("ERROR in compute_ip_Darcy_velocity: Cannot "
					"Compute Consistent Gravity.\n");
			return false;
		}
	}
	else
	{
		m_DensityFct(Rho, c);				// Rho = Density(c);
		VecScale(Vel, m_Gravity, Rho);		// Vel = Rho*Gravity
	}

// 	Compute Darcy velocity
	VecSubtract(Vel, Vel, grad_p);	// Vel := Rho*Gravity - Grad_p
	MatVecMult(DarcyVel, m_scvfPermeability[scvfIP], Vel);
	VecScale(DarcyVel, DarcyVel, 1./Viscosity);

//	we're done
	return true;
};

template<	template <class TElem, int TWorldDim> class TFVGeom,
			typename TDomain,
			typename TAlgebra>
inline
bool
DensityDrivenFlowElemDisc<TFVGeom, TDomain, TAlgebra>::
compute_D_ip_Darcy_velocity(	MathVector<dim>& DarcyVel,
								MathVector<dim> D_DarcyVel_c[],
								MathVector<dim> D_DarcyVel_p[],
								number c,
								const MathVector<dim>& grad_p,
								const MathMatrix<dim, dim>& JTInv,
								const std::vector<number>& vShape,
								const std::vector<MathVector<dim> >& vLocalGrad,
								const std::vector<MathVector<dim> >& vGlobalGrad,
								size_t scvfIP)
{
//	Number of shapes
	const size_t num_sh = vShape.size();

//	Variables
	number Rho; 				// Density
	number DRho; 				// Derivative of Density
	number Viscosity;			// Viscosity
	MathVector<dim> Vel;		// Velocity

//	Variables
	MathVector<dim> D_vel_c[m_MaxNumCorners], D_vel_p[m_MaxNumCorners];

//	Compute Gravity part
	if(m_UseConsistentGravity)
	{
	//	Vel
		if(!ComputeConsistentGravity<dim>(	Vel,
											num_sh,
											JTInv,
											&vLocalGrad[0],
											m_vConsGravity))
		{
			UG_LOG("ERROR in compute_ip_Darcy_velocity: Cannot "
					"Compute Consistent Gravity.\n");
			return false;
		}

	//	Derivative of Vel
		for(size_t sh = 0; sh < num_sh; ++sh)
		{
			if(!ComputeConsistentGravity<dim>(	D_vel_c[sh],
												num_sh,
												JTInv,
												&vLocalGrad[0],
												&m_vvDConsGravity[sh][0]))
			{
				UG_LOG("ERROR in compute_ip_Darcy_velocity: Cannot "
						"Compute Consistent Gravity.\n");
				return false;
			}
		}
	}
	else
	{
	//	Vel
		m_DensityFct(Rho, c);
		VecScale(Vel, m_Gravity, Rho);		// Rho*Gravity

	//	Evaluate Drho(c)
		m_DDensityFct(DRho, c);

	//	Derivative of Vel
		for(size_t sh = 0; sh < num_sh; ++sh)
		{
			// D_vel_c[sh] = Rho'(c) * phi_sh * Gravity
			VecScale(D_vel_c[sh], m_Gravity, DRho * vShape[sh]);
		}
	}

// 	Read in user data
	m_ViscosityFct(Viscosity, c);

// 	Compute Darcy velocity
	VecSubtract(Vel, Vel, grad_p);
	MatVecMult(DarcyVel, m_scvfPermeability[scvfIP], Vel);
	VecScale(DarcyVel, DarcyVel, 1./Viscosity);

// 	Compute derivative of Darcy Velocity with respect to _C_ and _P_
// 	Out of the parameters, only the density depends on c
	for(size_t sh = 0; sh < num_sh; ++sh)
	{
		VecScale(D_vel_p[sh], vGlobalGrad[sh], -1.);
		MatVecMult(D_DarcyVel_c[sh], m_scvfPermeability[scvfIP], D_vel_c[sh]);
		MatVecMult(D_DarcyVel_p[sh], m_scvfPermeability[scvfIP], D_vel_p[sh]);

		VecScale(D_DarcyVel_c[sh],D_DarcyVel_c[sh],1./Viscosity);
		VecScale(D_DarcyVel_p[sh],D_DarcyVel_p[sh],1./Viscosity);
	}

	// D_Viscosity == 0 !!!! since mu is constant

// 	we're done
	return true;
};


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
	m_vCornerCoords.resize(ref_elem_type::num_corners);

// 	Remember position attachement
	if(m_pDomain == NULL)
	{
		UG_LOG("ERROR in 'DensityDrivenFlowElemDisc::prepare_element_loop':"
				" Domain not set.");
		return false;
	}
	m_aaPos = m_pDomain->get_position_accessor();

//	set local positions for user data
	TFVGeom<TElem, dim>& geo = FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();
	m_scvPorosity.template set_local_ips<ref_elem_type::dim>(geo.scv_local_ips(),
	                                                 geo.num_scv_local_ips());
	m_scvfPorosity.template set_local_ips<ref_elem_type::dim>(geo.scvf_local_ips(),
	                                                 geo.num_scvf_local_ips());
	m_scvfPermeability.template set_local_ips<ref_elem_type::dim>(geo.scvf_local_ips(),
	                                                 geo.num_scvf_local_ips());
	m_scvfMolDiffusion.template set_local_ips<ref_elem_type::dim>(geo.scvf_local_ips(),
	                                                 geo.num_scvf_local_ips());

// 	User function (constant in whole domain)
	if(m_constGravity.constant_data())
	{
	//	cast gravity export
		ConstUserVector<dim>* pGrav
			= dynamic_cast<ConstUserVector<dim>*>(m_constGravity.get_data());

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
	// this loop will be performed inside the loop over the elements.
	// Therefore, it is TIME CRITICAL
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
	static const size_t num_co = ref_elem_type::num_corners;

	// load corners of this element
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
	m_scvPorosity.set_global_ips(geo.scv_global_ips(), geo.num_scv_global_ips());
	m_scvfPorosity.set_global_ips(geo.scvf_global_ips(), geo.num_scvf_global_ips());
	m_scvfPermeability.set_global_ips(geo.scvf_global_ips(), geo.num_scvf_global_ips());
	m_scvfMolDiffusion.set_global_ips(geo.scvf_global_ips(), geo.num_scvf_global_ips());

//	Compute density in corners
	for(size_t sh = 0; sh < num_co; sh++)
	{
		m_DensityFct(m_vDensity[sh], u(_C_, sh));
		m_DDensityFct(m_vDDensity[sh], u(_C_, sh));
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
assemble_JA(local_matrix_type& J, const local_vector_type& u, number time)
{
// 	Get finite volume geometry
	static const TFVGeom<TElem, dim>& geo =
				FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();

//	Number of Corners
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
	static const size_t num_co = ref_elem_type::num_corners;

//	Prepare Consistent Gravity if needed
	if(m_UseConsistentGravity)
	{
	// 	Prepare Density in Corners
		if(!PrepareConsistentGravity<dim>(	&m_vConsGravity[0],
		                                  	m_vCornerCoords.size(),
											&m_vCornerCoords[0],
											&m_vDensity[0],
											m_Gravity))
		{
			UG_LOG("ERROR in assemble_JA: Cannot "
					"prepare Consistent Gravity.\n");
			return false;
		}

	// 	Prepare DensityDerivative in Corners
		number DCoVal[num_co];
		memset(DCoVal, 0, sizeof(number)*num_co);

		for(size_t sh = 0; sh < num_co; sh++)
		{
			DCoVal[sh] = m_vDDensity[sh];
			if(!PrepareConsistentGravity<dim>(	&m_vvDConsGravity[sh][0],
			                                  	m_vCornerCoords.size(),
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

// 	Some variables
	MathMatrix<dim,dim> D;					// Diffusion Tensor
	number c_ip;
	MathVector<dim> grad_p_ip, grad_c_ip;	// Gradients
	MathVector<dim> DarcyVel;
	MathVector<dim> vDDarcyVel_c[num_co], vDDarcyVel_p[num_co];
	MathVector<dim> Dgrad;
	number vConvShape[num_co];
	MathVector<dim> vDConvShape_Vel[num_co];
	MathMatrix<dim,dim> vDConvShape_D[num_co];
	bool bNonZeroDerivD;
	number vDFlux_c[num_co];
	number vDFlux_p[num_co];

//	Loop Sub Control Volume Faces (SCVF)
	size_t scvfIP = 0;
	for(size_t i = 0; i < geo.num_scvf(); ++i)
	{
	// 	Get current SCVF
		const typename TFVGeom<TElem, dim>::SCVF& scvf = geo.scvf(i);

	// 	Loop integration point of SCVF
		for(size_t ip = 0; ip < scvf.num_ip(); ++ip, scvfIP++)
		{
			///////////////////////////////////////////
			// IP Values
			///////////////////////////////////////////

		//	Compute Gradients and concentration at ip
			VecSet(grad_p_ip, 0.0); VecSet(grad_c_ip, 0.0); c_ip = 0.0;
			for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
			{
				VecScaleAppend(grad_p_ip, u(_P_,sh), scvf.global_grad(sh, ip));
				VecScaleAppend(grad_c_ip, u(_C_,sh), scvf.global_grad(sh, ip));
				c_ip += u(_C_, sh) * scvf.shape(sh, ip);
			}

		//	Compute density and its derivative at ip
			number rho_ip, DRho_ip;
			if(!m_BoussinesqTransport || !m_BoussinesqFlow)
			{
				m_DensityFct(rho_ip, c_ip);
				m_DDensityFct(DRho_ip, c_ip);
			}

			///////////////////////////////////////////
			// Darcy Velocity at ip
			///////////////////////////////////////////

		//	Compute Darcy Velocity
			compute_D_ip_Darcy_velocity(	DarcyVel,
											vDDarcyVel_c, vDDarcyVel_p,
											c_ip, grad_p_ip,
											scvf.JTInv(ip),
											scvf.shape_vector(ip),
											scvf.local_grad_vector(ip),
											scvf.global_grad_vector(ip),
											scvfIP);

			// todo: Compute DarcyVelocity for Dispersion
			// todo: Compute Dispersion

			MatScale(m_Diffusion, m_scvfPorosity[scvfIP],
			         	 	 	  m_scvfMolDiffusion[scvfIP]);

			////////////////////////
			// Transport Equation
			////////////////////////

		//	Compute Convection Shapes
			switch(m_Upwind)
			{
			case NO_UPWIND:
				ConvectionShapes_NoUpwind(	geo, scvf, DarcyVel, m_Diffusion, vConvShape,
											vDConvShape_Vel, vDConvShape_D, &bNonZeroDerivD);
				break;
			case FULL_UPWIND:
				ConvectionShapes_FullUpwind(geo, scvf, DarcyVel, m_Diffusion, vConvShape,
											vDConvShape_Vel, vDConvShape_D, &bNonZeroDerivD);
				break;
			case PART_UPWIND:
				ConvectionShapes_PartialUpwind(geo, scvf, DarcyVel, m_Diffusion, vConvShape,
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
				MatVecMult(Dgrad, m_Diffusion, scvf.global_grad(sh, ip));
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
				MatVecMult(Dgrad, m_Diffusion,  grad_c_ip);
				flux -= VecDot(Dgrad, scvf.normal());

			//	Derivative of product
				for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
				{
					vDFlux_c[sh] = rho_ip * vDFlux_c[sh] +
									DRho_ip * scvf.shape(sh, ip) * flux;
					vDFlux_p[sh] *= rho_ip;
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
					vDFlux_c[sh] = rho_ip * vDFlux_c[sh] +
									DRho_ip * scvf.shape(sh, ip) * flux;
					vDFlux_p[sh] *= rho_ip;
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
	static const size_t num_co = ref_elem_type::num_corners;

//	Compute density in corners
	for(size_t sh = 0; sh < num_co; sh++)
	{
		m_DensityFct(m_vDensity[sh], u(_C_, sh));
	}

//	Prepare Consistent Gravity if needed
	if(m_UseConsistentGravity)
	{
	// 	Prepare Density in Corners
		if(!PrepareConsistentGravity<dim>(	&m_vConsGravity[0],
		                                  	m_vCornerCoords.size(),
											&m_vCornerCoords[0],
											&m_vDensity[0],
											m_Gravity))
		{
			UG_LOG("ERROR in assemble_A: Cannot "
					"prepare Consistent Gravity.\n");
			return false;
		}
	}

//	Some variables
	MathMatrix<dim,dim> D;
	number c_ip;
	MathVector<dim> grad_p_ip, grad_c_ip;
	MathVector<dim> Dgrad_c_ip;
	MathVector<dim> DarcyVel;
	number vConvShape[num_co];

// 	Loop Sub Control Volume Faces (SCVF)
	size_t scvfIP = 0;
	for(size_t i = 0; i < geo.num_scvf(); ++i)
	{
	// 	Get current SCVF
		const typename TFVGeom<TElem, dim>::SCVF& scvf = geo.scvf(i);

	// 	Loop integration points
		for(size_t ip = 0; ip < scvf.num_ip(); ++ip, ++scvfIP)
		{
			///////////////////////////////
			// Values at Integration Point
			///////////////////////////////

		//	Compute current solution and derivatives at ip
			VecSet(grad_p_ip, 0.0); VecSet(grad_c_ip, 0.0); c_ip = 0.0;
			for(size_t j = 0; j < scvf.num_sh(); ++j)
			{
				VecScaleAppend(grad_p_ip, u(_P_,j), scvf.global_grad(j, ip));
				VecScaleAppend(grad_c_ip, u(_C_,j), scvf.global_grad(j, ip));
				c_ip += u(_C_, j) * scvf.shape(j, ip);
			}

		//	Compute Density at ip iff no Boussinesq-Approximation used
			number rho_ip;
			if(!m_BoussinesqFlow || !m_BoussinesqTransport)
				m_DensityFct(rho_ip, c_ip);

			///////////////////
			// Darcy Velocity
			///////////////////

		//	Compute Darcy velocity
			compute_ip_Darcy_velocity(	DarcyVel,
										c_ip, grad_p_ip,
										scvf.JTInv(ip),
										scvf.local_grad_vector(ip),
										scvfIP);

			// todo: Compute Darcy velocity for Dispersion
			// todo: Compute Dispersion
			// todo: Use DiffDisp = m_Diffusion + Dispersion

			MatScale(m_Diffusion, m_scvfPorosity[scvfIP],
			         	 	 	  m_scvfMolDiffusion[scvfIP]);

			////////////////////////
			// Transport Equation
			////////////////////////

		//	Compute Convection Shapes
			switch(m_Upwind)
			{
			case NO_UPWIND:
				ConvectionShapes_NoUpwind<TFVGeom<TElem, dim> >
					(geo, scvf, DarcyVel, m_Diffusion, vConvShape, NULL, NULL, NULL); break;
			case FULL_UPWIND:
				ConvectionShapes_FullUpwind<TFVGeom<TElem, dim> >
					(geo, scvf, DarcyVel, m_Diffusion, vConvShape, NULL, NULL, NULL); break;
			case PART_UPWIND:
				ConvectionShapes_PartialUpwind<TFVGeom<TElem, dim> >
					(geo, scvf, DarcyVel, m_Diffusion, vConvShape, NULL, NULL, NULL); break;
			default: UG_LOG("Upwind not found.\n"); return false;
			}

		//	Compute Convective Flux
			number flux = 0;
			for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
			{
				flux += vConvShape[sh] * u(_C_, sh);
			}

		//	Compute Diffusive Flux
			MatVecMult(Dgrad_c_ip, m_Diffusion, grad_c_ip);
			const number diffFlux = VecDot(Dgrad_c_ip, scvf.normal());

		//	Sum total flux
			flux -= diffFlux;
			if(!m_BoussinesqTransport) flux *= rho_ip;

		//	Add contribution to transport equation
			d(_C_,scvf.from()) += flux;
			d(_C_,scvf.to()) -= flux;

			////////////////////
			// Flow Equation
			////////////////////

		//	Compute flux
			flux = VecDot(DarcyVel, scvf.normal());
			if(!m_BoussinesqFlow) flux *= rho_ip;

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
			J(_C_, co, _C_, co) += m_scvPorosity[ip] * scv.volume();
		else
			J(_C_, co, _C_, co) += m_scvPorosity[ip] * scv.volume() *
									(m_vDensity[co] + m_vDDensity[co] * u(_C_, co));

		if(!m_BoussinesqFlow)
			J(_P_, co, _C_, co) += m_scvPorosity[ip] * m_vDDensity[co] * scv.volume();
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
			d(_C_,co) += m_scvPorosity[ip] * u(_C_,co) * scv.volume();
		else
			d(_C_,co) += m_scvPorosity[ip] * m_vDensity[co] * u(_C_,co) * scv.volume();

		if(m_BoussinesqFlow)
			d(_P_,co) += m_scvPorosity[ip] * scv.volume();
		else
			d(_P_,co) += m_scvPorosity[ip] * m_vDensity[co] * scv.volume();
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


#endif /*__H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__DENSITY_DRIVEN_FLOW__FV1__DENSITY_DRIVEN_FLOW_IMPL__*/
