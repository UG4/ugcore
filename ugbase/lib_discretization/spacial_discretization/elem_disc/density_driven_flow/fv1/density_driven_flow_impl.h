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
DensityDrivenFlowElemDisc<TFVGeom, TDomain, TAlgebra>::
DensityDrivenFlowElemDisc(	TDomain& domain, number upwind_amount,
							Porosity_fct Porosity, Viscosity_fct Viscosity,
							Density_fct Density, D_Density_fct D_Density,
							Mol_Diff_Tensor_fct Mol_Diff,
							Permeability_Tensor_fct Permeability_Tensor,
							Gravity_fct Gravity)
: 	m_pDomain(&domain), m_upwindAmount(upwind_amount),
	m_PorosityFct(Porosity), m_ViscosityFct(Viscosity),
	m_DensityFct(Density), m_DDensityFct(D_Density),
	m_MolDiffTensorFct(Mol_Diff), m_PermeabilityTensorFct(Permeability_Tensor),
	m_GravityFct(Gravity)
{
	register_assemble_functions(Int2Type<dim>());
};


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
							const std::vector<MathVector<dim> >& vLocalGrad)
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
	MatVecMult(DarcyVel, m_Permeability, Vel);
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
								const std::vector<MathVector<dim> >& vGlobalGrad)
{
//	Number of shapes
	const size_t num_sh = vShape.size();

//	Variables
	number Rho; 				// Density
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
		m_DDensityFct(Rho, c);

	//	Derivative of Vel
		for(size_t sh = 0; sh < num_sh; ++sh)
		{
			// D_vel_c[sh] = Rho'(c) * phi_sh * Gravity
			VecScale(D_vel_c[sh], m_Gravity, Rho * vShape[sh]);
		}
	}

// 	Read in user data
	m_ViscosityFct(Viscosity, c);

// 	Compute Darcy velocity
	VecSubtract(Vel, Vel, grad_p);
	MatVecMult(DarcyVel, m_Permeability, Vel);
	VecScale(DarcyVel, DarcyVel, 1./Viscosity);

// 	Compute derivative of Darcy Velocity with respect to _C_ and _P_
// 	Out of the parameters, only the density depends on c
	for(size_t sh = 0; sh < num_sh; ++sh)
	{
		VecScale(D_vel_p[sh], vGlobalGrad[sh], -1);
		MatVecMult(D_DarcyVel_c[sh], m_Permeability, D_vel_c[sh]);
		MatVecMult(D_DarcyVel_p[sh], m_Permeability, D_vel_p[sh]);

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
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;

// 	Remember position attachement
	if(m_pDomain == NULL)
	{
		UG_LOG("ERROR in 'DensityDrivenFlowElemDisc::prepare_element_loop':"
				" Domain not set.");
		return false;
	}
	m_aaPos = m_pDomain->get_position_accessor();

// 	User function (constant in whole domain)
	m_PorosityFct(m_Porosity);
	m_GravityFct(m_Gravity);

	// todo: move this to prepare elem and evaluate at midpoint.
	m_PermeabilityTensorFct(m_Permeability);
	m_MolDiffTensorFct(m_Diffusion);
	MatScale(m_Diffusion, m_Porosity, m_Diffusion);

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
	for(size_t i = 0; i < num_co; ++i)
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
	D3F_PROFILE_FUNC();

	D3F_PROFILE_BEGIN(D3f_Setup);
// 	Get finite volume geometry
	static TFVGeom<TElem, dim>& geo =
				FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();

//	Number of Corners
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
	static const size_t num_co = ref_elem_type::num_corners;

//	Compute density in corners
	for(size_t sh = 0; sh < num_co; sh++)
	{
		m_DensityFct(m_vDensity[sh], u(_C_, sh));
		m_DDensityFct(m_vDDensity[sh], u(_C_, sh));
	}

//	Prepare Consistent Gravity if needed
	if(m_UseConsistentGravity)
	{
	// 	Prepare Density in Corners
		if(!PrepareConsistentGravity<dim>(	&m_vConsGravity[0],
											num_co,
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
												num_co,
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
	number flux, flux_c, flux_p; 			// Flux and its derivative at ip
	MathMatrix<dim,dim> D;					// Diffusion Tensor
	MathVector<dim> grad_p_ip, grad_c_ip;	// Gradients
	MathVector<dim> DarcyVel, D_DarcyVel_c[num_co], D_DarcyVel_p[num_co];
	number c_ip;
	MathVector<dim> Dgrad;
	D3F_PROFILE_END();

	D3F_PROFILE_BEGIN(D3f_SCVF_LOOP);
//	Loop Sub Control Volume Faces (SCVF)
	for(size_t i = 0; i < geo.num_scvf(); ++i)
	{
	// 	Get current SCVF
		const typename TFVGeom<TElem, dim>::SCVF& scvf = geo.scvf(i);

	// 	Loop integration point of SCVF
		for(size_t ip = 0; ip < scvf.num_ip(); ++ip)
		{
			D3F_PROFILE_BEGIN(D3f_ComputeIpSol);
		//	Resets
			VecSet(grad_p_ip, 0.0);
			VecSet(grad_c_ip, 0.0);
			c_ip = 0.0;

		//	Compute Gradients and concentration at ip
			for(size_t j = 0; j < scvf.num_sh(); ++j)
			{
				VecScaleAppend(grad_p_ip, u(_P_,j), scvf.global_grad(j, ip));
				VecScaleAppend(grad_c_ip, u(_C_,j), scvf.global_grad(j, ip));
				c_ip += u(_C_, j) * scvf.shape(j, ip);
			}
			D3F_PROFILE_END();

		//	Compute Darcy Velocity
			D3F_PROFILE_BEGIN(D3f_ip_Darcy_Velocity);
			compute_D_ip_Darcy_velocity(	DarcyVel,
											D_DarcyVel_c, D_DarcyVel_p,
											c_ip, grad_p_ip,
											scvf.JTInv(ip),
											scvf.shape_vector(ip),
											scvf.local_grad_vector(ip),
											scvf.global_grad_vector(ip));
			D3F_PROFILE_END();

		//	Loop Shape Functions
			D3F_PROFILE_BEGIN(D3f_CentralDiff);
			for(size_t j = 0; j < scvf.num_sh(); ++j)
			{
				///////////////////////////////////////////
				// Diffusive term (central discretization)
				///////////////////////////////////////////
				MatVecMult(Dgrad, m_Diffusion, scvf.global_grad(j, ip));
				flux = VecDot(Dgrad, scvf.normal());

				J(_C_, scvf.from(), _C_, j) -= flux;
				J(_C_, scvf.to(),   _C_, j) += flux;
			//  J(_C_, scvf.from(), _P_, j) -= 0.0;
			//  J(_C_, scvf.to(),   _P_, j) += 0.0;

				///////////////////////////////////////////////////
				// Convective term
				// (upwinding_amount == 1.0 -> full upwind;
				//  upwinding_amount == 0.0 -> full central disc)
				///////////////////////////////////////////////////

			// 	Central part convection
				if(m_upwindAmount != 1.0)
				{
					flux_c = (1.- m_upwindAmount) * (scvf.shape(j, ip) * VecDot(DarcyVel, scvf.normal())
													+ c_ip * VecDot(D_DarcyVel_c[j], scvf.normal()));
					flux_p = (1.- m_upwindAmount) * (											0.0
													+ c_ip * VecDot(D_DarcyVel_p[j], scvf.normal()));

					// coupling 'from' with j  (i.e. A[from][j]) and 'to' with j (i.e. A[to][j])
					J(_C_, scvf.from(), _C_, j) += flux_c;
					J(_C_, scvf.to(),   _C_, j) -= flux_c;

					J(_C_, scvf.from(), _P_, j) += flux_p;
					J(_C_, scvf.to(),   _P_, j) -= flux_p;
				}
			}
			PROFILE_END();
		// 	Upwind part Convection
			D3F_PROFILE_BEGIN(D3f_UpwindDiff);
			if(m_upwindAmount != 0.0)
			{
				size_t up;
				flux_c = m_upwindAmount * VecDot(DarcyVel, scvf.normal());
				if(flux_c >= 0.0) up = scvf.from(); else up = scvf.to();

				for(size_t j = 0; j < scvf.num_sh(); ++j)
				{
					if(j == up) flux_c = m_upwindAmount * ( VecDot(DarcyVel, scvf.normal()) );
					else flux_c = 0.0;
					flux_c += m_upwindAmount * u(_C_, up) * VecDot(D_DarcyVel_c[j], scvf.normal());

					J(_C_, scvf.from(), _C_, j) += flux_c;
					J(_C_, scvf.to(),   _C_, j) -= flux_c;

					flux_p =  m_upwindAmount * ( u(_C_, up) * VecDot(D_DarcyVel_p[j], scvf.normal()));
					J(_C_, scvf.from(), _P_, j) += flux_p;
					J(_C_, scvf.to(),   _P_, j) -= flux_p;
				}
			}
			D3F_PROFILE_END();

			///////////////////
			// Flow equation
			///////////////////
			D3F_PROFILE_BEGIN(D3f_FlowEq);
			for(size_t j = 0; j < scvf.num_sh(); ++j)
			{
				flux_c = VecDot(D_DarcyVel_c[j], scvf.normal());
				flux_p = VecDot(D_DarcyVel_p[j], scvf.normal());

				J(_P_, scvf.from(), _C_, j) += flux_c;
				J(_P_, scvf.to(),   _C_, j) -= flux_c;
				J(_P_, scvf.from(), _P_, j) += flux_p;
				J(_P_, scvf.to(),   _P_, j) -= flux_p;
			}
			D3F_PROFILE_END();
		}
	}
	D3F_PROFILE_END();

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
	// get finite volume geometry
	static TFVGeom<TElem, dim>& geo =
				FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();

	// loop Sub Control Volumes (SCV)
	for(size_t i = 0; i < geo.num_scv(); ++i)
	{
		// get current SCV
		const typename TFVGeom<TElem, dim>::SCV& scv = geo.scv(i);

		// get associated node
		const int co = scv.node_id();

		// Add to local matrix
		J(_C_, co, _C_, co) += m_Porosity * scv.volume();
		//J(_C_, _P_, co, co) += 0;
		//J(_P_, _C_, co, co) += 0;
		//J(_P_, _P_, co, co) += 0;
	}

	// we're done
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
	static TFVGeom<TElem, dim>& geo =
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
											num_co,
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
	number flux;
	MathVector<dim> grad_p_ip, grad_c_ip;
	number c_ip;
	MathMatrix<dim,dim> D;
	MathVector<dim> Dgrad_c_ip;
	MathVector<dim> DarcyVel;

// 	Loop Sub Control Volume Faces (SCVF)
	for(size_t i = 0; i < geo.num_scvf(); ++i)
	{
	// 	Get current SCVF
		const typename TFVGeom<TElem, dim>::SCVF& scvf = geo.scvf(i);

	// 	Loop integration points
		for(size_t ip = 0; ip < scvf.num_ip(); ++ip)
		{
		//	Resets
			VecSet(grad_p_ip, 0.0); VecSet(grad_c_ip, 0.0);
			c_ip = 0.0;

		//	Compute current solution and derivatives at ip
			for(size_t j = 0; j < scvf.num_sh(); ++j)
			{
				VecScaleAppend(grad_p_ip, u(_P_,j), scvf.global_grad(j, ip));
				VecScaleAppend(grad_c_ip, u(_C_,j), scvf.global_grad(j, ip));
				c_ip += u(_C_, j) * scvf.shape(j, ip);
			}

		//	Compute Darcy velocity
			compute_ip_Darcy_velocity(	DarcyVel,
										c_ip, grad_p_ip,
										scvf.JTInv(ip),
										scvf.local_grad_vector(ip));

			////////////////////////////////////////////
			// Diffusive term (central discretization)
			////////////////////////////////////////////
			MatVecMult(Dgrad_c_ip, m_Diffusion, grad_c_ip);
			flux = VecDot(Dgrad_c_ip, scvf.normal());

			d(_C_, scvf.from()) -= flux;
			d(_C_,scvf.to()) += flux;

			////////////////////////////////////
			// Convective term
			// (upwinding_amount == 1.0 -> full upwind;
			//  upwinding_amount == 0.0 -> full central disc)
			////////////////////////////////////////////

		// 	Central part convection
			if(m_upwindAmount != 1.0)
			{
				flux = (1.- m_upwindAmount) * c_ip * VecDot(DarcyVel, scvf.normal());

				d(_C_,scvf.from()) += flux;
				d(_C_,scvf.to()) -= flux;

			}

		// 	Upwind part convection
			if(m_upwindAmount != 0.0)
			{
				flux = m_upwindAmount * VecDot(DarcyVel, scvf.normal());
				if(flux >= 0.0) flux *= u(_C_,scvf.from());
				else 			flux *= u(_C_,scvf.to());

				d(_C_,scvf.from()) += flux;
				d(_C_,scvf.to()) -= flux;
			}

			flux = VecDot(DarcyVel, scvf.normal());
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
assemble_M(local_vector_type& d, const local_vector_type& u, number time)
{
// 	Get finite volume geometry
	static TFVGeom<TElem, dim>& geo =
				FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();

// 	Loop Sub Control Volumes (SCV)
	for(size_t i = 0; i < geo.num_scv(); ++i)
	{
	// 	Get current SCV
		const typename TFVGeom<TElem, dim>::SCV& scv = geo.scv(i);

	// 	Get associated node
		const int co = scv.node_id();

	// 	Add to local matrix
		d(_C_,co) += m_Porosity * u(_C_,co) * scv.volume();
		d(_P_,co) += m_Porosity * scv.volume();
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
