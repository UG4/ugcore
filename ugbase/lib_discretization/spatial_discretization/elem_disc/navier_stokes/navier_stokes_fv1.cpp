/*
 * navier_stokes_fv1.cpp
 *
 *  Created on: 20.09.2010
 *      Author: andreasvogel
 */

#include "navier_stokes.h"

#include "common/util/provider.h"
#include "lib_discretization/spatial_discretization/disc_util/finite_volume_geometry.h"

namespace ug{


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  Provide a generic implementation for all elements
//  (since this discretization can be implemented in a generic way)
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
bool
FVNavierStokesElemDisc<TDomain>::
prepare_element_loop()
{
// 	Only first order implementation
	if(!(TFVGeom<TElem, dim>::order == 1))
	{
		UG_LOG("ERROR in 'FVNavierStokesElemDisc::prepare_element_loop':"
				" Only first order implementation, but other Finite Volume"
				" Geometry set.\n");
		return false;
	}

//	check, that stabilization has been set
	if(m_pStab == NULL)
	{
		UG_LOG("ERROR in 'FVNavierStokesElemDisc::prepare_element_loop':"
				" Stabilization has not been set.\n");
		return false;
	}

//	init stabilization for element type
	if(!m_pStab->template set_geometry_type<TFVGeom<TElem, dim> >())
	{
		UG_LOG("ERROR in 'FVNavierStokesElemDisc::prepare_element_loop':"
				" Cannot init stabilization for element type.\n");
		return false;
	}

//	check, that convective upwinding has been set
	if(m_pConvStab == NULL  && m_pConvUpwind == NULL)
	{
		UG_LOG("ERROR in 'FVNavierStokesElemDisc::prepare_element_loop':"
				" Upwinding for convective Term in Momentum eq. not set.\n");
		return false;
	}

//	init convection stabilization for element type
	if(m_pConvStab != NULL)
		if(!m_pConvStab->template set_geometry_type<TFVGeom<TElem, dim> >())
		{
			UG_LOG("ERROR in 'FVNavierStokesElemDisc::prepare_element_loop':"
					" Cannot init upwind (PAC) for element type.\n");
			return false;
		}

//	init convection stabilization for element type
	if(m_pConvUpwind != NULL)
		if(!m_pConvUpwind->template set_geometry_type<TFVGeom<TElem, dim> >())
		{
			UG_LOG("ERROR in 'FVNavierStokesElemDisc::prepare_element_loop':"
					" Cannot init upwind for element type.\n");
			return false;
		}

//	check, that kinematic Viscosity has been set
	if(!m_imKinViscosity.data_given())
	{
		UG_LOG("ERROR in 'FVNavierStokesElemDisc::prepare_element_loop':"
				" Kinematic Viscosity has not been set, but is required.\n");
		return false;
	}

//	set local positions for imports
	typedef typename reference_element_traits<TElem>::reference_element_type
																ref_elem_type;
	static const int refDim = ref_elem_type::dim;

	if(!TFVGeom<TElem, dim>::usesHangingNodes)
	{
		TFVGeom<TElem, dim>& geo = Provider<TFVGeom<TElem,dim> >::get();
		m_imKinViscosity.template set_local_ips<refDim>(geo.scvf_local_ips(),
		                                                geo.num_scvf_ips());
		m_imSource.template set_local_ips<refDim>(geo.scv_local_ips(),
		                                          geo.num_scv_ips());
	}

//	we're done
	return true;
}

template<typename TDomain>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
bool
FVNavierStokesElemDisc<TDomain>::
finish_element_loop()
{
//	nothing to do
	return true;
}


template<typename TDomain>
bool
FVNavierStokesElemDisc<TDomain>::
time_point_changed(number time)
{
//	set new time point at imports
	m_imKinViscosity.set_time(time);
	m_imSource.set_time(time);

//	this disc needs the old time solutions, thus, return true
	return true;
}


template<typename TDomain>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
bool
FVNavierStokesElemDisc<TDomain>::
prepare_element(TElem* elem, const local_vector_type& u)
{
//	get corners
	m_vCornerCoords = this->template get_element_corners<TElem>(elem);

// 	Update Geometry for this element
	TFVGeom<TElem, dim>& geo = Provider<TFVGeom<TElem,dim> >::get();
	if(!geo.update(elem, &m_vCornerCoords[0], &(this->get_subset_handler())))
	{
		UG_LOG("FVNavierStokesElemDisc::prepare_element:"
				" Cannot update Finite Volume Geometry.\n"); return false;
	}

//	set local positions for imports
	if(TFVGeom<TElem, dim>::usesHangingNodes)
	{
	//	set local positions for imports
		typedef typename reference_element_traits<TElem>::reference_element_type
																	ref_elem_type;
		static const int refDim = ref_elem_type::dim;

	//	request ip series
		m_imKinViscosity.template set_local_ips<refDim>(geo.scvf_local_ips(),
		                                                geo.num_scvf_ips());
		m_imSource.template set_local_ips<refDim>(geo.scv_local_ips(),
		                                          geo.num_scv_ips());
	}

//	set global positions for imports
	m_imKinViscosity.set_global_ips(geo.scvf_global_ips(), geo.num_scvf_ips());
	m_imSource.set_global_ips(geo.scv_global_ips(), geo.num_scv_ips());

//	we're done
	return true;
}

template<typename TDomain>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
bool
FVNavierStokesElemDisc<TDomain>::
assemble_JA(local_matrix_type& J, const local_vector_type& u)
{
// 	Only first order implementation
	UG_ASSERT((TFVGeom<TElem, dim>::order == 1), "Only first order implemented.");

// 	get finite volume geometry
	static const TFVGeom<TElem, dim>& geo = Provider<TFVGeom<TElem,dim> >::get();

//	check for source term to pass to the stabilization
	const DataImport<MathVector<dim>, dim>* pSource = NULL;
	if(m_imSource.data_given())	pSource = &m_imSource;

//	check for solutions to pass to stabilization in time-dependent case
	const local_vector_type *pSol = &u, *pOldSol = NULL;
	number dt = 0.0;
	if(this->is_time_dependent())
	{
	//	get and check current and old solution
		const LocalVectorTimeSeries* vLocSol = this->local_time_solutions();
		if(vLocSol->size() != 2)
		{
			UG_LOG("ERROR in 'FVNavierStokesElemDisc::assemble_A': "
					" Stabilization needs exactly two time points.\n");
			return false;
		}

	//	remember local solutions
		pSol = &vLocSol->solution(0);
		pOldSol = &vLocSol->solution(1);
		dt = vLocSol->time(0) - vLocSol->time(1);
	}

//	compute stabilized velocities and shapes for continuity equation
	if(!m_pStab->update(&geo, *pSol, m_imKinViscosity, pSource, pOldSol, dt))
	{
		UG_LOG("ERROR in 'FVNavierStokesElemDisc::assemble_A': "
				"Cannot compute stabilized velocities and shapes.\n");
		return false;
	}

//	compute stabilized velocities and shapes for convection upwind
	if(m_pConvStab != NULL)
		if(m_pConvStab != m_pStab)
			if(!m_pConvStab->update(&geo, *pSol, m_imKinViscosity, pSource, pOldSol, dt))
			{
				UG_LOG("ERROR in 'FVNavierStokesElemDisc::assemble_A': "
						"Cannot compute upwind (PAC) velocities and shapes.\n");
				return false;
			}

//	compute upwind shapes
	if(m_pConvUpwind != NULL)
		if(m_pStab->get_upwind() != m_pConvUpwind)
			if(!m_pConvUpwind->update(&geo, *pSol))
			{
				UG_LOG("ERROR in 'FVNavierStokesElemDisc::assemble_A': "
						"Cannot compute upwind velocities and shapes.\n");
				return false;
			}

//	get a const (!!) reference to the stabilization
	const INavierStokesStabilization<dim>& stab
		= *const_cast<const INavierStokesStabilization<dim>*>(m_pStab);

//	get a const (!!) reference to the stabilization of Convective Term
	const INavierStokesStabilization<dim>& convStab
		= *const_cast<const INavierStokesStabilization<dim>*>(m_pConvStab);

//	get a const (!!) reference to the upwind of Convective Term
	const INavierStokesUpwind<dim>& upwind
		= *const_cast<const INavierStokesUpwind<dim>*>(m_pConvUpwind);

// 	loop Sub Control Volume Faces (SCVF)
	for(size_t i = 0; i < geo.num_scvf(); ++i)
	{
	// 	get current SCVF
		const typename TFVGeom<TElem, dim>::SCVF& scvf = geo.scvf(i);

	// 	loop shape functions
		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
		{
			////////////////////////////////////////////////////
			////////////////////////////////////////////////////
			// Momentum Equation (conservation of momentum)
			////////////////////////////////////////////////////
			////////////////////////////////////////////////////

			////////////////////////////////////////////////////
			// Diffusive Term (Momentum Equation)
			////////////////////////////////////////////////////

		// 	Compute flux derivative at IP
			const number flux_sh =  -1.0 * m_imKinViscosity[i]
									* VecDot(scvf.global_grad(sh), scvf.normal());

		// 	Add flux derivative  to local matrix
			for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
			{
				J(d1, scvf.from(), d1, sh) += flux_sh;
				J(d1, scvf.to()  , d1, sh) -= flux_sh;

				for(size_t d2 = 0; d2 < (size_t)dim; ++d2)
				{
					const number flux2_sh = -1.0 * m_imKinViscosity[i]
											* scvf.global_grad(sh)[d1] * scvf.normal()[d2];
					J(d1, scvf.from(), d2, sh) += flux2_sh;
					J(d1, scvf.to()  , d2, sh) -= flux2_sh;
				}
			}

			////////////////////////////////////////////////////
			// Pressure Term (Momentum Equation)
			////////////////////////////////////////////////////

		//	Add flux derivative for local matrix
			for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
			{
				const number flux_sh = scvf.shape(sh) * scvf.normal()[d1];
				J(d1, scvf.from(), _P_, sh) += flux_sh;
				J(d1, scvf.to()  , _P_, sh) -= flux_sh;
			}

			////////////////////////////////////////////////////
			// Convective Term (Momentum Equation)
			////////////////////////////////////////////////////

		//	compute upwind velocity
			MathVector<dim> UpwindVel;

		//	switch PAC
			if(m_pConvUpwind != NULL)  UpwindVel = m_pConvUpwind->upwind_vel(i);
			else if (m_pConvStab != NULL) UpwindVel = convStab.stab_vel(i);
			else
			{
				UG_LOG("ERROR in 'FVNavierStokesElemDisc::assemble_A': "
						" Cannot find upwind for convective term.\n");
				return false;
			}

		//	peclet blend
			number w = 1.0;
			if(m_bPecletBlend)
				w = peclet_blend(UpwindVel, scvf, u, m_imKinViscosity[i]);

		//	compute product of stabilized vel and normal
			const number prod = VecProd(UpwindVel, scvf.normal());

		///////////////////////////////////
		//	Add fixpoint linearization
		///////////////////////////////////

		//	Stabilization used as upwind
			if(m_pConvStab != NULL)
			{
			//	velocity derivatives
				if(stab.vel_comp_connected())
					for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
						for(size_t d2 = 0; d2 < (size_t)dim; ++d2)
						{
							const number convFlux_vel = prod * w * convStab.stab_shape_vel(i, d1, d2, sh);
							J(d1, scvf.from(), d2, sh) += convFlux_vel;
							J(d1, scvf.to()  , d2, sh) -= convFlux_vel;
						}
				else
					for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
					{
						const number convFlux_vel = prod * w * convStab.stab_shape_vel(i, d1, d1, sh);
						J(d1, scvf.from(), d1, sh) += convFlux_vel;
						J(d1, scvf.to()  , d1, sh) -= convFlux_vel;
					}

			//	pressure derivative
				for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
				{
					const number convFlux_p = prod * w * convStab.stab_shape_p(i, d1, sh);

					J(d1, scvf.from(), _P_, sh) += convFlux_p;
					J(d1, scvf.to()  , _P_, sh) -= convFlux_p;
				}
			}

		//	Upwind used as upwind
			if(m_pConvUpwind != NULL)
			{
				const number convFlux_vel = prod * w * upwind.upwind_shape_sh(i, sh);
				for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
				{
					J(d1, scvf.from(), d1, sh) += convFlux_vel;
					J(d1, scvf.to()  , d1, sh) -= convFlux_vel;
				}
			}

		//	derivative due to peclet blending
			if(m_bPecletBlend)
			{
				const number convFluxPe = prod * (1.0-w) * scvf.shape(sh);
				for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
				{
					J(d1, scvf.from(), d1, sh) += convFluxPe;
					J(d1, scvf.to()  , d1, sh) -= convFluxPe;
				}
			}

		/////////////////////////////////////////
		//	Add full jacobian (remaining part)
		/////////////////////////////////////////

		//	Add remaining term for exact jacobian
			if(m_bExactJacobian)
			{
			//	Stabilization used as upwind
				if(m_pConvStab != NULL)
				{
				//	loop defect components
					for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
					{
						for(size_t d2 = 0; d2 < (size_t)dim; ++d2)
						{
					//	derivatives w.r.t. velocity
					//	Compute n * derivs
						number prod_vel = 0.0;

					//	Compute sum_j n_j * \partial_{u_i^sh} u_j
						if(stab.vel_comp_connected())
							for(size_t k = 0; k < (size_t)dim; ++k)
								prod_vel += w * convStab.stab_shape_vel(i, k, d2, sh)
												* scvf.normal()[k];
						else
							prod_vel = convStab.stab_shape_vel(i, d1, d1, sh)
												* scvf.normal()[d1];

						J(d1, scvf.from(), d2, sh) += prod_vel * UpwindVel[d1];
						J(d1, scvf.to()  , d2, sh) -= prod_vel * UpwindVel[d1];
						}

					//	derivative w.r.t pressure
					//	Compute n * derivs
						number prod_p = 0.0;

					//	Compute sum_j n_j * \parial_{u_i^sh} u_j
						for(size_t k = 0; k < (size_t)dim; ++k)
							prod_p += convStab.stab_shape_p(i, k, sh)
												* scvf.normal()[k];

						J(d1, scvf.from(), _P_, sh) += prod_p * UpwindVel[d1];
						J(d1, scvf.to()  , _P_, sh) -= prod_p * UpwindVel[d1];
					}
				}

			//	Upwind used as upwind
				if(m_pConvUpwind != NULL)
				{
				//	loop defect components
					for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
					{
						for(size_t d2 = 0; d2 < (size_t)dim; ++d2)
						{
					//	derivatives w.r.t. velocity
						number prod_vel = w * upwind.upwind_shape_sh(i,sh)	* scvf.normal()[d2];

						J(d1, scvf.from(), d2, sh) += prod_vel * UpwindVel[d1];
						J(d1, scvf.to()  , d2, sh) -= prod_vel * UpwindVel[d1];
						}
					}
				}

			//	derivative due to peclet blending
				if(m_bPecletBlend)
				{
					for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
						for(size_t d2 = 0; d2 < (size_t)dim; ++d2)
						{
							const number convFluxPe = UpwindVel[d1] * (1.0-w)
													  * scvf.shape(sh)
													  * scvf.normal()[d2];
							J(d1, scvf.from(), d2, sh) += convFluxPe;
							J(d1, scvf.to()  , d2, sh) -= convFluxPe;
						}
				}
			} // end exact jacobian part

			////////////////////////////////////////////////////
			////////////////////////////////////////////////////
			// Continuity Equation (conservation of mass)
			////////////////////////////////////////////////////
			////////////////////////////////////////////////////

		//	Add derivative of stabilized flux w.r.t velocity comp to local matrix
			if(stab.vel_comp_connected())
			{
				for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
				{
					number contFlux_vel = 0.0;
					for(size_t d2 = 0; d2 < (size_t)dim; ++d2)
						contFlux_vel += stab.stab_shape_vel(i, d2, d1, sh)
										* scvf.normal()[d2];

					J(_P_, scvf.from(), d1, sh) += contFlux_vel;
					J(_P_, scvf.to()  , d1, sh) -= contFlux_vel;
				}
			}
			else
			{
				for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
				{
					const number contFlux_vel = stab.stab_shape_vel(i, d1, d1, sh)
													* scvf.normal()[d1];

					J(_P_, scvf.from(), d1, sh) += contFlux_vel;
					J(_P_, scvf.to()  , d1, sh) -= contFlux_vel;
				}
			}

		//	Add derivative of stabilized flux w.r.t pressure to local matrix
			number contFlux_p = 0.0;
			for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
				contFlux_p += stab.stab_shape_p(i, d1, sh) * scvf.normal()[d1];

			J(_P_, scvf.from(), _P_, sh) += contFlux_p;
			J(_P_, scvf.to()  , _P_, sh) -= contFlux_p;
		}
	}

	// we're done
	return true;
}

template<typename TDomain>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
bool
FVNavierStokesElemDisc<TDomain>::
assemble_A(local_vector_type& d, const local_vector_type& u)
{
// 	Only first order implemented
	UG_ASSERT((TFVGeom<TElem, dim>::order == 1), "Only first order implemented.");

// 	get finite volume geometry
	static const TFVGeom<TElem, dim>& geo = Provider<TFVGeom<TElem,dim> >::get();

//	check for source term to pass to the stabilization
	const DataImport<MathVector<dim>, dim>* pSource = NULL;
	if(m_imSource.data_given())	pSource = &m_imSource;

//	check for solutions to pass to stabilization in time-dependent case
	const local_vector_type *pSol = &u, *pOldSol = NULL;
	number dt = 0.0;
	if(this->is_time_dependent())
	{
	//	get and check current and old solution
		const LocalVectorTimeSeries* vLocSol = this->local_time_solutions();
		if(vLocSol->size() != 2)
		{
			UG_LOG("ERROR in 'FVNavierStokesElemDisc::assemble_A': "
					" Stabilization needs exactly two time points.\n");
			return false;
		}

	//	remember local solutions
		pSol = &vLocSol->solution(0);
		pOldSol = &vLocSol->solution(1);
		dt = vLocSol->time(0) - vLocSol->time(1);
	}

//	compute stabilized velocities and shapes for continuity equation
	// \todo: (optional) Here we can skip the computation of shapes, implement?
	if(!m_pStab->update(&geo, *pSol, m_imKinViscosity, pSource, pOldSol, dt))
	{
		UG_LOG("ERROR in 'FVNavierStokesElemDisc::assemble_A': "
				"Cannot compute stabilized velocities and shapes.\n");
		return false;
	}

//	compute stabilized velocities and shapes for convection upwind
	if(m_pConvStab != NULL)
		if(m_pConvStab != m_pStab)
			if(!m_pConvStab->update(&geo, *pSol, m_imKinViscosity, pSource, pOldSol, dt))
			{
				UG_LOG("ERROR in 'FVNavierStokesElemDisc::assemble_A': "
						"Cannot compute upwind (PAC) velocities and shapes.\n");
				return false;
			}

//	compute upwind shapes
	if(m_pConvUpwind != NULL)
		if(m_pStab->get_upwind() != m_pConvUpwind)
			if(!m_pConvUpwind->update(&geo, *pSol))
			{
				UG_LOG("ERROR in 'FVNavierStokesElemDisc::assemble_A': "
						"Cannot compute upwind velocities and shapes.\n");
				return false;
			}

//	get a const (!!) reference to the stabilization
	const INavierStokesStabilization<dim>& stab
		= *const_cast<const INavierStokesStabilization<dim>*>(m_pStab);

//	get a const (!!) reference to the stabilization of Convective Term
	const INavierStokesStabilization<dim>& convStab
		= *const_cast<const INavierStokesStabilization<dim>*>(m_pConvStab);

// 	loop Sub Control Volume Faces (SCVF)
	for(size_t i = 0; i < geo.num_scvf(); ++i)
	{
	// 	get current SCVF
		const typename TFVGeom<TElem, dim>::SCVF& scvf = geo.scvf(i);

		////////////////////////////////////////////////////
		////////////////////////////////////////////////////
		// Momentum Equation (conservation of momentum)
		////////////////////////////////////////////////////
		////////////////////////////////////////////////////

	//	\todo: we could also add all fluxes at once in order to save several
	//			accesses to the local defect, implement and loose clarity?

		////////////////////////////////////////////////////
		// Diffusive Term (Momentum Equation)
		////////////////////////////////////////////////////

	// 	1. Interpolate Functional Matrix of velocity at ip
		MathMatrix<dim, dim> gradVel;
		for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
			for(size_t d2 = 0; d2 < (size_t)dim; ++d2)
			{
			//	sum up contributions of each shape
				gradVel(d1, d2) = 0.0;
				for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
					gradVel(d1, d2) += scvf.global_grad(sh)[d2]
					                    * u(d1, sh);
			}

	//	2. Compute flux
		MathVector<dim> diffFlux;

	//	Add (\nabla u) \cdot \vec{n}
		MatVecMult(diffFlux, gradVel, scvf.normal());

	//	Add (\nabla u)^T \cdot \vec{n}
		TransposedMatVecMultAdd(diffFlux, gradVel, scvf.normal());

	//	scale by viscosity
		VecScale(diffFlux, diffFlux, (-1.0) * m_imKinViscosity[i]);

	//	3. Add flux to local defect
		for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
		{
			d(d1, scvf.from()) += diffFlux[d1];
			d(d1, scvf.to()  ) -= diffFlux[d1];
		}

		////////////////////////////////////////////////////
		// Convective Term (Momentum Equation)
		////////////////////////////////////////////////////

	//	find the upwind velocity at ip
		MathVector<dim> UpwindVel;

	//	switch PAC
		if(m_pConvUpwind != NULL)  UpwindVel = m_pConvUpwind->upwind_vel(i);
		else if (m_pConvStab != NULL) UpwindVel = convStab.stab_vel(i);
		else
		{
			UG_LOG("ERROR in 'FVNavierStokesElemDisc::assemble_A': "
					" Cannot find upwind for convective term.\n");
			return false;
		}

	//	Peclet Blend
		if(m_bPecletBlend)
			peclet_blend(UpwindVel, scvf, u, m_imKinViscosity[i]);


//	\todo: implement the linearization in a clean way
        MathVector<dim> StdVel; VecSet(StdVel, 0.0);
        for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
            for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
                StdVel[d1] += u(d1, sh) * scvf.shape(sh);


	//	compute product of upwinded (and blended) velocity and normal
		//const number prod = VecProd(UpwindVel, scvf.normal());
		const number prod = VecProd(StdVel, scvf.normal());

	//	Add contributions to local velocity components
		for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
		{
			d(d1, scvf.from()) += UpwindVel[d1] * prod;
			d(d1, scvf.to()  ) -= UpwindVel[d1] * prod;
		}

		////////////////////////////////////////////////////
		// Pressure Term (Momentum Equation)
		////////////////////////////////////////////////////

	//	1. Interpolate pressure at ip
		number pressure = 0.0;
		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
			pressure += scvf.shape(sh) * u(_P_, sh);

	//	2. Add contributions to local defect
		for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
		{
			d(d1, scvf.from()) += pressure * scvf.normal()[d1];
			d(d1, scvf.to()  ) -= pressure * scvf.normal()[d1];
		}

		////////////////////////////////////////////////////
		////////////////////////////////////////////////////
		// Continuity Equation (conservation of mass)
		////////////////////////////////////////////////////
		////////////////////////////////////////////////////

	//	compute flux at ip
		const number contFlux = VecProd(stab.stab_vel(i), scvf.normal());

	//	Add contributions to local defect
		d(_P_, scvf.from()) += contFlux;
		d(_P_, scvf.to()  ) -= contFlux;
	}

// 	we're done
	return true;
}


template<typename TDomain>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
bool
FVNavierStokesElemDisc<TDomain>::
assemble_JM(local_matrix_type& J, const local_vector_type& u)
{
// 	Only first order implementation
	UG_ASSERT((TFVGeom<TElem, dim>::order == 1), "Only first order implemented.");

// 	get finite volume geometry
	const static TFVGeom<TElem, dim>& geo = Provider<TFVGeom<TElem,dim> >::get();

// 	loop Sub Control Volumes (SCV)
	for(size_t i = 0; i < geo.num_scv(); ++i)
	{
	// 	get SCV
		const typename TFVGeom<TElem, dim>::SCV& scv = geo.scv(i);

	// 	get associated node
		const int sh = scv.node_id();

	// 	loop velocity components
		for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
		{
		// 	Add to local matrix
			J(d1, sh, d1, sh) += scv.volume();
		}
	}

// 	we're done
	return true;
}


template<typename TDomain>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
bool
FVNavierStokesElemDisc<TDomain>::
assemble_M(local_vector_type& d, const local_vector_type& u)
{
// 	Only first order implementation
	UG_ASSERT((TFVGeom<TElem, dim>::order == 1), "Only first order implemented.");

// 	get finite volume geometry
	const static TFVGeom<TElem, dim>& geo = Provider<TFVGeom<TElem,dim> >::get();

// 	loop Sub Control Volumes (SCV)
	for(size_t i = 0; i < geo.num_scv(); ++i)
	{
	// 	get current SCV
		const typename TFVGeom<TElem, dim>::SCV& scv = geo.scv(i);

	// 	get associated node
		const int sh = scv.node_id();

	// 	loop velocity components
		for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
		{
		// 	Add to local matrix
			d(d1, sh) += u(d1, sh) * scv.volume();
		}
	}

// 	we're done
	return true;
}


template<typename TDomain>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
bool
FVNavierStokesElemDisc<TDomain>::
assemble_f(local_vector_type& d)
{
// 	Only first order implementation
	UG_ASSERT((TFVGeom<TElem, dim>::order == 1), "Only first order implemented.");

//	if zero data given, return
	if(!m_imSource.data_given()) return true;

// 	get finite volume geometry
	const static TFVGeom<TElem, dim>& geo = Provider<TFVGeom<TElem,dim> >::get();

// 	loop Sub Control Volumes (SCV)
	for(size_t ep = 0; ep < geo.num_scv(); ++ep)
	{
	// 	get current SCV
		const typename TFVGeom<TElem, dim>::SCV& scv = geo.scv(ep);

	// 	get associated node
		const int sh = scv.node_id();

	// 	Add to local rhs
		for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
			d(d1, sh) += m_imSource[ep][d1] * scv.volume();
	}

// 	we're done
	return true;
}

template<typename TDomain>
template <typename SCVF>
inline
number
FVNavierStokesElemDisc<TDomain>::
peclet_blend(MathVector<dim>& UpwindVel, const SCVF& scvf,
             const local_vector_type& u, number kinVisco)
{
//	Interpolate velocity at ip
	MathVector<dim> InterpolVel; VecSet(InterpolVel, 0.0);
	for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
			InterpolVel[d1] += scvf.shape(sh) * u(d1, sh);

//	compute peclet number
	number Pe = VecProd(InterpolVel, scvf.normal());
	// \todo: ATTENTION: scvf.global_corner(0) gives edge midpoint of corresponding
	//					 edge, but this is not defined as an interface specification,
	//					 but only a feature of the special implementation (specify it)
	Pe *= VecDistance(scvf.global_ip(), scvf.global_corner(0));
	Pe /= kinVisco;

//	compute weight
	const number Pe2 = Pe * Pe;
	const number w = Pe2 / (5.0 + Pe2);

//	compute upwind vel
	VecScaleAdd(UpwindVel, w, UpwindVel, (1.0-w), InterpolVel);

	return w;
}


////////////////////////////////////////////////////////////////////////////////
//	Constructor
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
FVNavierStokesElemDisc<TDomain>::FVNavierStokesElemDisc()
: m_bExactJacobian(true), m_pStab(NULL),
  m_pConvStab(NULL), m_pConvUpwind(NULL)
{
//	set default options
	set_peclet_blend(false);

//	register imports
	register_import(m_imSource);
	register_import(m_imKinViscosity);

//	register assemble functions
	register_all_fv1_funcs(false);
}


////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////

// register for 1D
template<typename TDomain>
void
FVNavierStokesElemDisc<TDomain>::
register_all_fv1_funcs(bool bHang)
{
//	get all grid element types in this dimension and below
	typedef typename domain_traits<dim>::AllElemList ElemList;

//	switch assemble functions
	if(!bHang) boost::mpl::for_each<ElemList>( RegisterFV1<FV1Geometry>(this) );
	else throw(UGFatalError("Not implemented."));
}

template<typename TDomain>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
void
FVNavierStokesElemDisc<TDomain>::
register_fv1_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef this_type T;

	set_prep_elem_loop_fct(id, &T::template prepare_element_loop<TElem, TFVGeom>);
	set_prep_elem_fct(	 id, &T::template prepare_element<TElem, TFVGeom>);
	set_fsh_elem_loop_fct( id, &T::template finish_element_loop<TElem, TFVGeom>);
	set_ass_JA_elem_fct(		 id, &T::template assemble_JA<TElem, TFVGeom>);
	set_ass_JM_elem_fct(		 id, &T::template assemble_JM<TElem, TFVGeom>);
	set_ass_dA_elem_fct(		 id, &T::template assemble_A<TElem, TFVGeom>);
	set_ass_dM_elem_fct(		 id, &T::template assemble_M<TElem, TFVGeom>);
	set_ass_rhs_elem_fct(	 id, &T::template assemble_f<TElem, TFVGeom>);
}

////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

template class FVNavierStokesElemDisc<Domain1d>;
template class FVNavierStokesElemDisc<Domain2d>;
template class FVNavierStokesElemDisc<Domain3d>;


} // namespace ug
