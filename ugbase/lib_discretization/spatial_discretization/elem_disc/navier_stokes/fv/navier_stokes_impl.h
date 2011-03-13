/*
 * navier_stokes_impl.h
 *
 *  Created on: 20.09.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__NAVIER_STOKES__FV__NAVIER_STOKES_IMPL__
#define __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__NAVIER_STOKES__FV__NAVIER_STOKES_IMPL__

namespace ug{


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  Provide a generic implementation for all elements
//  (since this discretization can be implemented in a generic way)
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain, typename TAlgebra>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
bool
FVNavierStokesElemDisc<TDomain, TAlgebra>::
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
		TFVGeom<TElem, dim>& geo = FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();
		m_imKinViscosity.template set_local_ips<refDim>(geo.scvf_local_ips(),
		                                                geo.num_scvf_local_ips());
		m_imSource.template set_local_ips<refDim>(geo.scv_local_ips(),
		                                          geo.num_scv_local_ips());
	}

//	we're done
	return true;
}

template<typename TDomain, typename TAlgebra>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
bool
FVNavierStokesElemDisc<TDomain, TAlgebra>::
finish_element_loop()
{
//	nothing to do
	return true;
}

template<typename TDomain, typename TAlgebra>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
bool
FVNavierStokesElemDisc<TDomain, TAlgebra>::
prepare_element(TElem* elem, const local_vector_type& u,
								const local_index_type& glob_ind)
{
//	get corners
	m_vCornerCoords = this->template get_element_corners<TElem>(elem);

// 	Update Geometry for this element
	TFVGeom<TElem, dim>& geo = FVGeometryProvider::get_geom<TFVGeom, TElem, dim>();
	if(!geo.update(elem, this->get_subset_handler(), &m_vCornerCoords[0]))
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
		                                                geo.num_scvf_local_ips());
		m_imSource.template set_local_ips<refDim>(geo.scv_local_ips(),
		                                          geo.num_scv_local_ips());
	}

//	set global positions for imports
	m_imKinViscosity.set_global_ips(geo.scvf_global_ips(), geo.num_scvf_global_ips());
	m_imSource.set_global_ips(geo.scv_global_ips(), geo.num_scv_global_ips());

//	we're done
	return true;
}

template<typename TDomain, typename TAlgebra>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
bool
FVNavierStokesElemDisc<TDomain, TAlgebra>::
assemble_JA(local_matrix_type& J, const local_vector_type& u, number time)
{
// 	Only first order implementation
	UG_ASSERT((TFVGeom<TElem, dim>::order == 1), "Only first order implemented.");

//	for first order only one integration point on scvf
	static const size_t ip = 0;

// 	get finite volume geometry
	static const TFVGeom<TElem, dim>& geo = FVGeometryProvider::get_geom<TFVGeom, TElem, dim>();

//	check for source term to pass to the stabilization
	const DataImport<MathVector<dim>, dim, TAlgebra>* pSource = NULL;
	if(m_imSource.data_given())	pSource = &m_imSource;

//	compute stabilized velocities and shapes
	// \todo: handle time dependent case (i.e. not passing NULL)
	if(!m_pStab->update(&geo, u, m_imKinViscosity, pSource, NULL, 0.0))
	{
		UG_LOG("ERROR in 'FVNavierStokesElemDisc::assemble_A': "
				"Cannot compute stabilized velocities and shapes.\n");
		return false;
	}

//	get a const (!!) reference to the stabilization
	const INavierStokesStabilization<dim, algebra_type>& stab
		= *const_cast<const INavierStokesStabilization<dim, algebra_type>*>(m_pStab);

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
									* VecDot(scvf.global_grad(sh, ip), scvf.normal());

		// 	Add flux derivative  to local matrix
			for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
			{
				J(d1, scvf.from(), d1, sh) += flux_sh;
				J(d1, scvf.to()  , d1, sh) -= flux_sh;

				for(size_t d2 = 0; d2 < (size_t)dim; ++d2)
				{
					const number flux2_sh = -1.0 * m_imKinViscosity[i]
											* scvf.global_grad(sh, ip)[d1] * scvf.normal()[d2];
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
				const number flux_sh = scvf.shape(sh, ip) * scvf.normal()[d1];
				J(d1, scvf.from(), _P_, sh) += flux_sh;
				J(d1, scvf.to()  , _P_, sh) -= flux_sh;
			}

			////////////////////////////////////////////////////
			// Convective Term (Momentum Equation)
			////////////////////////////////////////////////////

		//	compute upwind velocity
			MathVector<dim> UpwindVel;

		//	switch PAC
			if(!m_bPAC)
			{
				UG_LOG("ERROR: NOPAC not implemented.\n"); return false;
			}
			else UpwindVel = stab.stab_vel(i);

		//	peclet blend
			number w = 1.0;
			if(m_bPecletBlend)
				w = peclet_blend(UpwindVel, scvf, u, m_imKinViscosity[i]);

		//	compute product of stabilized vel and normal
			const number prod = VecProd(UpwindVel, scvf.normal());

		//	Add fixpoint linearization
		//	velocity derivatives
			if(stab.vel_comp_connected())
				for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
					for(size_t d2 = 0; d2 < (size_t)dim; ++d2)
					{
						const number convFlux_vel = prod * w * stab.stab_shape_vel(i, d1, d2, sh);
						J(d1, scvf.from(), d2, sh) += convFlux_vel;
						J(d1, scvf.to()  , d2, sh) -= convFlux_vel;
					}
			else
				for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
				{
					const number convFlux_vel = prod * w * stab.stab_shape_vel(i, d1, d1, sh);
					J(d1, scvf.from(), d1, sh) += convFlux_vel;
					J(d1, scvf.to()  , d1, sh) -= convFlux_vel;
				}

		//	pressure derivative
			for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
			{
				const number convFlux_p = prod * w * stab.stab_shape_p(i, d1, sh);

				J(d1, scvf.from(), _P_, sh) += convFlux_p * UpwindVel[d1];
				J(d1, scvf.to()  , _P_, sh) -= convFlux_p * UpwindVel[d1];
			}

		//	derivative due to peclet blending
			if(m_bPecletBlend)
			{
				const number convFluxPe = prod * (1.0-w) * scvf.shape(sh, ip);
				for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
				{
					J(d1, scvf.from(), d1, sh) += convFluxPe;
					J(d1, scvf.to()  , d1, sh) -= convFluxPe;
				}
			}

		//	Add remaining term for exact jacobian
			if(m_bExactJacobian)
			{
			//	loop defect components
				for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
				{
				//	derivaztives w.r.t. velocity
				//	Compute n * derivs
					number prod_vel = 0.0;

				//	Compute sum_j n_j * \parial_{u_i^sh} u_j
					if(stab.vel_comp_connected())
						for(size_t d2 = 0; d2 < (size_t)dim; ++d2)
							prod_vel += stab.stab_shape_vel(i, d2, d1, sh)
											* scvf.normal()[d2];
					else
						prod_vel = stab.stab_shape_vel(i, d1, d1, sh)
											* scvf.normal()[d1];

					for(size_t d2 = 0; d2 < (size_t)dim; ++d2)
					{
						J(d1, scvf.from(), d2, sh) += prod_vel * UpwindVel[d1];
						J(d1, scvf.to()  , d2, sh) -= prod_vel * UpwindVel[d1];
					}

				//	derivative w.r.t pressure
				//	Compute n * derivs
					number prod_p = 0.0;

				//	Compute sum_j n_j * \parial_{u_i^sh} u_j
					for(size_t d2 = 0; d2 < (size_t)dim; ++d2)
						prod_p += stab.stab_shape_p(i, d2, sh)
											* scvf.normal()[d2];

					J(d1, scvf.from(), _P_, sh) += prod_p * UpwindVel[d1];
					J(d1, scvf.to()  , _P_, sh) -= prod_p * UpwindVel[d1];

				//	derivative due to peclet blending
					if(m_bPecletBlend)
					{
						for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
							for(size_t d2 = 0; d2 < (size_t)dim; ++d2)
							{
								const number convFluxPe = UpwindVel[d1] * (1.0-w)
														  * scvf.shape(sh, ip)
														  * scvf.normal()[d2];
								J(d1, scvf.from(), d2, sh) += convFluxPe;
								J(d1, scvf.to()  , d2, sh) -= convFluxPe;
							}
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

template<typename TDomain, typename TAlgebra>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
bool
FVNavierStokesElemDisc<TDomain, TAlgebra>::
assemble_A(local_vector_type& d, const local_vector_type& u, number time)
{
// 	Only first order implemented
	UG_ASSERT((TFVGeom<TElem, dim>::order == 1), "Only first order implemented.");

//	for first order only one integration point on scvf
	static const size_t ip = 0;

// 	get finite volume geometry
	static const TFVGeom<TElem, dim>& geo = FVGeometryProvider::get_geom<TFVGeom, TElem, dim>();

//	check for source term to pass to the stabilization
	const DataImport<MathVector<dim>, dim, TAlgebra>* pSource = NULL;
	if(m_imSource.data_given())	pSource = &m_imSource;

//	compute stabilized velocities and shapes
	// \todo: (optional) Here we can skip the computation of shapes, implement?
	// \todo: handle time dependent case (i.e. not passing NULL)
	if(!m_pStab->update(&geo, u, m_imKinViscosity, pSource, NULL, 0.0))
	{
		UG_LOG("ERROR in 'FVNavierStokesElemDisc::assemble_A': "
				"Cannot compute stabilized velocities and shapes.\n");
		return false;
	}

//	get a const (!!) reference to the stabilization
	const INavierStokesStabilization<dim, algebra_type>& stab
		= *const_cast<const INavierStokesStabilization<dim, algebra_type>*>(m_pStab);

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
					gradVel(d1, d2) += scvf.global_grad(sh, ip)[d2]
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
		if(!m_bPAC)
		{
			UG_LOG("ERROR: NOPAC not implemented.\n"); return false;
		}
		else UpwindVel = stab.stab_vel(i);

	//	Peclet Blend
		if(m_bPecletBlend)
			peclet_blend(UpwindVel, scvf, u, m_imKinViscosity[i]);

	//	compute product of upwinded (and blended) velocity and normal
		const number prod = VecProd(UpwindVel, scvf.normal());

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
			pressure += scvf.shape(sh, ip) * u(_P_, sh);

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


template<typename TDomain, typename TAlgebra>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
bool
FVNavierStokesElemDisc<TDomain, TAlgebra>::
assemble_JM(local_matrix_type& J, const local_vector_type& u, number time)
{
// 	Only first order implementation
	UG_ASSERT((TFVGeom<TElem, dim>::order == 1), "Only first order implemented.");

// 	get finite volume geometry
	static TFVGeom<TElem, dim>& geo = FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();

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


template<typename TDomain, typename TAlgebra>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
bool
FVNavierStokesElemDisc<TDomain, TAlgebra>::
assemble_M(local_vector_type& d, const local_vector_type& u, number time)
{
// 	Only first order implementation
	UG_ASSERT((TFVGeom<TElem, dim>::order == 1), "Only first order implemented.");

// 	get finite volume geometry
	static TFVGeom<TElem, dim>& geo = FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();

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


template<typename TDomain, typename TAlgebra>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
bool
FVNavierStokesElemDisc<TDomain, TAlgebra>::
assemble_f(local_vector_type& d, number time)
{
// 	Only first order implementation
	UG_ASSERT((TFVGeom<TElem, dim>::order == 1), "Only first order implemented.");

//	if zero data given, return
	if(!m_imSource.data_given()) return true;

// 	get finite volume geometry
	static TFVGeom<TElem, dim>& geo
		= FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();

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

template<typename TDomain, typename TAlgebra>
template <typename SCVF>
inline
number
FVNavierStokesElemDisc<TDomain, TAlgebra>::
peclet_blend(MathVector<dim>& UpwindVel, const SCVF& scvf,
             const local_vector_type& u, number kinVisco)
{
//	for first order only one integration point on scvf
	static const size_t ip = 0;

//	Interpolate velocity at ip
	MathVector<dim> InterpolVel; VecSet(InterpolVel, 0.0);
	for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
			InterpolVel[d1] += scvf.shape(sh, ip) * u(d1, sh);

//	compute peclet number
	number Pe = VecProd(InterpolVel, scvf.normal());
	// \todo: ATTENTION: scvf.global_corner(0) gives edge midpoint of corresponding
	//					 edge, but this is not defined as an interface specification,
	//					 but only a feature of the special implementation (specify it)
	Pe *= VecDistance(scvf.global_ip(ip), scvf.global_corner(0));
	Pe /= kinVisco;

//	compute weight
	const number Pe2 = Pe * Pe;
	const number w = Pe2 / (5.0 + Pe2);

//	compute upwind vel
	VecScaleAdd(UpwindVel, w, UpwindVel, (1.0-w), InterpolVel);

	return w;
}

} // namespace ug


#endif /*__H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__NAVIER_STOKES__FV__NAVIER_STOKES_IMPL__*/
