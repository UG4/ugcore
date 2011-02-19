/*
 * navier_stokes_impl.h
 *
 *  Created on: 20.09.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__NAVIER_STOKES__FV__NAVIER_STOKES_IMPL__
#define __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__NAVIER_STOKES__FV__NAVIER_STOKES_IMPL__

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
template<typename TElem >
inline
bool
FVNavierStokesElemDisc<TFVGeom, TDomain, TAlgebra>::
prepare_element_loop()
{
	// all this will be performed outside of the loop over the elements.
	// Therefore it is not time critical.

	// resize corner coordinates
	typedef typename reference_element_traits<TElem>::reference_element_type
																ref_elem_type;
	m_vCornerCoords.resize(ref_elem_type::num_corners);

	// remember position attachement
	if(m_pDomain == NULL)
	{
		UG_LOG("ERROR in 'FVNavierStokesElemDisc::prepare_element_loop':"
				" Domain not set.");
		return false;
	}
	m_aaPos = m_pDomain->get_position_accessor();

	return true;
}

template<	template <class TElem, int TWorldDim> class TFVGeom,
			typename TDomain,
			typename TAlgebra>
template<typename TElem >
inline
bool
FVNavierStokesElemDisc<TFVGeom, TDomain, TAlgebra>::
finish_element_loop()
{
	// all this will be performed outside of the loop over the elements.
	// Therefore it is not time critical.

	return true;
}

template<	template <class TElem, int TWorldDim> class TFVGeom,
			typename TDomain,
			typename TAlgebra>
template<typename TElem >
inline
bool
FVNavierStokesElemDisc<TFVGeom, TDomain, TAlgebra>::
prepare_element(TElem* elem, const local_vector_type& u,
								const local_index_type& glob_ind)
{
	// this loop will be performed inside the loop over the elements.
	// Therefore, it is TIME CRITICAL

// 	Load corners of this element
	for(size_t i = 0; i < m_vCornerCoords.size(); ++i)
	{
		VertexBase* vert = elem->vertex(i);
		m_vCornerCoords[i] = m_aaPos[vert];
	}

// 	Update Geometry for this element
	TFVGeom<TElem, dim>& geo = FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();
	if(!geo.update(elem, m_pDomain->get_subset_handler(), &m_vCornerCoords[0]))
	{
		UG_LOG("FVNavierStokesElemDisc::prepare_element:"
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
FVNavierStokesElemDisc<TFVGeom, TDomain, TAlgebra>::
assemble_JA(local_matrix_type& J, const local_vector_type& u, number time)
{
//	number of corners == number of scv
	static const size_t numCo = TFVGeom<TElem, dim>::m_numSCV;

//	number of SubControlVolumeFaces == num ips
	static const size_t numSCVF = TFVGeom<TElem, dim>::m_numSCVF;

	// get finite volume geometry
	static const TFVGeom<TElem, dim>& geo = FVGeometryProvider::get_geom<TFVGeom, TElem, dim>();

	// Only first order implemented
	UG_ASSERT((TFVGeom<TElem, dim>::order == 1), "Only first order implemented.");

	// Some Variables

    MathVector<dim> vCornerVels[numSCVF];
	MathVector<dim> vCornerVelsOld[numSCVF];
    number          vCornerPress[numCo];
    number          vConvLength[numSCVF];
    // MathVector<dim> vIPVelUpwindShapesMomEq[numSCVF][numCo][dim];
    number vIPVelUpwindShapesContiEq[numSCVF][numCo];
    number vIPVelUpwindDependencies[numSCVF][numSCVF];
    MathVector<dim> vIPStabVelShapesContiEq[numSCVF][numCo][dim+1];

   	// Store Corner values in vCornerVels variable
    for(size_t co = 0; co < numCo; ++co)
        for (size_t component = 0; component < dim; ++component)
            vCornerVels[co][component]=u(component,co);

	// Compute Upwind Shapes at Ip's and ConvectionLength here fot the Momentum Equation
    GetUpwindShapes(geo, vCornerVels, m_UpwindMethod, vIPVelUpwindShapesContiEq,vIPVelUpwindDependencies, vConvLength);

	// todo: Implement for timedependent
	bool bTimeDependent = false;
	number dt = 0.0;

	// todo: switch
	// Compute Stabilized Velocities at IP's here (depending on Upwind Velocities)
	GetStabilizedShapes(	geo, vCornerVels, vCornerPress, m_StabOptions, vIPVelUpwindShapesContiEq,vIPVelUpwindDependencies, vConvLength,
									dt, bTimeDependent, vCornerVelsOld,
									m_Viscosity,
                                    vIPStabVelShapesContiEq);

	// loop Sub Control Volume Faces (SCVF)
	for(size_t i = 0; i < geo.num_scvf(); ++i)
	{
		// get current SCVF
		const typename TFVGeom<TElem, dim>::SCVF& scvf = geo.scvf(i);

		// loop integration point of SCVF
		for(size_t ip = 0; ip < scvf.num_ip(); ++ip)
		{
			// loop shape functions
			for(size_t co = 0; co < scvf.num_sh(); ++co)
			{
				////////////////////////////////////////////////////
				// Diffusiv Term (Momentum Equation)
				////////////////////////////////////////////////////

				// Compute flux at IP
				number flux = VecDot(scvf.global_grad(co, ip), scvf.normal());

				// Add flux term to local matrix
				for(size_t vel1 = 0; vel1 < dim; ++vel1)
				{
					J(vel1, scvf.from(), vel1, co) -= m_Viscosity * flux;
					J(vel1, scvf.to()  , vel1, co) += m_Viscosity * flux;

					for(size_t vel2 = 0; vel2 < dim; ++vel2)
					{
						const number flux2 = m_Viscosity * scvf.global_grad(co, ip)[vel1] * scvf.normal()[vel2];
						J(vel1, scvf.from(), vel2, co) -= flux2;
						J(vel1, scvf.to()  , vel2, co) += flux2;
					}
				}

				////////////////////////////////////////////////////
				// Pressure Term (Momentum Equation)
				////////////////////////////////////////////////////

                for(size_t vel = 0; vel < dim; ++vel)
				{
					J(vel, scvf.from(), _P_, co) += scvf.shape(co, ip) * scvf.normal()[vel];
					J(vel, scvf.to()  , _P_, co) -= scvf.shape(co, ip) * scvf.normal()[vel];
				}

				////////////////////////////////////////////////////
				// Convective Term (Momentum Equation)
				////////////////////////////////////////////////////

				// Compute flux at IP
				flux = VecDot(scvf.global_grad(co, ip), scvf.normal());

                // Add flux term to local matrix
				for(size_t vel1 = 0; vel1 < dim; ++vel1)
				{
					J(vel1, scvf.from(), vel1, co) += m_Viscosity * flux;
					J(vel1, scvf.to()  , vel1, co) -= m_Viscosity * flux;

					for(size_t vel2 = 0; vel2 < dim; ++vel2)
					{
						const number flux2 = m_Viscosity * scvf.global_grad(co, ip)[vel1] * scvf.normal()[vel2];
						J(vel1, scvf.from(), vel2, co) += flux2;
						J(vel1, scvf.to()  , vel2, co) -= flux2;
					}
				}

				// todo: implement convective term (momentum equation) use vIPVelStab

				// todo: implement mass equations use vIPVelStabMassEq
                //J(_P_, ,vel1, ) +=
                }
          }
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
FVNavierStokesElemDisc<TFVGeom, TDomain, TAlgebra>::
assemble_JM(local_matrix_type& J, const local_vector_type& u, number time)
{
	// get finite volume geometry
	TFVGeom<TElem, dim>& geo = FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();

	// loop Sub Control Volumes (SCV)
	for(size_t i = 0; i < geo.num_scv(); ++i)
	{
		// get current SCV
		const typename TFVGeom<TElem, dim>::SCV& scv = geo.scv(i);

		// get associated node
		const int co = scv.node_id();

		// loop velocity components
		for(size_t vel1 = 0; vel1 < dim; ++vel1)
		{
			// Add to local matrix
			J(vel1, co, vel1, co) += scv.volume();
		}
	}

	// we're done
	return true;
}


template<template <class TElem, int TWorldDim> class TFVGeom, typename TDomain, typename TAlgebra>
template<typename TElem >
inline
bool
FVNavierStokesElemDisc<TFVGeom, TDomain, TAlgebra>::
assemble_A(local_vector_type& d, const local_vector_type& u, number time)
{
//	number of corners == number of scv
	static const size_t numCo = TFVGeom<TElem, dim>::m_numSCV;

//	number of SubControlVolumeFaces == num ips
	static const size_t numSCVF = TFVGeom<TElem, dim>::m_numSCVF;

	// get finite volume geometry
	static const TFVGeom<TElem, dim>& geo = FVGeometryProvider::get_geom<TFVGeom, TElem, dim>();

	// Only first order implemented
	UG_ASSERT((TFVGeom<TElem, dim>::order == 1), "Only first order implemented.");

	// Some Variables

    MathVector<dim> vCornerVels[numCo];
	MathVector<dim> vCornerVelsOld[numCo];
    number          vCornerPress[numCo];
    number          vConvLength[numSCVF];
//   number vIPVelUpwindShapesMomEq[numSCVF][numCo];
    number vIPVelUpwindShapesContiEq[numSCVF][numCo];
    number vIPVelUpwindDependenciesContiEq[numSCVF][numSCVF];
    MathVector<dim> vIPStabVelShapesContiEq[numSCVF][numCo][dim+1];

   	// Store Corner values in vCornerVels variable
    for(size_t co = 0; co < numCo; ++co){
        for (size_t component = 0; component < dim; ++component)
            vCornerVels[co][component]=u(component,co);
        vCornerPress[co]=u(_P_,co);
    }
	// Compute Upwind Shapes at Ip's and ConvectionLength here fot the Momentum Equation
    GetUpwindShapes(geo, vCornerVels, m_UpwindMethod, vIPVelUpwindShapesContiEq, vIPVelUpwindDependenciesContiEq, vConvLength);

	// todo: Implement for timedependent
	bool bTimeDependent = false;
	number dt = 0.0;

	// todo: switch
	// Compute Stabilized Velocities at IP's here (depending on Upwind Velocities)
	GetStabilizedShapes(geo, vCornerVels, vCornerPress, m_StabOptions, vIPVelUpwindShapesContiEq,vIPVelUpwindDependenciesContiEq, vConvLength,
									dt, bTimeDependent, vCornerVelsOld,
									m_Viscosity,
                                    vIPStabVelShapesContiEq);

	// loop Sub Control Volume Faces (SCVF)
	for(size_t i = 0; i < geo.num_scvf(); ++i)
	{
		// get current SCVF
		const typename TFVGeom<TElem, dim>::SCVF& scvf = geo.scvf(i);

		// loop integration point of SCVF
		for(size_t ip = 0; ip < scvf.num_ip(); ++ip)
		{
			// loop shape functions
			for(size_t co = 0; co < scvf.num_sh(); ++co)
			{
				////////////////////////////////////////////////////
				// Diffusiv Term (Momentum Equation)
				////////////////////////////////////////////////////

				// Compute flux at IP
				number flux = VecDot(scvf.global_grad(co, ip), scvf.normal());

				// Add flux term to local matrix
				for(size_t vel1 = 0; vel1 < dim; ++vel1)
				{
					//d(vel1, scvf.from(), vel1, co) -= m_Viscosity * flux;
					//d(vel1, scvf.to()  , vel1, co) += m_Viscosity * flux;

					for(size_t vel2 = 0; vel2 < dim; ++vel2)
					{
						//const number flux2 = m_Viscosity * scvf.global_grad(co, ip)[vel1] * scvf.normal()[vel2];
						//d(vel1, scvf.from(), vel2, co) -= flux2;
						//d(vel1, scvf.to()  , vel2, co) += flux2;
					}
				}

				////////////////////////////////////////////////////
				// Pressure Term (Momentum Equation)
				////////////////////////////////////////////////////

                for(size_t vel = 0; vel < dim; ++vel)
				{
					//d(vel, scvf.from(), _P_, co) += scvf.shape(co, ip) * scvf.normal()[vel];
					//d(vel, scvf.to()  , _P_, co) -= scvf.shape(co, ip) * scvf.normal()[vel];
				}

				////////////////////////////////////////////////////
				// Convective Term (Momentum Equation)
				////////////////////////////////////////////////////

				// Compute flux at IP
				flux = VecDot(scvf.global_grad(co, ip), scvf.normal());

                // Add flux term to local matrix
				for(size_t vel1 = 0; vel1 < dim; ++vel1)
				{
					//d(vel1, scvf.from(), vel1, co) += m_Viscosity * flux;
					//d(vel1, scvf.to()  , vel1, co) -= m_Viscosity * flux;

					for(size_t vel2 = 0; vel2 < dim; ++vel2)
					{
						//const number flux2 = m_Viscosity * scvf.global_grad(co, ip)[vel1] * scvf.normal()[vel2];
						//d(vel1, scvf.from(), vel2, co) += flux2;
						//d(vel1, scvf.to()  , vel2, co) -= flux2;
					}
				}

				// todo: implement convective term (momentum equation) use vIPVelStab

				// todo: implement mass equations use vIPVelStabMassEq
                //J(_P_, ,vel1, ) +=
                }
          }
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
FVNavierStokesElemDisc<TFVGeom, TDomain, TAlgebra>::
assemble_M(local_vector_type& d, const local_vector_type& u, number time)
{
	// get finite volume geometry
	TFVGeom<TElem, dim>& geo = FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();

	// loop Sub Control Volumes (SCV)
	for(size_t i = 0; i < geo.num_scv(); ++i)
	{
		// get current SCV
		const typename TFVGeom<TElem, dim>::SCV& scv = geo.scv(i);

		// get associated node
		const int co = scv.node_id();

		// loop velocity components
		for(size_t vel1 = 0; vel1 < dim; ++vel1)
		{
			// Add to local matrix
			d(vel1, co) += u(vel1, co) * scv.volume();
		}
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
FVNavierStokesElemDisc<TFVGeom, TDomain, TAlgebra>::
assemble_f(local_vector_type& d, number time)
{
	// we're done
	return true;
}

} // namespace ug


#endif /*__H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__NAVIER_STOKES__FV__NAVIER_STOKES_IMPL__*/
