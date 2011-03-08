/*
 * convection_diffusion_impl.h
 *
 *  Created on: 26.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__CONVECTION_DIFFUSION__FV1__CONVECTION_DIFFUSION_IMPL__
#define __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__CONVECTION_DIFFUSION__FV1__CONVECTION_DIFFUSION_IMPL__

#include "convection_diffusion.h"

namespace ug{


////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
//  Provide a generic implementation for all elements
//  (since this discretization can be implemented in a generic way)
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

template<typename TDomain, typename TAlgebra>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
inline
bool
FVConvectionDiffusionElemDisc<TDomain, TAlgebra>::
prepare_element_loop()
{
	// all this will be performed outside of the loop over the elements.
	// Therefore it is not time critical.

	typedef typename reference_element_traits<TElem>::reference_element_type
																ref_elem_type;
	static const int refDim = ref_elem_type::dim;

//	set local positions for rhs
	if(!TFVGeom<TElem, dim>::usesHangingNodes)
	{
		TFVGeom<TElem, dim>& geo = FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();
		m_Diff.template 	set_local_ips<refDim>(geo.scvf_local_ips(),
		                	                      geo.num_scvf_local_ips());
		m_ConvVel.template 	set_local_ips<refDim>(geo.scvf_local_ips(),
		                   	                      geo.num_scvf_local_ips());
		m_Rhs.template 		set_local_ips<refDim>(geo.scv_local_ips(),
		               		                      geo.num_scv_local_ips());
		m_Reaction.template set_local_ips<refDim>(geo.scv_local_ips(),
		                                          geo.num_scv_local_ips());
		m_MassScale.template set_local_ips<refDim>(geo.scv_local_ips(),
		                                           geo.num_scv_local_ips());
	}

	return true;
}

template<typename TDomain, typename TAlgebra>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
inline
bool
FVConvectionDiffusionElemDisc<TDomain, TAlgebra>::
finish_element_loop()
{
	// all this will be performed outside of the loop over the elements.
	// Therefore it is not time critical.

	return true;
}

template<typename TDomain, typename TAlgebra>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
inline
bool
FVConvectionDiffusionElemDisc<TDomain, TAlgebra>::
prepare_element(TElem* elem, const local_vector_type& u,
								const local_index_type& glob_ind)
{
	// this loop will be performed inside the loop over the elements.
	// Therefore, it is TIME CRITICAL
	typedef typename reference_element_traits<TElem>::reference_element_type
																ref_elem_type;
	static const int refDim = ref_elem_type::dim;

//	get corners
	m_vCornerCoords = this->template get_element_corners<TElem>(elem);

// 	Update Geometry for this element
	TFVGeom<TElem, dim>& geo = FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();
	if(!geo.update(elem, this->get_subset_handler(), &m_vCornerCoords[0]))
	{
		UG_LOG("FVConvectionDiffusionElemDisc::prepare_element:"
				" Cannot update Finite Volume Geometry.\n"); return false;
	}

//	set local positions for rhs
	if(TFVGeom<TElem, dim>::usesHangingNodes)
	{
		m_Diff.template 	set_local_ips<refDim>(geo.scvf_local_ips(),
												  geo.num_scvf_local_ips());
		m_ConvVel.template 	set_local_ips<refDim>(geo.scvf_local_ips(),
												  geo.num_scvf_local_ips());
		m_Rhs.template 		set_local_ips<refDim>(geo.scv_local_ips(),
												  geo.num_scv_local_ips());
		m_Reaction.template set_local_ips<refDim>(geo.scv_local_ips(),
												  geo.num_scv_local_ips());
		m_MassScale.template set_local_ips<refDim>(geo.scv_local_ips(),
												   geo.num_scv_local_ips());
	}

//	set global positions for rhs
	m_Diff.set_global_ips(geo.scvf_global_ips(), geo.num_scvf_global_ips());
	m_ConvVel.set_global_ips(geo.scvf_global_ips(), geo.num_scvf_global_ips());
	m_Rhs.set_global_ips(geo.scv_global_ips(), geo.num_scv_global_ips());
	m_Reaction.set_global_ips(geo.scv_global_ips(), geo.num_scv_global_ips());
	m_MassScale.set_global_ips(geo.scv_global_ips(), geo.num_scv_global_ips());

//	we're done
	return true;
}

template<typename TDomain, typename TAlgebra>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
inline
bool
FVConvectionDiffusionElemDisc<TDomain, TAlgebra>::
assemble_JA(local_matrix_type& J, const local_vector_type& u, number time)
{
// get finite volume geometry
	static TFVGeom<TElem, dim>& geo =
			FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();

//	Diff. Tensor times Gradient
	MathVector<dim> Dgrad;

//	Diffusion and Velocity Term
	size_t ip = 0;
	if(m_Diff.data_given() || m_ConvVel.data_given())
	{
	// 	loop Sub Control Volume Faces (SCVF)
		for(size_t i = 0; i < geo.num_scvf(); ++i)
		{
		// 	get current SCVF
			const typename TFVGeom<TElem, dim>::SCVF& scvf = geo.scvf(i);

		// 	loop integration point of SCVF
			for(size_t i = 0; i < scvf.num_ip(); ++i, ++ip)
			{
			// 	loop shape functions
				for(size_t j = 0; j < scvf.num_sh(); ++j)
				{
				////////////////////////////////////////////////////
				// Diffusive Term
				// (central discretization)
				////////////////////////////////////////////////////
					if(m_Diff.data_given())
					{
					// 	Compute Diffusion Tensor times Gradient
						MatVecMult(Dgrad, m_Diff[ip], scvf.global_grad(j, i));

					//	Compute flux at IP
						const number flux = VecDot(Dgrad, scvf.normal());

					// 	Add flux term to local matrix
						J(_C_, scvf.from(), _C_, j) -= flux;
						J(_C_, scvf.to()  , _C_, j) += flux;
					}

				////////////////////////////////////////////////////
				// Convective Term
				// (upwinding_amount == 1.0 -> full upwind;
				//  upwinding_amount == 0.0 -> full central disc)
				////////////////////////////////////////////////////

				// central part convection
					if(m_upwindAmount != 1.0)
					{
						if(m_ConvVel.data_given())
						{
						// 	Compute flux
							const number flux = (1.- m_upwindAmount) * scvf.shape(j, i)
												* VecDot(m_ConvVel[ip], scvf.normal());

						// 	Add flux to local matrix
							J(_C_, scvf.from(), _C_, j) += flux;
							J(_C_, scvf.to()  , _C_, j) -= flux;
						}
					}
				}

			// 	upwind part convection (does only depend on one shape function)
				if(m_upwindAmount != 0.0)
				{
					if(m_ConvVel.data_given())
					{
					// corner for upwind switch
						int up;

					//	Compute flux
						const number flux = m_upwindAmount * VecDot(m_ConvVel[ip], scvf.normal());

					// 	switch upwind direction
						if(flux >= 0.0) up = scvf.from(); else up = scvf.to();

					// 	Add flux to local matrix
						J(_C_, scvf.from(), _C_, up) += flux;
						J(_C_, scvf.to()  , _C_, up) -= flux;
					}
				}
			}
		}
	}

////////////////////////////////////////////////////
// Reaction Term
// (using lumping)
////////////////////////////////////////////////////

//	if no data for reaction, return
	if(!m_Reaction.data_given()) return true;

// 	loop Sub Control Volume (SCV)
	ip = 0;
	for(size_t i = 0; i < geo.num_scv(); ++i)
	{
	// 	get current SCV
		const typename TFVGeom<TElem, dim>::SCV& scv = geo.scv(i);

		number val = m_Reaction[ip++];

	// 	loop other ip positions
		for(size_t ip = 1; ip < scv.num_ip(); ++ip)
		{
			// TODO: Add here scaling factor for ip point
			val += m_Reaction[ip++];
		}

	// 	scale with volume of SCV
		val *= scv.volume();

	// 	get associated node
		const int co = scv.node_id();

	// 	Add to local matrix
		J(_C_, co, _C_, co) += val;
	}

// 	we're done
	return true;
}


template<typename TDomain, typename TAlgebra>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
inline
bool
FVConvectionDiffusionElemDisc<TDomain, TAlgebra>::
assemble_JM(local_matrix_type& J, const local_vector_type& u, number time)
{
// 	get finite volume geometry
	static TFVGeom<TElem, dim>& geo
		= FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();

// 	loop Sub Control Volumes (SCV)
	size_t ip = 0;
	for(size_t i = 0; i < geo.num_scv(); ++i)
	{
	// 	get current SCV
		const typename TFVGeom<TElem, dim>::SCV& scv = geo.scv(i);

	// 	get associated node
		const int co = scv.node_id();

	//	get value
		number val = scv.volume();

	//	multiply by scaling
		if(m_MassScale.data_given())
			val *= m_MassScale[ip++];

	// 	Add to local matrix
		J(_C_, co, _C_, co) += val;
	}

//	 we're done
	return true;
}


template<typename TDomain, typename TAlgebra>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
inline
bool
FVConvectionDiffusionElemDisc<TDomain, TAlgebra>::
assemble_A(local_vector_type& d, const local_vector_type& u, number time)
{
// 	get finite volume geometry
	static TFVGeom<TElem, dim>& geo
			= FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();

	MathVector<dim> grad_u;			// gradient of solution at ip
	number shape_u;					// solution at ip
	MathVector<dim> Dgrad_u;		// Diff.Tensor times gradient of solution

	size_t ip = 0;
	if(m_Diff.data_given() || m_ConvVel.data_given())
	{
	// 	loop Sub Control Volume Faces (SCVF)
		for(size_t i = 0; i < geo.num_scvf(); ++i)
		{
		// 	get current SCVF
			const typename TFVGeom<TElem, dim>::SCVF& scvf = geo.scvf(i);

		// 	loop integration points
			for(size_t i = 0; i < scvf.num_ip(); ++i, ++ip)
			{
			// 	reset values
				VecSet(grad_u, 0.0); shape_u = 0.0;

			// 	compute gradient and shape at ip
				for(size_t j = 0; j < scvf.num_sh(); ++j)
				{
					VecScaleAppend(grad_u, u(_C_,j), scvf.global_grad(j, i));
					shape_u += u(_C_,j) * scvf.shape(j, i);
				}

			/////////////////////////////////////////////////////
			// Diffusive Term
			// (central discretization)
			/////////////////////////////////////////////////////
				if(m_Diff.data_given())
				{
					MatVecMult(Dgrad_u, m_Diff[ip], grad_u);

				// 	Compute flux
					const number flux = VecDot(Dgrad_u, scvf.normal());

				// 	Add to local defect
					d(_C_, scvf.from()) -= flux;
					d(_C_, scvf.to()  ) += flux;
				}

			/////////////////////////////////////////////////////
			// Convective Term
			// (upwinding_amount == 1.0 -> full upwind;
			//  upwinding_amount == 0.0 -> full central disc)
			/////////////////////////////////////////////////////
				if(m_ConvVel.data_given())
				{
				// 	central part convection
					if(m_upwindAmount != 1.0)
					{
					// 	Compute flux at ip
						const number flux = (1.- m_upwindAmount) * shape_u
											* VecDot(m_ConvVel[ip], scvf.normal());

					// 	Add to local defect
						d(_C_, scvf.from()) += flux;
						d(_C_, scvf.to()  ) -= flux;
					}

				// 	upwind part convection
					if(m_upwindAmount != 0.0)
					{
					// 	compute flux at ip
						number flux = m_upwindAmount * VecDot(m_ConvVel[ip], scvf.normal());

					// 	Upwind switch
						if(flux >= 0.0) flux *= u(_C_, scvf.from());
						else flux *= u(_C_, scvf.to());

					// 	Add to local defect
						d(_C_, scvf.from()) += flux;
						d(_C_, scvf.to()  ) -= flux;
					}
				}
			}
		}
	}

//	if no reaction data given, return
	if(!m_Reaction.data_given()) return true;

// 	loop Sub Control Volumes (SCV)
	ip = 0;
	for(size_t i = 0; i < geo.num_scv(); ++i)
	{
	// 	get current SCV
		const typename TFVGeom<TElem, dim>::SCV& scv = geo.scv(i);

	// 	first ip
		number val = m_Reaction[ip++];

	// 	loop other ip
		for(size_t ip = 1; ip < scv.num_ip(); ++ip)
		{
			// TODO: Add weight for integration
			val += m_Reaction[ip++];
		}

	// 	scale with volume of SCV
		val *= scv.volume();

	// 	get associated node
		const int co = scv.node_id();

	// 	Add to local defect
		d(_C_, co) += u(_C_, co) * val;
	}

// 	we're done
	return true;
}


template<typename TDomain, typename TAlgebra>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
inline
bool
FVConvectionDiffusionElemDisc<TDomain, TAlgebra>::
assemble_M(local_vector_type& d, const local_vector_type& u, number time)
{
// 	get finite volume geometry
	static TFVGeom<TElem, dim>& geo
		= FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();

// 	loop Sub Control Volumes (SCV)
	size_t ip = 0;
	for(size_t i = 0; i < geo.num_scv(); ++i)
	{
	// 	get current SCV
		const typename TFVGeom<TElem, dim>::SCV& scv = geo.scv(i);

	// 	get associated node
		const int co = scv.node_id();

	//	mass value
		number val = u(_C_, co) * scv.volume();

	//	multiply by scaling
		if(m_MassScale.data_given())
			val *= m_MassScale[ip++];

	// 	Add to local defect
		d(_C_, co) += val;
	}

// 	we're done
	return true;
}


template<typename TDomain, typename TAlgebra>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
inline
bool
FVConvectionDiffusionElemDisc<TDomain, TAlgebra>::
assemble_f(local_vector_type& d, number time)
{
//	if zero data given, return
	if(!m_Rhs.data_given()) return true;

// 	get finite volume geometry
	static TFVGeom<TElem, dim>& geo
		= FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();

// 	loop Sub Control Volumes (SCV)
	size_t ip = 0;
	for(size_t i = 0; i < geo.num_scv(); ++i)
	{
	// 	get current SCV
		const typename TFVGeom<TElem, dim>::SCV& scv = geo.scv(i);

	// 	first value
		number val = m_Rhs[ip++];

	// 	other values
		for(size_t i = 1; i < scv.num_ip(); ++i)
		{
			number ip_val;
			ip_val = m_Rhs[ip++];

			// TODO: add weights for integration
			val += ip_val;
		}

	// 	scale with volume of SCV
		val *= scv.volume();

	// 	get associated node
		const int co = scv.node_id();

	// 	Add to local rhs
		d(_C_, co) += val;
	}

// 	we're done
	return true;
}


//	computes the linearized defect w.r.t to the velocity
template<typename TDomain, typename TAlgebra>
template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
bool
FVConvectionDiffusionElemDisc<TDomain, TAlgebra>::
lin_defect_velocity(const local_vector_type& u)
{
//  get finite volume geometry
	static TFVGeom<TElem, dim>& geo = FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();

//  loop Sub Control Volume Faces (SCVF)
	size_t ip = 0;
	for(size_t i = 0; i < geo.num_scvf(); ++i)
	{
	// get current SCVF
		const typename TFVGeom<TElem, dim>::SCVF& scvf = geo.scvf(i);

	// loop integration point of SCVF
		for(size_t i = 0; i < scvf.num_ip(); ++i, ++ip)
		{
		//  reset
			for(size_t j = 0; j < scvf.num_sh(); ++j)
			{
				m_ConvVel.lin_defect(ip, _C_, j) = 0.0;
			}

		// central part convection
			number scale = (1.- m_upwindAmount);

			number shape_u = 0.0;
			for(size_t j = 0; j < scvf.num_sh(); ++j)
				shape_u += u(_C_,j) * scvf.shape(j, i);

			scale *= shape_u;

		// upwind part convection
			const number flux = m_upwindAmount
								* VecDot(m_ConvVel[ip], scvf.normal());

			if(flux >= 0.0) scale += m_upwindAmount* u(_C_, scvf.from());
			else scale += m_upwindAmount * u(_C_, scvf.to());

			MathVector<dim> linDefect = scvf.normal();
			linDefect *= scale;

			m_ConvVel.lin_defect(ip, _C_, scvf.from()) += linDefect;
			m_ConvVel.lin_defect(ip, _C_, scvf.to()) -= linDefect;
		}
	}

//	we're done
	return true;
}

//	computes the linearized defect w.r.t to the velocity
template<typename TDomain, typename TAlgebra>
template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
bool
FVConvectionDiffusionElemDisc<TDomain, TAlgebra>::
lin_defect_diffusion(const local_vector_type& u)
{
//  get finite volume geometry
	static TFVGeom<TElem, dim>& geo = FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();

//  loop Sub Control Volume Faces (SCVF)
	size_t ip = 0;
	for(size_t i = 0; i < geo.num_scvf(); ++i)
	{
	// get current SCVF
		const typename TFVGeom<TElem, dim>::SCVF& scvf = geo.scvf(i);

	// loop integration point of SCVF
		for(size_t i = 0; i < scvf.num_ip(); ++i, ++ip)
		{
			// 	reset values
			for(size_t j = 0; j < scvf.num_sh(); ++j)
			{
				m_Diff.lin_defect(ip, _C_, j) = 0.0;
			}

			MathVector<dim> grad_u;
		// 	compute gradient at ip
			VecSet(grad_u, 0.0);
			for(size_t j = 0; j < scvf.num_sh(); ++j)
			{
				VecScaleAppend(grad_u, u(_C_,j), scvf.global_grad(j, i));
			}

			MathMatrix<dim,dim> linDefect;

			for(size_t k=0; k < (size_t)dim; ++k)
				for(size_t j = 0; j < (size_t)dim; ++j)
					linDefect(j,k) = (scvf.normal())[j] * grad_u[k];

			m_Diff.lin_defect(ip, _C_, scvf.from()) -= linDefect;
			m_Diff.lin_defect(ip, _C_, scvf.to()) += linDefect;
		}
	}

//	we're done
	return true;
}


//	computes the linearized defect w.r.t to the velocity
template<typename TDomain, typename TAlgebra>
template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
bool
FVConvectionDiffusionElemDisc<TDomain, TAlgebra>::
compute_concentration_export(const local_vector_type& u, bool compDeriv)
{
//  get finite volume geometry
	static TFVGeom<TElem, dim>& geo
		= FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();

	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
	static const size_t refDim = ref_elem_type::dim;
	static const size_t num_co = ref_elem_type::num_corners;

//	loop all ip series
	for(size_t s = 0; s < m_exConcentration.num_series(); ++s)
	{
	//	FV1 SCVF ip
		if(m_exConcentration.template local_ips<refDim>(s)
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
				//	compute concentration at ip
					number& cIP = m_exConcentration.value(s, ip);
					for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
						cIP += u(_C_, sh) * scvf.shape(sh, j);

				//	compute derivative w.r.t. to unknowns iff needed
					if(compDeriv)
					{
					//	get field of derivatives
						number* cIP_c = m_exConcentration.deriv(s, ip, _C_);

					//	set shapes
						for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
							cIP_c[sh] = scvf.shape(sh, j);
					}
				}
			}

		//	this position type is done
			continue;
		}

	//	FV1 SCV ip
		if(m_exConcentration.template local_ips<refDim>(s)
				== geo.scv_local_ips())
		{
		//	Loop Corners
			for(size_t ip = 0; ip < num_co; ip++)
			{
			//	Compute Gradients and concentration at ip
				m_exConcentration.value(s, ip) = u(_C_, ip);

			//	set derivatives if needed
				if(compDeriv)
				{
					number* cIP_c = m_exConcentration.deriv(s, ip, _C_);

					for(size_t sh = 0; sh < num_co; ++sh)
						cIP_c[sh] = (ip==sh) ? 1.0 : 0.0;
				}
			}

		//	this position type is done
			continue;
		}

	// 	others not implemented
		UG_LOG("Evaluation not implemented.");
		return false;
	}

//	we're done
	return true;
}

//	computes the linearized defect w.r.t to the velocity
template<typename TDomain, typename TAlgebra>
template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
bool
FVConvectionDiffusionElemDisc<TDomain, TAlgebra>::
compute_concentration_grad_export(const local_vector_type& u, bool compDeriv)
{
// 	Get finite volume geometry
	static const TFVGeom<TElem, dim>& geo =
				FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();

	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
	static const size_t refDim = ref_elem_type::dim;

	for(size_t s = 0; s < m_exConcentrationGrad.num_series(); ++s)
	{
	//	FV1 SCVF ip
		if(m_exConcentrationGrad.template local_ips<refDim>(s)
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
					MathVector<dim>& cIP = m_exConcentrationGrad.value(s, ip);

					VecSet(cIP, 0.0);
					for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
						VecScaleAppend(cIP, u(_C_, sh), scvf.global_grad(sh, j));

					if(compDeriv)
					{
						MathVector<dim>* cIP_c = m_exConcentrationGrad.deriv(s, ip, _C_);

						for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
							cIP_c[sh] = scvf.global_grad(sh, j);
					}
				}
			}

		//	this position type is done
			continue;
		}

		// others not implemented
		UG_LOG("Evaluation not implemented.");
		return false;
	}

//	we're done
	return true;
}


} // namespace ug


#endif /*__H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__CONVECTION_DIFFUSION__FV1__CONVECTION_DIFFUSION_IMPL__*/
