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
		m_imDiffusion.template 	set_local_ips<refDim>(geo.scvf_local_ips(),
		                	                      geo.num_scvf_local_ips());
		m_imVelocity.template 	set_local_ips<refDim>(geo.scvf_local_ips(),
		                   	                      geo.num_scvf_local_ips());
		m_imSource.template 		set_local_ips<refDim>(geo.scv_local_ips(),
		               		                      geo.num_scv_local_ips());
		m_imReaction.template set_local_ips<refDim>(geo.scv_local_ips(),
		                                          geo.num_scv_local_ips());
		m_imMassScale.template set_local_ips<refDim>(geo.scv_local_ips(),
		                                           geo.num_scv_local_ips());
	}

//	check, that upwind has been set
	if(m_pConvShape == NULL)
	{
		UG_LOG("ERROR in 'FVConvectionDiffusionElemDisc::prepare_element_loop':"
				" Upwind has not been set.\n");
		return false;
	}

//	init upwind for element type
	if(!m_pConvShape->template set_geometry_type<TFVGeom<TElem, dim> >())
	{
		UG_LOG("ERROR in 'FVConvectionDiffusionElemDisc::prepare_element_loop':"
				" Cannot init upwind for element type.\n");
		return false;
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
bool
FVConvectionDiffusionElemDisc<TDomain, TAlgebra>::
time_point_changed(number time)
{
//	set new time point at imports
	m_imDiffusion.set_time(time);
	m_imVelocity.set_time(time);
	m_imSource.set_time(time);
	m_imReaction.set_time(time);
	m_imMassScale.set_time(time);

//	this disc does not need the old time solutions, thus, return false
	return false;
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
		m_imDiffusion.template 	set_local_ips<refDim>(geo.scvf_local_ips(),
												  geo.num_scvf_local_ips());
		m_imVelocity.template 	set_local_ips<refDim>(geo.scvf_local_ips(),
												  geo.num_scvf_local_ips());
		m_imSource.template 		set_local_ips<refDim>(geo.scv_local_ips(),
												  geo.num_scv_local_ips());
		m_imReaction.template set_local_ips<refDim>(geo.scv_local_ips(),
												  geo.num_scv_local_ips());
		m_imMassScale.template set_local_ips<refDim>(geo.scv_local_ips(),
												   geo.num_scv_local_ips());
	}

//	set global positions for rhs
	m_imDiffusion.set_global_ips(geo.scvf_global_ips(), geo.num_scvf_global_ips());
	m_imVelocity.set_global_ips(geo.scvf_global_ips(), geo.num_scvf_global_ips());
	m_imSource.set_global_ips(geo.scv_global_ips(), geo.num_scv_global_ips());
	m_imReaction.set_global_ips(geo.scv_global_ips(), geo.num_scv_global_ips());
	m_imMassScale.set_global_ips(geo.scv_global_ips(), geo.num_scv_global_ips());

//	we're done
	return true;
}

template<typename TDomain, typename TAlgebra>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
inline
bool
FVConvectionDiffusionElemDisc<TDomain, TAlgebra>::
assemble_JA(local_matrix_type& J, const local_vector_type& u)
{
// get finite volume geometry
	static TFVGeom<TElem, dim>& geo =
			FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();

//	Diff. Tensor times Gradient
	MathVector<dim> Dgrad;

//	get conv shapes
	const IConvectionShapes<dim>& convShape = get_updated_conv_shapes(geo);

//	Diffusion and Velocity Term
	if(m_imDiffusion.data_given() || m_imVelocity.data_given())
	{
	// 	loop Sub Control Volume Faces (SCVF)
		for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
		{
		// 	get current SCVF
			const typename TFVGeom<TElem, dim>::SCVF& scvf = geo.scvf(ip);

		////////////////////////////////////////////////////
		// Diffusive Term
		////////////////////////////////////////////////////
			if(m_imDiffusion.data_given())
			{
			// 	loop shape functions
				for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
				{
				// 	Compute Diffusion Tensor times Gradient
					MatVecMult(Dgrad, m_imDiffusion[ip], scvf.global_grad(sh));

				//	Compute flux at IP
					const number D_diff_flux = VecDot(Dgrad, scvf.normal());

				// 	Add flux term to local matrix
					J(_C_, scvf.from(), _C_, sh) -= D_diff_flux;
					J(_C_, scvf.to()  , _C_, sh) += D_diff_flux;
				}
			}

		////////////////////////////////////////////////////
		// Convective Term
		////////////////////////////////////////////////////
			if(m_imVelocity.data_given())
			{
			//	Add Flux contribution
				for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
				{
					const number D_conv_flux = convShape.conv_shape(ip, sh);

				//	Add fkux term to local matrix
					J(_C_, scvf.from(), _C_, sh) += D_conv_flux;
					J(_C_, scvf.to(),   _C_, sh) -= D_conv_flux;
				}
			}
		}
	}

////////////////////////////////////////////////////
// Reaction Term (using lumping)
////////////////////////////////////////////////////

//	if no data for reaction, return
	if(!m_imReaction.data_given()) return true;

// 	loop Sub Control Volume (SCV)
	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
	// 	get current SCV
		const typename TFVGeom<TElem, dim>::SCV& scv = geo.scv(ip);

	// 	get associated node
		const int co = scv.node_id();

	// 	Add to local matrix
		J(_C_, co, _C_, co) += m_imReaction[ip] * scv.volume();
	}

// 	we're done
	return true;
}


template<typename TDomain, typename TAlgebra>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
inline
bool
FVConvectionDiffusionElemDisc<TDomain, TAlgebra>::
assemble_JM(local_matrix_type& J, const local_vector_type& u)
{
// 	get finite volume geometry
	static TFVGeom<TElem, dim>& geo
		= FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();

// 	loop Sub Control Volumes (SCV)
	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
	// 	get current SCV
		const typename TFVGeom<TElem, dim>::SCV& scv = geo.scv(ip);

	// 	get associated node
		const int co = scv.node_id();

	//	get value
		number val = scv.volume();

	//	multiply by scaling
		if(m_imMassScale.data_given())
			val *= m_imMassScale[ip];

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
assemble_A(local_vector_type& d, const local_vector_type& u)
{
// 	get finite volume geometry
	static TFVGeom<TElem, dim>& geo
			= FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();

//	get conv shapes
	const IConvectionShapes<dim>& convShape = get_updated_conv_shapes(geo);

	if(m_imDiffusion.data_given() || m_imVelocity.data_given())
	{
	// 	loop Sub Control Volume Faces (SCVF)
		for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
		{
		// 	get current SCVF
			const typename TFVGeom<TElem, dim>::SCVF& scvf = geo.scvf(ip);

		/////////////////////////////////////////////////////
		// Diffusive Term
		/////////////////////////////////////////////////////
			if(m_imDiffusion.data_given())
			{
			//	to compute D \nabla c
				MathVector<dim> Dgrad_c, grad_c;

			// 	compute gradient and shape at ip
				VecSet(grad_c, 0.0);
				for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
					VecScaleAppend(grad_c, u(_C_,sh), scvf.global_grad(sh));

			//	scale by diffusion tensor
				MatVecMult(Dgrad_c, m_imDiffusion[ip], grad_c);

			// 	Compute flux
				const number diff_flux = VecDot(Dgrad_c, scvf.normal());

			// 	Add to local defect
				d(_C_, scvf.from()) -= diff_flux;
				d(_C_, scvf.to()  ) += diff_flux;
			}

		/////////////////////////////////////////////////////
		// Convective Term
		/////////////////////////////////////////////////////
			if(m_imVelocity.data_given())
			{
			//	sum up convective flux using convection shapes
				number conv_flux = 0.0;
				for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
					conv_flux += u(_C_, sh) * convShape.conv_shape(ip, sh);

			//  add to local defect
				d(_C_, scvf.from()) += conv_flux;
				d(_C_, scvf.to()  ) -= conv_flux;
			}
		}
	}

//	if no reaction data given, return
	if(!m_imReaction.data_given()) return true;

// 	loop Sub Control Volumes (SCV)
	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
	// 	get current SCV
		const typename TFVGeom<TElem, dim>::SCV& scv = geo.scv(ip);

	// 	get associated node
		const int co = scv.node_id();

	// 	Add to local defect
		d(_C_, co) += u(_C_, co) * m_imReaction[ip] * scv.volume();
	}

// 	we're done
	return true;
}


template<typename TDomain, typename TAlgebra>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
inline
bool
FVConvectionDiffusionElemDisc<TDomain, TAlgebra>::
assemble_M(local_vector_type& d, const local_vector_type& u)
{
// 	get finite volume geometry
	static TFVGeom<TElem, dim>& geo
		= FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();

// 	loop Sub Control Volumes (SCV)
	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
	// 	get current SCV
		const typename TFVGeom<TElem, dim>::SCV& scv = geo.scv(ip);

	// 	get associated node
		const int co = scv.node_id();

	//	mass value
		number val = u(_C_, co) * scv.volume();

	//	multiply by scaling
		if(m_imMassScale.data_given())
			val *= m_imMassScale[ip];

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
assemble_f(local_vector_type& d)
{
//	if zero data given, return
	if(!m_imSource.data_given()) return true;

// 	get finite volume geometry
	static TFVGeom<TElem, dim>& geo
		= FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();

// 	loop Sub Control Volumes (SCV)
	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
	// 	get current SCV
		const typename TFVGeom<TElem, dim>::SCV& scv = geo.scv(ip);

	// 	get associated node
		const int co = scv.node_id();

	// 	Add to local rhs
		d(_C_, co) += m_imSource[ip] * scv.volume();
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

//	get conv shapes
	const IConvectionShapes<dim>& convShape = get_updated_conv_shapes(geo);

//	reset the values for the linearized defect
	m_imVelocity.clear_lin_defect();

//  loop Sub Control Volume Faces (SCVF)
	for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
	{
	// get current SCVF
		const typename TFVGeom<TElem, dim>::SCVF& scvf = geo.scvf(ip);

	//	get derivatives to shapes at ip
		MathVector<dim>* vLinDef = m_imVelocity.lin_defect(ip, _C_);

	//	sum up contributions of convection shapes
		MathVector<dim> linDefect;
		VecSet(linDefect, 0.0);
		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
			VecScaleAppend(linDefect, u(_C_,sh), convShape.conv_shape_vel(ip, sh));

	//	add parts for both sides of scvf
		vLinDef[scvf.from()] += linDefect;
		vLinDef[scvf.to()] -= linDefect;
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

//	get conv shapes
	const IConvectionShapes<dim>& convShape = get_updated_conv_shapes(geo);

//	cleat the linearized defect
	m_imDiffusion.clear_lin_defect();

//  loop Sub Control Volume Faces (SCVF)
	for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
	{
	// get current SCVF
		const typename TFVGeom<TElem, dim>::SCVF& scvf = geo.scvf(ip);

	// 	compute gradient at ip
		MathVector<dim> grad_u;	VecSet(grad_u, 0.0);
		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
			VecScaleAppend(grad_u, u(_C_,sh), scvf.global_grad(sh));

	//	compute the lin defect at this ip
		MathMatrix<dim,dim> linDefect;

	//	part coming from -\nabla u * \vec{n}
		for(size_t k=0; k < (size_t)dim; ++k)
			for(size_t j = 0; j < (size_t)dim; ++j)
				linDefect(j,k) = (scvf.normal())[j] * grad_u[k];

	//	add contribution from convection shapes
		if(convShape.non_zero_deriv_diffusion())
			for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
				MatAdd(linDefect, convShape.conv_shape_diffusion(ip, sh), u(_C_, sh));

	//	add contributions
		m_imDiffusion.lin_defect(ip, _C_, scvf.from()) -= linDefect;
		m_imDiffusion.lin_defect(ip, _C_, scvf.to()) += linDefect;
	}

//	we're done
	return true;
}

//	computes the linearized defect w.r.t to the reaction
template<typename TDomain, typename TAlgebra>
template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
bool
FVConvectionDiffusionElemDisc<TDomain, TAlgebra>::
lin_defect_reaction(const local_vector_type& u)
{
//  get finite volume geometry
	static TFVGeom<TElem, dim>& geo = FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();

// 	loop Sub Control Volumes (SCV)
	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
	// 	get current SCV
		const typename TFVGeom<TElem, dim>::SCV& scv = geo.scv(ip);

	// 	get associated node
		const int co = scv.node_id();

	// 	Add to local defect
		m_imReaction.lin_defect(ip, _C_, co) = u(_C_, co) * scv.volume();
	}

//	we're done
	return true;
}

//	computes the linearized defect w.r.t to the source
template<typename TDomain, typename TAlgebra>
template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
bool
FVConvectionDiffusionElemDisc<TDomain, TAlgebra>::
lin_defect_source(const local_vector_type& u)
{
//  get finite volume geometry
	static TFVGeom<TElem, dim>& geo = FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();

// 	loop Sub Control Volumes (SCV)
	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
	// 	get current SCV
		const typename TFVGeom<TElem, dim>::SCV& scv = geo.scv(ip);

	// 	get associated node
		const int co = scv.node_id();

	// 	Add to local defect
		m_imSource.lin_defect(ip, _C_, co) = scv.volume();
	}

//	we're done
	return true;
}

//	computes the linearized defect w.r.t to the mass scale
template<typename TDomain, typename TAlgebra>
template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
bool
FVConvectionDiffusionElemDisc<TDomain, TAlgebra>::
lin_defect_mass_scale(const local_vector_type& u)
{
//  get finite volume geometry
	static TFVGeom<TElem, dim>& geo = FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();

// 	loop Sub Control Volumes (SCV)
	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
	// 	get current SCV
		const typename TFVGeom<TElem, dim>::SCV& scv = geo.scv(ip);

	// 	get associated node
		const int co = scv.node_id();

	// 	Add to local defect
		m_imMassScale.lin_defect(ip, _C_, co) = u(_C_, co) * scv.volume();
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

//	reference element
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
	static const size_t refDim = ref_elem_type::dim;
	static const size_t numSh = ref_elem_type::num_corners;

//	loop all ip series
	for(size_t s = 0; s < m_exConcentration.num_series(); ++s)
	{
	//	get position array
		const MathVector<refDim>* vLocalIP
							= m_exConcentration.template local_ips<refDim>(s);

	//	FV1 SCVF ip
		if(vLocalIP	== geo.scvf_local_ips())
		{
		//	Loop Sub Control Volume Faces (SCVF)
			for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
			{
			// 	Get current SCVF
				const typename TFVGeom<TElem, dim>::SCVF& scvf = geo.scvf(ip);

			//	compute concentration at ip
				number& cIP = m_exConcentration.value(s, ip);
				cIP = 0.0;
				for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
					cIP += u(_C_, sh) * scvf.shape(sh);

			//	compute derivative w.r.t. to unknowns iff needed
				if(compDeriv)
				{
				//	get field of derivatives
					number* cIP_c = m_exConcentration.deriv(s, ip, _C_);

				//	set shapes
					for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
						cIP_c[sh] = scvf.shape(sh);
				}
			}
		}
	//	FV1 SCV ip
		else if(vLocalIP == geo.scv_local_ips())
		{
		//	Loop Corners
			for(size_t sh = 0; sh < numSh; ++sh)
			{
			//	Compute Gradients and concentration at ip
				m_exConcentration.value(s, sh) = u(_C_, sh);

			//	set derivatives if needed
				if(compDeriv)
				{
					number* cIP_c = m_exConcentration.deriv(s, sh, _C_);

					for(size_t sh2 = 0; sh2 < numSh; ++sh2)
						cIP_c[sh2] = (sh==sh2) ? 1.0 : 0.0;
				}
			}
		}
	// 	others not implemented
		else
		{
			UG_LOG("Evaluation not implemented.");
			return false;
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
compute_concentration_grad_export(const local_vector_type& u, bool compDeriv)
{
// 	Get finite volume geometry
	static const TFVGeom<TElem, dim>& geo =
				FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();

//	get reference element dimension
	static const size_t refDim = TFVGeom<TElem, dim>::dim;

//	loop all requested ip series
	for(size_t s = 0; s < m_exConcentrationGrad.num_series(); ++s)
	{
	//	get position array
		const MathVector<refDim>* vLocalIP
						= m_exConcentrationGrad.template local_ips<refDim>(s);

	//	FV1 SCVF ip
		if(vLocalIP == geo.scvf_local_ips())
		{
		//	Loop Sub Control Volume Faces (SCVF)
			for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
			{
			// 	Get current SCVF
				const typename TFVGeom<TElem, dim>::SCVF& scvf = geo.scvf(ip);

			//	Compute Gradients and concentration at ip
				MathVector<dim>& cIP = m_exConcentrationGrad.value(s, ip);

				VecSet(cIP, 0.0);
				for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
					VecScaleAppend(cIP, u(_C_, sh), scvf.global_grad(sh));

				if(compDeriv)
				{
					MathVector<dim>* cIP_c = m_exConcentrationGrad.deriv(s, ip, _C_);

					for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
						cIP_c[sh] = scvf.global_grad(sh);
				}
			}
		}
	// others not implemented
		else
		{
			UG_LOG("Evaluation not implemented.");
			return false;
		}
	}

//	we're done
	return true;
};

//	computes the linearized defect w.r.t to the velocity
template<typename TDomain, typename TAlgebra>
const typename FVConvectionDiffusionElemDisc<TDomain, TAlgebra>::conv_shape_type&
FVConvectionDiffusionElemDisc<TDomain, TAlgebra>::
get_updated_conv_shapes(const FVGeometryBase& geo)
{
//	compute upwind shapes for transport equation
//	\todo: we should move this computation into the preparation part of the
//			disc, to only compute the shapes once, reusing them several times.
	if(m_imVelocity.data_given())
	{
	//	get diffusion at ips
		const MathMatrix<dim, dim>* vDiffusion = NULL;
		if(m_imDiffusion.data_given()) vDiffusion = m_imDiffusion.values();

	//	update convection shapes
		if(!m_pConvShape->update(&geo, m_imVelocity.values(), vDiffusion, true))
		{
			UG_LOG("ERROR in 'ConvectionDiffusionElemDisc::assemble_JA': "
					"Cannot compute convection shapes.\n");
		}
	}

//	return a const (!!) reference to the upwind
	return *const_cast<const IConvectionShapes<dim>*>(m_pConvShape);
}

} // namespace ug


#endif /*__H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__CONVECTION_DIFFUSION__FV1__CONVECTION_DIFFUSION_IMPL__*/
