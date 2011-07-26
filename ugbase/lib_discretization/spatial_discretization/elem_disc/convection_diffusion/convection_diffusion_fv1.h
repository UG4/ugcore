/*
 * convection_diffusion_fv1.cpp
 *
 *  Created on: 26.02.2010
 *      Author: andreasvogel
 */

#include "convection_diffusion.h"

#include "common/util/provider.h"
#include "lib_discretization/spatial_discretization/disc_util/finite_volume_geometry.h"
#include "lib_discretization/spatial_discretization/disc_util/hanging_finite_volume_geometry.h"

namespace ug{


////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
//  Provide a generic implementation for all elements
//  (since this discretization can be implemented in a generic way)
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
template<typename TElem, typename TFVGeom>
bool ConvectionDiffusionElemDisc<TDomain>::
elem_loop_prepare_fv1()
{
//	reference dimension
	static const int refDim = reference_element_traits<TElem>::dim;

//	set local positions
	if(!TFVGeom::usesHangingNodes)
	{
		TFVGeom& geo = Provider::get<TFVGeom>();
		m_imDiffusion.template 	set_local_ips<refDim>(geo.scvf_local_ips(),
		                       	                      geo.num_scvf_ips());
		m_imVelocity.template 	set_local_ips<refDim>(geo.scvf_local_ips(),
		                      	                      geo.num_scvf_ips());
		m_imSource.template 	set_local_ips<refDim>(geo.scv_local_ips(),
		                    	                      geo.num_scv_ips());
		m_imReaction.template 	set_local_ips<refDim>(geo.scv_local_ips(),
		                      	                      geo.num_scv_ips());
		m_imMassScale.template 	set_local_ips<refDim>(geo.scv_local_ips(),
		                       	                      geo.num_scv_ips());
	}

//	check, that upwind has been set
	if(m_pConvShape == NULL)
	{
		UG_LOG("ERROR in 'ConvectionDiffusionElemDisc::prepare_element_loop':"
				" Upwind has not been set.\n");
		return false;
	}

//	init upwind for element type
	if(!m_pConvShape->template set_geometry_type<TFVGeom>())
	{
		UG_LOG("ERROR in 'ConvectionDiffusionElemDisc::prepare_element_loop':"
				" Cannot init upwind for element type.\n");
		return false;
	}

//	done
	return true;
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
bool ConvectionDiffusionElemDisc<TDomain>::
elem_loop_finish_fv1()
{
//	nothing to do
	return true;
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
bool ConvectionDiffusionElemDisc<TDomain>::
elem_prepare_fv1(TElem* elem, const local_vector_type& u)
{
//	get reference elements
	static const int refDim = reference_element_traits<TElem>::dim;

//	get corners
	m_vCornerCoords = this->template get_element_corners<TElem>(elem);

// 	Update Geometry for this element
	static TFVGeom& geo = Provider::get<TFVGeom>();

	if(!geo.update(elem, this->get_subset_handler(), &m_vCornerCoords[0]))
	{
		UG_LOG("ConvectionDiffusionElemDisc::prepare_element:"
				" Cannot update Finite Volume Geometry.\n"); return false;
	}

//	set local positions
	if(TFVGeom::usesHangingNodes)
	{
		m_imDiffusion.template 	set_local_ips<refDim>(geo.scvf_local_ips(),
		                       	                      geo.num_scvf_ips());
		m_imVelocity.template 	set_local_ips<refDim>(geo.scvf_local_ips(),
		                      	                      geo.num_scvf_ips());
		m_imSource.template 	set_local_ips<refDim>(geo.scv_local_ips(),
		                    	                      geo.num_scv_ips());
		m_imReaction.template 	set_local_ips<refDim>(geo.scv_local_ips(),
		                      	                      geo.num_scv_ips());
		m_imMassScale.template 	set_local_ips<refDim>(geo.scv_local_ips(),
		                       	                      geo.num_scv_ips());
	}

//	set global positions
	m_imDiffusion.	set_global_ips(geo.scvf_global_ips(), geo.num_scvf_ips());
	m_imVelocity.	set_global_ips(geo.scvf_global_ips(), geo.num_scvf_ips());
	m_imSource.		set_global_ips(geo.scv_global_ips(), geo.num_scv_ips());
	m_imReaction.	set_global_ips(geo.scv_global_ips(), geo.num_scv_ips());
	m_imMassScale.	set_global_ips(geo.scv_global_ips(), geo.num_scv_ips());

//	we're done
	return true;
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
bool ConvectionDiffusionElemDisc<TDomain>::
elem_JA_fv1(local_matrix_type& J, const local_vector_type& u)
{
// get finite volume geometry
	const static TFVGeom& geo = Provider::get<TFVGeom>();

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
			const typename TFVGeom::SCVF& scvf = geo.scvf(ip);

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
					const number D_conv_flux = convShape(ip, sh);

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
		const typename TFVGeom::SCV& scv = geo.scv(ip);

	// 	get associated node
		const int co = scv.node_id();

	// 	Add to local matrix
		J(_C_, co, _C_, co) += m_imReaction[ip] * scv.volume();
	}

// 	we're done
	return true;
}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
bool ConvectionDiffusionElemDisc<TDomain>::
elem_JM_fv1(local_matrix_type& J, const local_vector_type& u)
{
// 	get finite volume geometry
	const static TFVGeom& geo = Provider::get<TFVGeom>();

// 	loop Sub Control Volumes (SCV)
	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
	// 	get current SCV
		const typename TFVGeom::SCV& scv = geo.scv(ip);

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


template<typename TDomain>
template<typename TElem, typename TFVGeom>
bool ConvectionDiffusionElemDisc<TDomain>::
elem_dA_fv1(local_vector_type& d, const local_vector_type& u)
{
// 	get finite volume geometry
	const static TFVGeom& geo = Provider::get<TFVGeom>();

//	get conv shapes
	const IConvectionShapes<dim>& convShape = get_updated_conv_shapes(geo);

	if(m_imDiffusion.data_given() || m_imVelocity.data_given())
	{
	// 	loop Sub Control Volume Faces (SCVF)
		for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
		{
		// 	get current SCVF
			const typename TFVGeom::SCVF& scvf = geo.scvf(ip);

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
					conv_flux += u(_C_, sh) * convShape(ip, sh);

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
		const typename TFVGeom::SCV& scv = geo.scv(ip);

	// 	get associated node
		const int co = scv.node_id();

	// 	Add to local defect
		d(_C_, co) += u(_C_, co) * m_imReaction[ip] * scv.volume();
	}

// 	we're done
	return true;
}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
bool ConvectionDiffusionElemDisc<TDomain>::
elem_dM_fv1(local_vector_type& d, const local_vector_type& u)
{
// 	get finite volume geometry
	const static TFVGeom& geo = Provider::get<TFVGeom>();

// 	loop Sub Control Volumes (SCV)
	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
	// 	get current SCV
		const typename TFVGeom::SCV& scv = geo.scv(ip);

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


template<typename TDomain>
template<typename TElem, typename TFVGeom>
bool ConvectionDiffusionElemDisc<TDomain>::
elem_rhs_fv1(local_vector_type& d)
{
//	if zero data given, return
	if(!m_imSource.data_given()) return true;

// 	get finite volume geometry
	const static TFVGeom& geo = Provider::get<TFVGeom>();

// 	loop Sub Control Volumes (SCV)
	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
	// 	get current SCV
		const typename TFVGeom::SCV& scv = geo.scv(ip);

	// 	get associated node
		const int co = scv.node_id();

	// 	Add to local rhs
		d(_C_, co) += m_imSource[ip] * scv.volume();
	}

// 	we're done
	return true;
}


//	computes the linearized defect w.r.t to the velocity
template<typename TDomain>
template <typename TElem, typename TFVGeom>
bool ConvectionDiffusionElemDisc<TDomain>::
lin_defect_velocity_fv1(const local_vector_type& u)
{
// 	get finite volume geometry
	const static TFVGeom& geo = Provider::get<TFVGeom>();

//	get conv shapes
	const IConvectionShapes<dim>& convShape = get_updated_conv_shapes(geo);

//	reset the values for the linearized defect
	m_imVelocity.clear_lin_defect();

//  loop Sub Control Volume Faces (SCVF)
	for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
	{
	// get current SCVF
		const typename TFVGeom::SCVF& scvf = geo.scvf(ip);

	//	get derivatives to shapes at ip
		MathVector<dim>* vLinDef = m_imVelocity.lin_defect(ip, _C_);

	//	sum up contributions of convection shapes
		MathVector<dim> linDefect;
		VecSet(linDefect, 0.0);
		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
			VecScaleAppend(linDefect, u(_C_,sh), convShape.D_vel(ip, sh));

	//	add parts for both sides of scvf
		vLinDef[scvf.from()] += linDefect;
		vLinDef[scvf.to()] -= linDefect;
	}

//	we're done
	return true;
}

//	computes the linearized defect w.r.t to the velocity
template<typename TDomain>
template <typename TElem, typename TFVGeom>
bool ConvectionDiffusionElemDisc<TDomain>::
lin_defect_diffusion_fv1(const local_vector_type& u)
{
//  get finite volume geometry
	const static TFVGeom& geo = Provider::get<TFVGeom>();

//	get conv shapes
	const IConvectionShapes<dim>& convShape = get_updated_conv_shapes(geo);

//	cleat the linearized defect
	m_imDiffusion.clear_lin_defect();

//  loop Sub Control Volume Faces (SCVF)
	for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
	{
	// get current SCVF
		const typename TFVGeom::SCVF& scvf = geo.scvf(ip);

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
				MatAdd(linDefect, convShape.D_diffusion(ip, sh), u(_C_, sh));

	//	add contributions
		m_imDiffusion.lin_defect(ip, _C_, scvf.from()) -= linDefect;
		m_imDiffusion.lin_defect(ip, _C_, scvf.to()) += linDefect;
	}

//	we're done
	return true;
}

//	computes the linearized defect w.r.t to the reaction
template<typename TDomain>
template <typename TElem, typename TFVGeom>
bool ConvectionDiffusionElemDisc<TDomain>::
lin_defect_reaction_fv1(const local_vector_type& u)
{
//  get finite volume geometry
	const static TFVGeom& geo = Provider::get<TFVGeom>();

// 	loop Sub Control Volumes (SCV)
	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
	// 	get current SCV
		const typename TFVGeom::SCV& scv = geo.scv(ip);

	// 	get associated node
		const int co = scv.node_id();

	// 	Add to local defect
		m_imReaction.lin_defect(ip, _C_, co) = u(_C_, co) * scv.volume();
	}

//	we're done
	return true;
}

//	computes the linearized defect w.r.t to the source
template<typename TDomain>
template <typename TElem, typename TFVGeom>
bool ConvectionDiffusionElemDisc<TDomain>::
lin_defect_source_fv1(const local_vector_type& u)
{
//  get finite volume geometry
	const static TFVGeom& geo = Provider::get<TFVGeom>();

// 	loop Sub Control Volumes (SCV)
	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
	// 	get current SCV
		const typename TFVGeom::SCV& scv = geo.scv(ip);

	// 	get associated node
		const int co = scv.node_id();

	// 	Add to local defect
		m_imSource.lin_defect(ip, _C_, co) = scv.volume();
	}

//	we're done
	return true;
}

//	computes the linearized defect w.r.t to the mass scale
template<typename TDomain>
template <typename TElem, typename TFVGeom>
bool ConvectionDiffusionElemDisc<TDomain>::
lin_defect_mass_scale_fv1(const local_vector_type& u)
{
//  get finite volume geometry
	const static TFVGeom& geo = Provider::get<TFVGeom>();

// 	loop Sub Control Volumes (SCV)
	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
	// 	get current SCV
		const typename TFVGeom::SCV& scv = geo.scv(ip);

	// 	get associated node
		const int co = scv.node_id();

	// 	Add to local defect
		m_imMassScale.lin_defect(ip, _C_, co) = u(_C_, co) * scv.volume();
	}

//	we're done
	return true;
}

//	computes the linearized defect w.r.t to the velocity
template<typename TDomain>
template <typename TElem, typename TFVGeom>
bool ConvectionDiffusionElemDisc<TDomain>::
comp_export_concentration_fv1(const local_vector_type& u, bool bDeriv)
{
//  get finite volume geometry
	const static TFVGeom& geo = Provider::get<TFVGeom>();

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
				const typename TFVGeom::SCVF& scvf = geo.scvf(ip);

			//	compute concentration at ip
				number& cIP = m_exConcentration.value(s, ip);
				cIP = 0.0;
				for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
					cIP += u(_C_, sh) * scvf.shape(sh);

			//	compute derivative w.r.t. to unknowns iff needed
				if(bDeriv)
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
				if(bDeriv)
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
template<typename TDomain>
template <typename TElem, typename TFVGeom>
bool ConvectionDiffusionElemDisc<TDomain>::
comp_export_concentration_grad_fv1(const local_vector_type& u, bool bDeriv)
{
// 	Get finite volume geometry
	static const TFVGeom& geo = Provider::get<TFVGeom>();

//	get reference element dimension
	static const size_t refDim = TFVGeom::dim;

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
				const typename TFVGeom::SCVF& scvf = geo.scvf(ip);

			//	Compute Gradients and concentration at ip
				MathVector<dim>& cIP = m_exConcentrationGrad.value(s, ip);

				VecSet(cIP, 0.0);
				for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
					VecScaleAppend(cIP, u(_C_, sh), scvf.global_grad(sh));

				if(bDeriv)
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
template<typename TDomain>
const typename ConvectionDiffusionElemDisc<TDomain>::conv_shape_type&
ConvectionDiffusionElemDisc<TDomain>::
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


////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////

// register for all dim
template<typename TDomain>
void ConvectionDiffusionElemDisc<TDomain>::
register_all_fv1_funcs(bool bHang)
{
//	get all grid element types in this dimension and below
	typedef typename domain_traits<dim>::DimElemList ElemList;

//	switch assemble functions
	if(!bHang) boost::mpl::for_each<ElemList>( RegisterFV1<FV1Geometry>(this) );
	else boost::mpl::for_each<ElemList>( RegisterFV1<HFV1Geometry>(this) );

//	if(!bHang) boost::mpl::for_each<ElemList>( RegisterDimFV1<DimFV1Geometry>(this) );
//	else throw(UGFatalError("Hanging - DimFV1Geometry not implemented"));
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionElemDisc<TDomain>::
register_fv1_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef this_type T;

	reg_prepare_elem_loop_fct(id, &T::template elem_loop_prepare_fv1<TElem, TFVGeom>);
	reg_prepare_elem_fct(	  id, &T::template elem_prepare_fv1<TElem, TFVGeom>);
	reg_finish_elem_loop_fct( id, &T::template elem_loop_finish_fv1<TElem, TFVGeom>);
	reg_ass_JA_elem_fct(	  id, &T::template elem_JA_fv1<TElem, TFVGeom>);
	reg_ass_JM_elem_fct(	  id, &T::template elem_JM_fv1<TElem, TFVGeom>);
	reg_ass_dA_elem_fct(	  id, &T::template elem_dA_fv1<TElem, TFVGeom>);
	reg_ass_dM_elem_fct(	  id, &T::template elem_dM_fv1<TElem, TFVGeom>);
	reg_ass_rhs_elem_fct(	  id, &T::template elem_rhs_fv1<TElem, TFVGeom>);

//	set computation of linearized defect w.r.t velocity
	m_imVelocity. reg_lin_defect_fct(id, this, &T::template lin_defect_velocity_fv1<TElem, TFVGeom>);
	m_imDiffusion.reg_lin_defect_fct(id, this, &T::template lin_defect_diffusion_fv1<TElem, TFVGeom>);
	m_imReaction. reg_lin_defect_fct(id, this, &T::template lin_defect_reaction_fv1<TElem, TFVGeom>);
	m_imSource.	  reg_lin_defect_fct(id, this, &T::template lin_defect_source_fv1<TElem, TFVGeom>);
	m_imMassScale.reg_lin_defect_fct(id, this, &T::template lin_defect_mass_scale_fv1<TElem, TFVGeom>);

//	exports
	m_exConcentration.	  reg_export_fct(id, this, &T::template comp_export_concentration_fv1<TElem, TFVGeom>);
	m_exConcentrationGrad.reg_export_fct(id, this, &T::template comp_export_concentration_grad_fv1<TElem, TFVGeom>);
}

} // namespace ug

