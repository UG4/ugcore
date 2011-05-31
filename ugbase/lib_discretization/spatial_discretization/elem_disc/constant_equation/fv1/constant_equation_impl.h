/*
 * constant_equation_impl.h
 *
 *  Created on: 13.05.2011
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__CONSTANT_EQUATION__FV1__CONSTANT_EQUATION_IMPL__
#define __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__CONSTANT_EQUATION__FV1__CONSTANT_EQUATION_IMPL__

#include "constant_equation.h"

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
FVConstantEquationElemDisc<TDomain, TAlgebra>::
prepare_element_loop()
{
	// all this will be performed outside of the loop over the elements.
	// Therefore it is not time critical.
	static const int refDim = TFVGeom<TElem, dim>::dim;

//	set local positions for rhs
	if(!TFVGeom<TElem, dim>::usesHangingNodes)
	{
		TFVGeom<TElem, dim>& geo = GeomProvider::get<TFVGeom<TElem,dim> >();
		m_imVelocity.template 	set_local_ips<refDim>(geo.scvf_local_ips(),
		                   	                      geo.num_scvf_ips());
		m_imSource.template 		set_local_ips<refDim>(geo.scv_local_ips(),
		               		                      geo.num_scv_ips());
		m_imMassScale.template set_local_ips<refDim>(geo.scv_local_ips(),
		                                           geo.num_scv_ips());
	}

//	we're done
	return true;
}

template<typename TDomain, typename TAlgebra>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
inline
bool
FVConstantEquationElemDisc<TDomain, TAlgebra>::
finish_element_loop()
{
//	nothing to do
	return true;
}

template<typename TDomain, typename TAlgebra>
bool
FVConstantEquationElemDisc<TDomain, TAlgebra>::
time_point_changed(number time)
{
//	set new time point at imports
	m_imVelocity.set_time(time);
	m_imSource.set_time(time);
	m_imMassScale.set_time(time);

//	this disc does not need the old time solutions, thus, return false
	return false;
}


template<typename TDomain, typename TAlgebra>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
inline
bool
FVConstantEquationElemDisc<TDomain, TAlgebra>::
prepare_element(TElem* elem, const local_vector_type& u,
								const local_index_type& glob_ind)
{
//	get dimension of reference element
	static const int refDim = TFVGeom<TElem, dim>::dim;

//	get corners
	m_vCornerCoords = this->template get_element_corners<TElem>(elem);

// 	Update Geometry for this element
	TFVGeom<TElem, dim>& geo = GeomProvider::get<TFVGeom<TElem,dim> >();
	if(!geo.update(elem, this->get_subset_handler(), &m_vCornerCoords[0]))
	{
		UG_LOG("FVConstantEquationElemDisc::prepare_element:"
				" Cannot update Finite Volume Geometry.\n"); return false;
	}

//	set local positions for rhs
	if(TFVGeom<TElem, dim>::usesHangingNodes)
	{
		m_imVelocity.template 	set_local_ips<refDim>(geo.scvf_local_ips(),
												  geo.num_scvf_ips());
		m_imSource.template 		set_local_ips<refDim>(geo.scv_local_ips(),
												  geo.num_scv_ips());
		m_imMassScale.template set_local_ips<refDim>(geo.scv_local_ips(),
												   geo.num_scv_ips());
	}

//	set global positions for rhs
	m_imVelocity.set_global_ips(geo.scvf_global_ips(), geo.num_scvf_ips());
	m_imSource.set_global_ips(geo.scv_global_ips(), geo.num_scv_ips());
	m_imMassScale.set_global_ips(geo.scv_global_ips(), geo.num_scv_ips());

//	we're done
	return true;
}

template<typename TDomain, typename TAlgebra>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
inline
bool
FVConstantEquationElemDisc<TDomain, TAlgebra>::
assemble_JA(local_matrix_type& J, const local_vector_type& u)
{
//	no own contribution to jacobian (only due to linearized defect)
	return true;
}


template<typename TDomain, typename TAlgebra>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
inline
bool
FVConstantEquationElemDisc<TDomain, TAlgebra>::
assemble_JM(local_matrix_type& J, const local_vector_type& u)
{
//	no own contribution to jacobian (only due to linearized defect)
	return true;
}


template<typename TDomain, typename TAlgebra>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
inline
bool
FVConstantEquationElemDisc<TDomain, TAlgebra>::
assemble_A(local_vector_type& d, const local_vector_type& u)
{
// 	get finite volume geometry
	const static TFVGeom<TElem, dim>& geo	= GeomProvider::get<TFVGeom<TElem,dim> >();

//	check if data given
	if(!m_imVelocity.data_given()) return true;

// 	loop Sub Control Volume Faces (SCVF)
	for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
	{
	// 	get current SCVF
		const typename TFVGeom<TElem, dim>::SCVF& scvf = geo.scvf(ip);

	//	sum up convective flux using convection shapes
		const number conv_flux = VecDot(m_imVelocity[ip], scvf.normal());

	//  add to local defect
		d(_C_, scvf.from()) += conv_flux;
		d(_C_, scvf.to()  ) -= conv_flux;
	}

// 	we're done
	return true;
}


template<typename TDomain, typename TAlgebra>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
inline
bool
FVConstantEquationElemDisc<TDomain, TAlgebra>::
assemble_M(local_vector_type& d, const local_vector_type& u)
{
// 	get finite volume geometry
	const static TFVGeom<TElem, dim>& geo = GeomProvider::get<TFVGeom<TElem,dim> >();

// 	loop Sub Control Volumes (SCV)
	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
	// 	get current SCV
		const typename TFVGeom<TElem, dim>::SCV& scv = geo.scv(ip);

	// 	get associated node
		const int co = scv.node_id();

	//	mass value
		number val = scv.volume();

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
FVConstantEquationElemDisc<TDomain, TAlgebra>::
assemble_f(local_vector_type& d)
{
//	if zero data given, return
	if(!m_imSource.data_given()) return true;

// 	get finite volume geometry
	const static TFVGeom<TElem, dim>& geo
		= GeomProvider::get<TFVGeom<TElem,dim> >();

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
FVConstantEquationElemDisc<TDomain, TAlgebra>::
lin_defect_velocity(const local_vector_type& u)
{
//  get finite volume geometry
	const static TFVGeom<TElem, dim>& geo = GeomProvider::get<TFVGeom<TElem,dim> >();

//	reset the values for the linearized defect
	m_imVelocity.clear_lin_defect();

//  loop Sub Control Volume Faces (SCVF)
	for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
	{
	// get current SCVF
		const typename TFVGeom<TElem, dim>::SCVF& scvf = geo.scvf(ip);

	//	get derivatives to shapes at ip
		MathVector<dim>* vLinDef = m_imVelocity.lin_defect(ip, _C_);

	//	add parts for both sides of scvf
		vLinDef[scvf.from()] += scvf.normal();
		vLinDef[scvf.to()] -= scvf.normal();
	}

//	we're done
	return true;
}

//	computes the linearized defect w.r.t to the source
template<typename TDomain, typename TAlgebra>
template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
bool
FVConstantEquationElemDisc<TDomain, TAlgebra>::
lin_defect_source(const local_vector_type& u)
{
//  get finite volume geometry
	const static TFVGeom<TElem, dim>& geo = GeomProvider::get<TFVGeom<TElem,dim> >();

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
FVConstantEquationElemDisc<TDomain, TAlgebra>::
lin_defect_mass_scale(const local_vector_type& u)
{
// 	get finite volume geometry
	const static TFVGeom<TElem, dim>& geo	= GeomProvider::get<TFVGeom<TElem,dim> >();

//	reset all values
	m_imMassScale.clear_lin_defect();

// 	loop Sub Control Volumes (SCV)
	for(size_t co = 0; co < geo.num_scv(); ++co)
	{
	// 	get current SCV
		const typename TFVGeom<TElem, dim>::SCV& scv = geo.scv(co);

	// 	Check associated node
		UG_ASSERT(co == scv.node_id(), "Only one shape per SCV");

	// 	Add to local defect
		m_imMassScale.lin_defect(co, _C_, co) = scv.volume();
	}

//	we're done
	return true;
}

//	computes the linearized defect w.r.t to the velocity
template<typename TDomain, typename TAlgebra>
template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
bool
FVConstantEquationElemDisc<TDomain, TAlgebra>::
compute_concentration_export(const local_vector_type& u, bool compDeriv)
{
//  get finite volume geometry
	const static TFVGeom<TElem, dim>& geo = GeomProvider::get<TFVGeom<TElem,dim> >();

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
FVConstantEquationElemDisc<TDomain, TAlgebra>::
compute_concentration_grad_export(const local_vector_type& u, bool compDeriv)
{
// 	Get finite volume geometry
	static const TFVGeom<TElem, dim>& geo = GeomProvider::get<TFVGeom<TElem,dim> >();

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


} // namespace ug


#endif /*__H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__CONSTANT_EQUATION__FV1__CONSTANT_EQUATION_IMPL__*/
