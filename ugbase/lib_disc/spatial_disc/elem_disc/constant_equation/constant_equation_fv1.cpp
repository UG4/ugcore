/*
 * constant_equation_fv1.cpp
 *
 *  Created on: 13.05.2011
 *      Author: andreasvogel
 */

#include "constant_equation.h"

#include "common/util/provider.h"
#include "lib_disc/spatial_disc/disc_util/finite_volume_geometry.h"
#include "lib_disc/spatial_disc/disc_util/hanging_finite_volume_geometry.h"
#include "lib_disc/spatial_disc/ip_data/const_user_data.h"
#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

namespace ug{


////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
//  Provide a generic implementation for all elements
//  (since this discretization can be implemented in a generic way)
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

//////// Velocity

template<typename TDomain>
void ConstantEquation<TDomain>::
set_velocity(SmartPtr<IPData<MathVector<dim>, dim> > user) {m_imVelocity.set_data(user);}

template<typename TDomain>
void ConstantEquation<TDomain>::
set_velocity(number vel_x)
{
	UG_THROW_FATAL("ConstantEquation: Setting velocity vector of dimension 1"
					" to a Discretization for world dim " << dim);
}

template<>
void ConstantEquation<Domain1d>::
set_velocity(number vel_x)
{
	SmartPtr<ConstUserVector<dim> > vel(new ConstUserVector<dim>());
	vel->set_entry(0, vel_x);
	set_velocity(vel);
}

template<typename TDomain>
void ConstantEquation<TDomain>::
set_velocity(number vel_x, number vel_y)
{
	UG_THROW_FATAL("ConstantEquation: Setting velocity vector of dimension 2"
					" to a Discretization for world dim " << dim);
}

template<>
void ConstantEquation<Domain2d>::
set_velocity(number vel_x, number vel_y)
{
	SmartPtr<ConstUserVector<dim> > vel(new ConstUserVector<dim>());
	vel->set_entry(0, vel_x);
	vel->set_entry(1, vel_y);
	set_velocity(vel);
}

template<typename TDomain>
void ConstantEquation<TDomain>::
set_velocity(number vel_x, number vel_y, number vel_z)
{
	UG_THROW_FATAL("ConstantEquation: Setting velocity vector of dimension 3"
					" to a Discretization for world dim " << dim);
}

template<>
void ConstantEquation<Domain3d>::
set_velocity(number vel_x, number vel_y, number vel_z)
{
	SmartPtr<ConstUserVector<dim> > vel(new ConstUserVector<dim>());
	vel->set_entry(0, vel_x);
	vel->set_entry(1, vel_y);
	vel->set_entry(2, vel_z);
	set_velocity(vel);
}


#ifdef UG_FOR_LUA
template<typename TDomain>
void ConstantEquation<TDomain>::
set_velocity(const char* fctName)
{
	set_velocity(LuaUserDataFactory<MathVector<dim>, dim>::create(fctName));
}
#endif

//////// Source

template<typename TDomain>
void ConstantEquation<TDomain>::
set_source(SmartPtr<IPData<number, dim> > user)	{m_imSource.set_data(user);}

template<typename TDomain>
void ConstantEquation<TDomain>::
set_source(number val)
{
	set_source(CreateSmartPtr(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ConstantEquation<TDomain>::
set_source(const char* fctName)
{
	set_source(LuaUserDataFactory<number, dim>::create(fctName));
}
#endif

//////// Mass

template<typename TDomain>
void ConstantEquation<TDomain>::
set_mass(SmartPtr<IPData<number, dim> > user)	{m_imMass.set_data(user);}

template<typename TDomain>
void ConstantEquation<TDomain>::
set_mass(number val)
{
	set_mass(CreateSmartPtr(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ConstantEquation<TDomain>::
set_mass(const char* fctName)
{
	set_mass(LuaUserDataFactory<number, dim>::create(fctName));
}
#endif


template<typename TDomain>
ConstantEquation<TDomain>::
ConstantEquation(const char* functions, const char* subsets)
:IDomainElemDisc<TDomain>(functions,subsets),
 m_exConcentration(new DataExport<number, dim>),
 m_exConcentrationGrad(new DataExport<MathVector<dim>, dim>)
{
//	check number of functions
	if(this->num_fct() != 1)
		UG_THROW_FATAL("Wrong number of functions: The ElemDisc 'ConstantEquation'"
					   " needs exactly "<<1<<" symbolic function.");

//	register assemling functions
	register_all_fv1_funcs(false);

//	register exports
	this->register_export(m_exConcentration);
	this->register_export(m_exConcentrationGrad);

//	register imports
	this->register_import(m_imVelocity);
	this->register_import(m_imSource);
	this->register_import(m_imMass);

	m_imMass.set_mass_part(true);
	m_imSource.set_rhs_part(true);
}

template<typename TDomain>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
void ConstantEquation<TDomain>::prepare_element_loop()
{
	// all this will be performed outside of the loop over the elements.
	// Therefore it is not time critical.
	static const int refDim = TFVGeom<TElem, dim>::dim;

//	set local positions for rhs
	if(!TFVGeom<TElem, dim>::usesHangingNodes)
	{
		TFVGeom<TElem, dim>& geo = Provider<TFVGeom<TElem,dim> >::get();
		m_imVelocity.template  set_local_ips<refDim>(geo.scvf_local_ips(), geo.num_scvf_ips(), false);
		m_imSource.template    set_local_ips<refDim>(geo.scv_local_ips(), geo.num_scv_ips(), false);
		m_imMass.template set_local_ips<refDim>(geo.scv_local_ips(), geo.num_scv_ips(), false);
	}
}

template<typename TDomain>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
void ConstantEquation<TDomain>::finish_element_loop()
{}


template<typename TDomain>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
void ConstantEquation<TDomain>::
prepare_element(TElem* elem, const LocalVector& u)
{
//	get dimension of reference element
	static const int refDim = TFVGeom<TElem, dim>::dim;

//	get corners
	m_vCornerCoords = this->template element_corners<TElem>(elem);

// 	Update Geometry for this element
	TFVGeom<TElem, dim>& geo = Provider<TFVGeom<TElem,dim> >::get();
	if(!geo.update(elem, &m_vCornerCoords[0], &(this->subset_handler())))
		UG_THROW_FATAL("ConstantEquation::prepare_element:"
						" Cannot update Finite Volume Geometry.");

//	set local positions for rhs
	if(TFVGeom<TElem, dim>::usesHangingNodes)
	{
		m_imVelocity.template  set_local_ips<refDim>(geo.scvf_local_ips(), geo.num_scvf_ips());
		m_imSource.template    set_local_ips<refDim>(geo.scv_local_ips(), geo.num_scv_ips());
		m_imMass.template set_local_ips<refDim>(geo.scv_local_ips(), geo.num_scv_ips());
	}

//	set global positions for rhs
	m_imVelocity.set_global_ips(geo.scvf_global_ips(), geo.num_scvf_ips());
	m_imSource.set_global_ips(geo.scv_global_ips(), geo.num_scv_ips());
	m_imMass.set_global_ips(geo.scv_global_ips(), geo.num_scv_ips());
}

template<typename TDomain>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
void ConstantEquation<TDomain>::
ass_JA_elem(LocalMatrix& J, const LocalVector& u)
{
//	no own contribution to jacobian (only due to linearized defect)
}


template<typename TDomain>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
void ConstantEquation<TDomain>::
ass_JM_elem(LocalMatrix& J, const LocalVector& u)
{
//	no own contribution to jacobian (only due to linearized defect)
}


template<typename TDomain>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
void ConstantEquation<TDomain>::
ass_dA_elem(LocalVector& d, const LocalVector& u)
{
// 	get finite volume geometry
	const static TFVGeom<TElem, dim>& geo	= Provider<TFVGeom<TElem,dim> >::get();

//	check if data given
	if(!m_imVelocity.data_given()) return;

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
}


template<typename TDomain>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
void ConstantEquation<TDomain>::
ass_dM_elem(LocalVector& d, const LocalVector& u)
{
// 	get finite volume geometry
	const static TFVGeom<TElem, dim>& geo = Provider<TFVGeom<TElem,dim> >::get();

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
		if(m_imMass.data_given())
			val *= m_imMass[ip];

	// 	Add to local defect
		d(_C_, co) += val;
	}
}


template<typename TDomain>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
void ConstantEquation<TDomain>::
ass_rhs_elem(LocalVector& d)
{
//	if zero data given, return
	if(!m_imSource.data_given()) return;

// 	get finite volume geometry
	const static TFVGeom<TElem, dim>& geo
		= Provider<TFVGeom<TElem,dim> >::get();

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
}


//	computes the linearized defect w.r.t to the velocity
template<typename TDomain>
template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
void ConstantEquation<TDomain>::
lin_def_velocity(const LocalVector& u,
                 std::vector<std::vector<MathVector<dim> > > vvvLinDef[],
                 const size_t nip)
{
//  get finite volume geometry
	const static TFVGeom<TElem, dim>& geo = Provider<TFVGeom<TElem,dim> >::get();

//	reset the values for the linearized defect
	for(size_t ip = 0; ip < nip; ++ip)
		for(size_t c = 0; c < vvvLinDef[ip].size(); ++c)
			for(size_t sh = 0; sh < vvvLinDef[ip][c].size(); ++sh)
				vvvLinDef[ip][c][sh] = 0.0;

//  loop Sub Control Volume Faces (SCVF)
	for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
	{
	// get current SCVF
		const typename TFVGeom<TElem, dim>::SCVF& scvf = geo.scvf(ip);

	//	add parts for both sides of scvf
		vvvLinDef[ip][_C_][scvf.from()] += scvf.normal();
		vvvLinDef[ip][_C_][scvf.to()] -= scvf.normal();
	}
}

//	computes the linearized defect w.r.t to the source
template<typename TDomain>
template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
void ConstantEquation<TDomain>::
lin_def_source(const LocalVector& u,
               std::vector<std::vector<number> > vvvLinDef[],
               const size_t nip)
{
//  get finite volume geometry
	const static TFVGeom<TElem, dim>& geo = Provider<TFVGeom<TElem,dim> >::get();

// 	loop Sub Control Volumes (SCV)
	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
	// 	get current SCV
		const typename TFVGeom<TElem, dim>::SCV& scv = geo.scv(ip);

	// 	get associated node
		const int co = scv.node_id();

	// 	set lin defect
		vvvLinDef[ip][_C_][co] = scv.volume();
	}
}

//	computes the linearized defect w.r.t to the mass scale
template<typename TDomain>
template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
void ConstantEquation<TDomain>::
lin_def_mass(const LocalVector& u,
             std::vector<std::vector<number> > vvvLinDef[],
             const size_t nip)
{
// 	get finite volume geometry
	const static TFVGeom<TElem, dim>& geo	= Provider<TFVGeom<TElem,dim> >::get();

// 	loop Sub Control Volumes (SCV)
	for(size_t co = 0; co < geo.num_scv(); ++co)
	{
	// 	get current SCV
		const typename TFVGeom<TElem, dim>::SCV& scv = geo.scv(co);

	// 	Check associated node
		UG_ASSERT(co == scv.node_id(), "Only one shape per SCV");

	// 	set lin defect
		vvvLinDef[co][_C_][co] = scv.volume();
	}
}

//	computes the linearized defect w.r.t to the velocity
template<typename TDomain>
template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
void ConstantEquation<TDomain>::
ex_concentration(const LocalVector& u,
                 const MathVector<dim> vGlobIP[],
                 const MathVector<TFVGeom<TElem, dim>::dim> vLocIP[],
                 const size_t nip,
                 number vValue[],
                 bool bDeriv,
                 std::vector<std::vector<number> > vvvDeriv[])
{
//  get finite volume geometry
	const static TFVGeom<TElem, dim>& geo = Provider<TFVGeom<TElem,dim> >::get();

//	reference element
	typedef typename reference_element_traits<TElem>::reference_element_type
			ref_elem_type;

//	number of shape functions
	static const size_t numSH =	ref_elem_type::num_corners;

//	FV1 SCVF ip
	if(vLocIP == geo.scvf_local_ips())
	{
	//	Loop Sub Control Volume Faces (SCVF)
		for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
		{
		// 	Get current SCVF
			const typename TFVGeom<TElem, dim>::SCVF& scvf = geo.scvf(ip);

		//	compute concentration at ip
			vValue[ip] = 0.0;
			for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
				vValue[ip] += u(_C_, sh) * scvf.shape(sh);

		//	compute derivative w.r.t. to unknowns iff needed
			if(bDeriv)
				for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
					vvvDeriv[ip][_C_][sh] = scvf.shape(sh);
		}
	}
//	FV1 SCV ip
	else if(vLocIP == geo.scv_local_ips())
	{
	//	solution at ip
		for(size_t sh = 0; sh < numSH; ++sh)
			vValue[sh] = u(_C_, sh);

	//	set derivatives if needed
		if(bDeriv)
			for(size_t sh = 0; sh < numSH; ++sh)
				for(size_t sh2 = 0; sh2 < numSH; ++sh2)
					vvvDeriv[sh][_C_][sh2] = (sh==sh2) ? 1.0 : 0.0;
	}
// 	general case
	else
	{
	//	get trial space
		LagrangeP1<ref_elem_type> rTrialSpace = Provider<LagrangeP1<ref_elem_type> >::get();

	//	storage for shape function at ip
		number vShape[numSH];

	//	loop ips
		for(size_t ip = 0; ip < nip; ++ip)
		{
		//	evaluate at shapes at ip
			rTrialSpace.shapes(vShape, vLocIP[ip]);

		//	compute concentration at ip
			vValue[ip] = 0.0;
			for(size_t sh = 0; sh < numSH; ++sh)
				vValue[ip] += u(_C_, sh) * vShape[sh];

		//	compute derivative w.r.t. to unknowns iff needed
		//	\todo: maybe store shapes directly in vvvDeriv
			if(bDeriv)
				for(size_t sh = 0; sh < numSH; ++sh)
					vvvDeriv[ip][_C_][sh] = vShape[sh];
		}
	}
}

//	computes the linearized defect w.r.t to the velocity
template<typename TDomain>
template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
void ConstantEquation<TDomain>::
ex_concentration_grad(const LocalVector& u,
                      const MathVector<dim> vGlobIP[],
                      const MathVector<TFVGeom<TElem, dim>::dim> vLocIP[],
                      const size_t nip,
                      MathVector<dim> vValue[],
                      bool bDeriv,
                      std::vector<std::vector<MathVector<dim> > > vvvDeriv[])
{
// 	Get finite volume geometry
	static const TFVGeom<TElem, dim>& geo = Provider<TFVGeom<TElem,dim> >::get();

//	reference element
	typedef typename reference_element_traits<TElem>::reference_element_type
			ref_elem_type;

//	reference dimension
	static const int refDim = ref_elem_type::dim;

//	number of shape functions
	static const size_t numSH =	ref_elem_type::num_corners;

//	FV1 SCVF ip
	if(vLocIP == geo.scvf_local_ips())
	{
	//	Loop Sub Control Volume Faces (SCVF)
		for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
		{
		// 	Get current SCVF
			const typename TFVGeom<TElem, dim>::SCVF& scvf = geo.scvf(ip);

			VecSet(vValue[ip], 0.0);
			for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
				VecScaleAppend(vValue[ip], u(_C_, sh), scvf.global_grad(sh));

			if(bDeriv)
				for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
					vvvDeriv[ip][_C_][sh] = scvf.global_grad(sh);
		}
	}
// 	general case
	else
	{
	//	get trial space
		LagrangeP1<ref_elem_type>& rTrialSpace = Provider<LagrangeP1<ref_elem_type> >::get();

	//	storage for shape function at ip
		MathVector<refDim> vLocGrad[numSH];
		MathVector<refDim> locGrad;

	//	Reference Mapping
		MathMatrix<dim, refDim> JTInv;
		ReferenceMapping<ref_elem_type, dim> mapping(&m_vCornerCoords[0]);

	//	loop ips
		for(size_t ip = 0; ip < nip; ++ip)
		{
		//	evaluate at shapes at ip
			rTrialSpace.grads(vLocGrad, vLocIP[ip]);

		//	compute grad at ip
			VecSet(locGrad, 0.0);
			for(size_t sh = 0; sh < numSH; ++sh)
				VecScaleAppend(locGrad, u(_C_, sh), vLocGrad[sh]);

		//	compute global grad
			mapping.jacobian_transposed_inverse(JTInv, vLocIP[ip]);
			MatVecMult(vValue[ip], JTInv, locGrad);

		//	compute derivative w.r.t. to unknowns iff needed
			if(bDeriv)
				for(size_t sh = 0; sh < numSH; ++sh)
					MatVecMult(vvvDeriv[ip][_C_][sh], JTInv, vLocGrad[sh]);
		}
	}
};


////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////

// register for 1D
template<typename TDomain>
void
ConstantEquation<TDomain>::
register_all_fv1_funcs(bool bHang)
{
//	get all grid element types in this dimension and below
	typedef typename domain_traits<dim>::AllElemList ElemList;

//	switch assemble functions
	if(!bHang) boost::mpl::for_each<ElemList>( RegisterFV1<FV1Geometry>(this) );
	else boost::mpl::for_each<ElemList>( RegisterFV1<HFV1Geometry>(this) );
}

template<typename TDomain>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
void
ConstantEquation<TDomain>::
register_fv1_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef this_type T;
	static const int refDim = reference_element_traits<TElem>::dim;

	this->enable_fast_ass_elem(true);
	this->set_prep_elem_loop_fct(	id, &T::template prepare_element_loop<TElem, TFVGeom>);
	this->set_prep_elem_fct(		id, &T::template prepare_element<TElem, TFVGeom>);
	this->set_fsh_elem_loop_fct( 	id, &T::template finish_element_loop<TElem, TFVGeom>);
	this->set_ass_JA_elem_fct(	id, &T::template ass_JA_elem<TElem, TFVGeom>);
	this->set_ass_JM_elem_fct(	id, &T::template ass_JM_elem<TElem, TFVGeom>);
	this->set_ass_dA_elem_fct(	id, &T::template ass_dA_elem<TElem, TFVGeom>);
	this->set_ass_dM_elem_fct(	id, &T::template ass_dM_elem<TElem, TFVGeom>);
	this->set_ass_rhs_elem_fct(	id, &T::template ass_rhs_elem<TElem, TFVGeom>);

//	set computation of linearized defect w.r.t velocity
	m_imVelocity. set_fct(id, this, &T::template lin_def_velocity<TElem, TFVGeom>);
	m_imSource.	  set_fct(id, this, &T::template lin_def_source<TElem, TFVGeom>);
	m_imMass.set_fct(id, this, &T::template lin_def_mass<TElem, TFVGeom>);

//	exports
	m_exConcentration->	   template set_fct<T,refDim>(id, this, &T::template ex_concentration<TElem, TFVGeom>);
	m_exConcentrationGrad->template set_fct<T,refDim>(id, this, &T::template ex_concentration_grad<TElem, TFVGeom>);
}

////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

template class ConstantEquation<Domain1d>;
template class ConstantEquation<Domain2d>;
template class ConstantEquation<Domain3d>;


} // namespace ug

