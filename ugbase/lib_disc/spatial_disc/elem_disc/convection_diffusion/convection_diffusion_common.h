/*
 * convection_diffusion_common.cpp
 *
 *  Created on: 26.02.2010
 *      Author: andreasvogel
 */

#include "convection_diffusion.h"
#include "lib_disc/spatial_disc/ip_data/const_user_data.h"
#include "lib_disc/spatial_disc/disc_util/conv_shape.h"
#include "bindings/lua/lua_user_data.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
//	user data
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
void ConvectionDiffusion<TDomain>::
set_upwind(SmartPtr<IConvectionShapes<dim> > shapes) {m_spConvShape = shapes;}

//////// Diffusion

template<typename TDomain>
void ConvectionDiffusion<TDomain>::
set_diffusion(SmartPtr<IPData<MathMatrix<dim, dim>, dim> > user) {m_imDiffusion.set_data(user);}

template<typename TDomain>
void ConvectionDiffusion<TDomain>::
set_diffusion(number val)
{
	set_diffusion(CreateSmartPtr(new ConstUserMatrix<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ConvectionDiffusion<TDomain>::
set_diffusion(const char* fctName)
{
	set_diffusion(CreateSmartPtr(new LuaUserData<MathMatrix<dim,dim>, dim>(fctName)));
}
#endif

//////// Velocity

template<typename TDomain>
void ConvectionDiffusion<TDomain>::
set_velocity(SmartPtr<IPData<MathVector<dim>, dim> > user) {m_imVelocity.set_data(user);}

template<typename TDomain>
void ConvectionDiffusion<TDomain>::
set_velocity(number vel_x)
{
	UG_THROW_FATAL("ConvectionDiffusion: Setting velocity vector of dimension 1"
					" to a Discretization for world dim " << dim);
}

template<>
void ConvectionDiffusion<Domain1d>::
set_velocity(number vel_x)
{
	SmartPtr<ConstUserVector<dim> > vel(new ConstUserVector<dim>());
	vel->set_entry(0, vel_x);
	set_velocity(vel);
}

template<typename TDomain>
void ConvectionDiffusion<TDomain>::
set_velocity(number vel_x, number vel_y)
{
	UG_THROW_FATAL("ConvectionDiffusion: Setting velocity vector of dimension 2"
					" to a Discretization for world dim " << dim);
}

template<>
void ConvectionDiffusion<Domain2d>::
set_velocity(number vel_x, number vel_y)
{
	SmartPtr<ConstUserVector<dim> > vel(new ConstUserVector<dim>());
	vel->set_entry(0, vel_x);
	vel->set_entry(1, vel_y);
	set_velocity(vel);
}

template<typename TDomain>
void ConvectionDiffusion<TDomain>::
set_velocity(number vel_x, number vel_y, number vel_z)
{
	UG_THROW_FATAL("ConvectionDiffusion: Setting velocity vector of dimension 3"
					" to a Discretization for world dim " << dim);
}

template<>
void ConvectionDiffusion<Domain3d>::
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
void ConvectionDiffusion<TDomain>::
set_velocity(const char* fctName)
{
	set_velocity(CreateSmartPtr(new LuaUserData<MathVector<dim>, dim>(fctName)));
}
#endif

//////// Reaction Rate

template<typename TDomain>
void ConvectionDiffusion<TDomain>::
set_reaction_rate(SmartPtr<IPData<number, dim> > user) {m_imReactionRate.set_data(user);}

template<typename TDomain>
void ConvectionDiffusion<TDomain>::
set_reaction_rate(number val)
{
	set_reaction_rate(CreateSmartPtr(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ConvectionDiffusion<TDomain>::
set_reaction_rate(const char* fctName)
{
	set_reaction_rate(CreateSmartPtr(new LuaUserData<number, dim>(fctName)));
}
#endif

//////// Reaction

template<typename TDomain>
void ConvectionDiffusion<TDomain>::
set_reaction(SmartPtr<IPData<number, dim> > user) {m_imReaction.set_data(user);}

template<typename TDomain>
void ConvectionDiffusion<TDomain>::
set_reaction(number val)
{
	set_reaction(CreateSmartPtr(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ConvectionDiffusion<TDomain>::
set_reaction(const char* fctName)
{
	set_reaction(CreateSmartPtr(new LuaUserData<number, dim>(fctName)));
}
#endif

//////// Source

template<typename TDomain>
void ConvectionDiffusion<TDomain>::
set_source(SmartPtr<IPData<number, dim> > user)	{m_imSource.set_data(user);}

template<typename TDomain>
void ConvectionDiffusion<TDomain>::
set_source(number val)
{
	set_source(CreateSmartPtr(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ConvectionDiffusion<TDomain>::
set_source(const char* fctName)
{
	set_source(CreateSmartPtr(new LuaUserData<number, dim>(fctName)));
}
#endif

//////// Mass Scale

template<typename TDomain>
void ConvectionDiffusion<TDomain>::
set_mass_scale(SmartPtr<IPData<number, dim> > user)	{m_imMassScale.set_data(user);}

template<typename TDomain>
void ConvectionDiffusion<TDomain>::
set_mass_scale(number val)
{
	set_mass_scale(CreateSmartPtr(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ConvectionDiffusion<TDomain>::
set_mass_scale(const char* fctName)
{
	set_mass_scale(CreateSmartPtr(new LuaUserData<number, dim>(fctName)));
}
#endif

//////// Mass

template<typename TDomain>
void ConvectionDiffusion<TDomain>::
set_mass(SmartPtr<IPData<number, dim> > user)	{m_imMass.set_data(user);}

template<typename TDomain>
void ConvectionDiffusion<TDomain>::
set_mass(number val)
{
	set_mass(CreateSmartPtr(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ConvectionDiffusion<TDomain>::
set_mass(const char* fctName)
{
	set_mass(CreateSmartPtr(new LuaUserData<number, dim>(fctName)));
}
#endif



template<typename TDomain>
bool ConvectionDiffusion<TDomain>::
time_point_changed(number time)
{
//	set new time point at imports
	m_imDiffusion.set_time(time);
	m_imVelocity.set_time(time);
	m_imSource.set_time(time);
	m_imReactionRate.set_time(time);
	m_imReaction.set_time(time);
	m_imMassScale.set_time(time);
	m_imMass.set_time(time);

//	this disc does not need the old time solutions, thus, return false
	return false;
}

template <typename TDomain>
typename ConvectionDiffusion<TDomain>::NumberExport
ConvectionDiffusion<TDomain>::
value() {return m_exConcentration;}


template <typename TDomain>
typename ConvectionDiffusion<TDomain>::GradExport
ConvectionDiffusion<TDomain>::
gradient() {return m_exConcentrationGrad;}

////////////////////////////////////////////////////////////////////////////////
//	general
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
ConvectionDiffusion<TDomain>::
ConvectionDiffusion(const char* functions, const char* subsets)
 : IDomainElemDisc<TDomain>(functions,subsets), m_spConvShape(new ConvectionShapesNoUpwind<dim>),
   m_exConcentration(new DataExport<number, dim>),
   m_exConcentrationGrad(new DataExport<MathVector<dim>, dim>)
{
//	check number of functions
	if(this->num_fct() != 1)
		UG_THROW_FATAL("Wrong number of functions: The ElemDisc 'ConvectionDiffusion'"
					   " needs exactly "<<1<<" symbolic function.");

//	register exports
	this->register_export(m_exConcentration);
	this->register_export(m_exConcentrationGrad);

//	register imports
	this->register_import(m_imDiffusion);
	this->register_import(m_imVelocity);
	this->register_import(m_imReactionRate);
	this->register_import(m_imReaction);
	this->register_import(m_imSource);
	this->register_import(m_imMassScale);
	this->register_import(m_imMass);

	m_imMassScale.set_mass_part(true);
	m_imMass.set_mass_part(true);
	m_imSource.set_rhs_part(true);

//	set defaults
	m_order = 1;
	m_bQuadOrderUserDef = false;
	m_quadOrder = -1;
	m_quadOrderSCV = -1;
	m_quadOrderSCVF = -1;
	m_bNonRegularGrid = false;
	m_discScheme = "fe";

//	update assemble functions
	set_ass_funcs();
}

template<typename TDomain>
bool ConvectionDiffusion<TDomain>::
request_finite_element_id(const std::vector<LFEID>& vLfeID)
{
//	check number
	if(vLfeID.size() != 1)
	{
		UG_LOG("ERROR in 'ConvectionDiffusion::request_finite_element_id':"
				" Wrong number of functions given. Need exactly "<<1<<"\n");
		return false;
	}

//	check that Lagrange order
	if(vLfeID[0].type() != LFEID::LAGRANGE)
	{
		UG_LOG("ERROR in 'ConvectionDiffusion::request_finite_element_id':"
				" Lagrange trial space needed.\n");
		return false;
	}

//	for fv only 1st order
	if(m_discScheme == "fv1" && vLfeID[0].order() != 1)
	{
		UG_LOG("ERROR in 'ConvectionDiffusion::request_finite_element_id':"
				" FV Scheme only implemented for 1st order.\n");
		return false;
	}

//	check that not ADAPTIVE
	if(vLfeID[0].order() < 1)
	{
		UG_LOG("ERROR in 'ConvectionDiffusion::request_finite_element_id':"
				" Adaptive or invalid order not implemented.\n");
		return false;
	}

//	remember lfeID;
	m_lfeID = vLfeID[0];

//	set order
	m_order = vLfeID[0].order();

//	update assemble functions
	set_ass_funcs();

//	is supported
	return true;
}

template<typename TDomain>
bool ConvectionDiffusion<TDomain>::
request_non_regular_grid(bool bNonRegular)
{
//	remember
	m_bNonRegularGrid = bNonRegular;

//	update assemble functions
	set_ass_funcs();

//	this disc supports both grids
	return true;
}

template<typename TDomain>
bool ConvectionDiffusion<TDomain>::
use_hanging() const
{
	if(m_discScheme == "fv1") return true;
	else if(m_discScheme == "fv") return false;
	else if(m_discScheme == "fe") return false;
	else throw(UGFatalError("Disc Scheme not recognized. Internal error."));
}

template<typename TDomain>
void ConvectionDiffusion<TDomain>::
set_disc_scheme(const char* c_scheme)
{
//	convert to string
	std::string scheme = c_scheme;

//	check
	if(scheme != std::string("fe") &&
	   scheme != std::string("fv1") &&
	   scheme != std::string("fv"))
	{
		UG_LOG("ERROR in 'ConvectionDiffusion::set_disc_scheme':"
				" Only 'fe', 'fv' and 'fv1' supported.\n");
	}

//	remember
	m_discScheme = scheme;

//	update assemble functions
	set_ass_funcs();
}

template<typename TDomain>
void ConvectionDiffusion<TDomain>::
set_ass_funcs()
{
//	set default quadrature order if not set by user
	if(!m_bQuadOrderUserDef)
	{
	//	FE
		m_quadOrder = 2* m_order + 1;

	//	FV
		m_quadOrderSCV = m_order;
		m_quadOrderSCVF = m_order;
	}
//	set all non-set orders
	else
	{
		if(m_quadOrder < 0) m_quadOrder = 2 * m_order + 1;
		if(m_quadOrderSCV < 0) m_quadOrderSCV = m_order;
		if(m_quadOrderSCVF < 0) m_quadOrderSCVF = m_order;
	}

//	switch, which assemble functions to use; both supported.
	if(m_discScheme == "fv1") register_all_fv1_funcs(m_bNonRegularGrid);
	else if(m_discScheme == "fv") register_all_fvho_funcs(m_order, m_quadOrderSCV, m_quadOrderSCVF);
	else if(m_discScheme == "fe") register_all_fe_funcs(m_order, m_quadOrder);
	else throw(UGFatalError("Disc Scheme not recognized. Internal error."));
}

} // namespace ug

