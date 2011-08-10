/*
 * convection_diffusion_common.cpp
 *
 *  Created on: 26.02.2010
 *      Author: andreasvogel
 */

#include "convection_diffusion.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
//	user data
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
void ConvectionDiffusionElemDisc<TDomain>::
set_upwind(IConvectionShapes<dim>& shapes) {m_pConvShape = &shapes;}

template<typename TDomain>
void ConvectionDiffusionElemDisc<TDomain>::
set_diffusion(IPData<MathMatrix<dim, dim>, dim>& user) {m_imDiffusion.set_data(user);}

template<typename TDomain>
void ConvectionDiffusionElemDisc<TDomain>::
set_velocity(IPData<MathVector<dim>, dim>& user) {m_imVelocity.set_data(user);}

template<typename TDomain>
void ConvectionDiffusionElemDisc<TDomain>::
set_reaction(IPData<number, dim>& user) {m_imReaction.set_data(user);}

template<typename TDomain>
void ConvectionDiffusionElemDisc<TDomain>::
set_source(IPData<number, dim>& user)	{m_imSource.set_data(user);}

template<typename TDomain>
void ConvectionDiffusionElemDisc<TDomain>::
set_mass_scale(IPData<number, dim>& user)	{m_imMassScale.set_data(user);}

template<typename TDomain>
bool ConvectionDiffusionElemDisc<TDomain>::
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

template <typename TDomain>
typename ConvectionDiffusionElemDisc<TDomain>::NumberExport &
ConvectionDiffusionElemDisc<TDomain>::
get_concentration() {return m_exConcentration;}


template <typename TDomain>
typename ConvectionDiffusionElemDisc<TDomain>::GradExport &
ConvectionDiffusionElemDisc<TDomain>::
get_concentration_grad() {return m_exConcentrationGrad;}

////////////////////////////////////////////////////////////////////////////////
//	general
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
ConvectionDiffusionElemDisc<TDomain>::
ConvectionDiffusionElemDisc()
 : m_pConvShape(NULL)
{
//	register exports
	register_export(m_exConcentration);
	register_export(m_exConcentrationGrad);

//	register imports
	register_import(m_imDiffusion);
	register_import(m_imVelocity);
	register_import(m_imReaction);
	register_import(m_imSource);
	register_import(m_imMassScale);

	m_imMassScale.set_mass_part(true);

//	set defaults
	m_order = 1;
	m_bQuadOrderUserDef = false;
	m_quadOrder = -1;
	m_quadOrderSCV = -1;
	m_quadOrderSCVF = -1;
	m_bNonRegularGrid = false;
	m_discScheme = "fe";

//	update assemble functions
	set_assemble_funcs();
}

template<typename TDomain>
size_t ConvectionDiffusionElemDisc<TDomain>::
num_fct(){return 1;}

template<typename TDomain>
bool ConvectionDiffusionElemDisc<TDomain>::
request_finite_element_id(const std::vector<LFEID>& vLfeID)
{
//	check number
	if(vLfeID.size() != num_fct())
	{
		UG_LOG("ERROR in 'ConvectionDiffusionElemDisc::request_finite_element_id':"
				" Wrong number of functions given. Need exactly "<<num_fct()<<"\n");
		return false;
	}

//	check that Lagrange order
	if(vLfeID[0].type() != LFEID::LAGRANGE)
	{
		UG_LOG("ERROR in 'ConvectionDiffusionElemDisc::request_finite_element_id':"
				" Lagrange trial space needed.\n");
		return false;
	}

//	for fv only 1st order
	if(m_discScheme == "fv1" && vLfeID[0].order() != 1)
	{
		UG_LOG("ERROR in 'ConvectionDiffusionElemDisc::request_finite_element_id':"
				" FV Scheme only implemented for 1st order.\n");
		return false;
	}

//	check that not ADAPTIVE
	if(vLfeID[0].order() < 1)
	{
		UG_LOG("ERROR in 'ConvectionDiffusionElemDisc::request_finite_element_id':"
				" Adaptive or invalid order not implemented.\n");
		return false;
	}

//	remember lfeID;
	m_lfeID = vLfeID[0];

//	set order
	m_order = vLfeID[0].order();

//	update assemble functions
	set_assemble_funcs();

//	is supported
	return true;
}

template<typename TDomain>
bool ConvectionDiffusionElemDisc<TDomain>::
treat_non_regular_grid(bool bNonRegular)
{
//	remember
	m_bNonRegularGrid = bNonRegular;

//	update assemble functions
	set_assemble_funcs();

//	this disc supports both grids
	return true;
}

template<typename TDomain>
bool ConvectionDiffusionElemDisc<TDomain>::
use_hanging() const
{
	if(m_discScheme == "fv1") return true;
	else if(m_discScheme == "fv") return false;
	else if(m_discScheme == "fe") return false;
	else throw(UGFatalError("Disc Scheme not recognized. Internal error."));
}

template<typename TDomain>
void ConvectionDiffusionElemDisc<TDomain>::
set_disc_scheme(const char* c_scheme)
{
//	convert to string
	std::string scheme = c_scheme;

//	check
	if(scheme != std::string("fe") &&
	   scheme != std::string("fv1") &&
	   scheme != std::string("fv"))
	{
		UG_LOG("ERROR in 'ConvectionDiffusionElemDisc::set_disc_scheme':"
				" Only 'fe', 'fv' and 'fv1' supported.\n");
	}

//	remember
	m_discScheme = scheme;

//	update assemble functions
	set_assemble_funcs();
}

template<typename TDomain>
void ConvectionDiffusionElemDisc<TDomain>::
set_assemble_funcs()
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

