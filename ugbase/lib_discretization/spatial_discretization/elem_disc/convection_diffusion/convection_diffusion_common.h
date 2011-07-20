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
//	register assemling functions
	register_all_fv1_funcs(false);

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

//	set default assembling to fe
	set_disc_scheme("fe");
}

template<typename TDomain>
size_t ConvectionDiffusionElemDisc<TDomain>::
num_fct(){return 1;}

template<typename TDomain>
bool ConvectionDiffusionElemDisc<TDomain>::
request_finite_element_id(const std::vector<LFEID>& vLfeID)
{
//	check number
	if(vLfeID.size() != num_fct()) return false;

//	check that Lagrange 1st order
	return vLfeID[0] == LFEID(LFEID::LAGRANGE, 1);
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
	if(m_discScheme == "fv") return true;
	else if(m_discScheme == "fe") return true;
	else throw(UGFatalError("Disc Scheme not recognized. Internal error."));
}

template<typename TDomain>
void ConvectionDiffusionElemDisc<TDomain>::
set_disc_scheme(const char* c_scheme)
{
//	convert to string
	std::string scheme = c_scheme;

//	check
	if(scheme != std::string("fe") && scheme != std::string("fv"))
	{
		UG_LOG("ERROR in 'ConvectionDiffusionElemDisc::set_disc_scheme':"
				" Only 'fe' and 'fv' supported.\n");
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
//	switch, which assemble functions to use; both supported.
	if(m_discScheme == "fv") register_all_fv1_funcs(m_bNonRegularGrid);
	else if(m_discScheme == "fe") register_all_fe1_funcs();
	else throw(UGFatalError("Disc Scheme not recognized. Internal error."));
}

} // namespace ug

