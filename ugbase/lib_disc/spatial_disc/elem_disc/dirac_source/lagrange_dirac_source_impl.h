/*
 * Copyright (c) 2019-2020:  G-CSC, Goethe University Frankfurt
 * Author: Arne Naegel
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef __H__UG__LIB_DISC__SPACIAL_DISCRETIZATION__ELEM_DISC__DIRAC_SOURCE__LAGRANGE_DIRAC_SOURCE_IMPL__
#define __H__UG__LIB_DISC__SPACIAL_DISCRETIZATION__ELEM_DISC__DIRAC_SOURCE__LAGRANGE_DIRAC_SOURCE_IMPL__

#include "lagrange_dirac_source.h"
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"
#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"



namespace ug
{

template<typename TDomain>
void DiracSourceDisc<TDomain>::
add_source(number scale, MathVector<dim> &coord)
{
	add_source(make_sp(new ConstUserNumber<dim>(scale)), coord);
}


template<typename TDomain>
void DiracSourceDisc<TDomain>::
add_source(SmartPtr<UserData<number, dim> > srcData, MathVector<dim> &srcCoord)
{
	m_srcData = srcData;
	m_srcCoord = srcCoord;
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void DiracSourceDisc<TDomain>::
add_source(const char* luaFctName, MathVector<dim> &coord)
{
	add_source(LuaUserDataFactory<number,dim>::create(luaFctName), coord);
}
#endif



template<typename TDomain>
void DiracSourceDisc<TDomain>::
prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid)
{
	// check that Lagrange 1st order
	for(size_t i = 0; i < vLfeID.size(); ++i)
		if(vLfeID[i].type() != LFEID::LAGRANGE )
			UG_THROW("DiracSourceDisc: Lagrange polynomials expected.");

	// update assemble functions
	m_bNonRegularGrid = bNonRegularGrid;
	register_all_funcs();
}

template<typename TDomain>
bool DiracSourceDisc<TDomain>::
use_hanging() const
{
	return true;
}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void DiracSourceDisc<TDomain>::
add_rhs_elem(LocalVector& rhs, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	UG_LOG("'DiracSourceDisc::add_rhs_elem' called for " << vCornerCoords[0]);
	UG_LOG("elem:"<<  elem << "\n");

	UG_ASSERT(vCornerCoords[0] == m_srcCoord, "Source must be located at element corner!");

	// Request source strength form user data.
	const int co = 0;
	double srcVal;
	(*m_srcData)(srcVal, m_srcCoord, this->time (), -1);
	 rhs(0, co) += srcVal;

}





////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////

// register for 1D
template<typename TDomain>
void DiracSourceDisc<TDomain>::
register_all_funcs()
{


//	switch assemble functions
	// boost::mpl::for_each<ElemList>(RegisterFV1(this));
	register_func<RegularVertex, FV1Geometry<RegularVertex, dim> >();
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void DiracSourceDisc<TDomain>::
register_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef this_type T;

	this->clear_add_fct(id);
	this->set_prep_elem_loop_fct(	id, &T::template prep_elem_loop<TElem, TFVGeom>);
	this->set_prep_elem_fct(	 	id, &T::template prep_elem<TElem, TFVGeom>);
	this->set_fsh_elem_loop_fct( 	id, &T::template fsh_elem_loop<TElem, TFVGeom>);
	this->set_add_jac_A_elem_fct(	id, &T::template add_jac_A_elem<TElem, TFVGeom>);
	this->set_add_jac_M_elem_fct(	id, &T::template add_jac_M_elem<TElem, TFVGeom>);
	this->set_add_def_A_elem_fct(	id, &T::template add_def_A_elem<TElem, TFVGeom>);
	this->set_add_def_M_elem_fct(	id, &T::template add_def_M_elem<TElem, TFVGeom>);
	this->set_add_rhs_elem_fct(	 	id, &T::template add_rhs_elem<TElem, TFVGeom>);

	// error estimator parts
	this->set_prep_err_est_elem_loop(	id, &T::template prep_err_est_elem_loop<TElem, TFVGeom>);
	this->set_prep_err_est_elem(		id, &T::template prep_err_est_elem<TElem, TFVGeom>);
	this->set_compute_err_est_A_elem(	id, &T::template compute_err_est_A_elem<TElem, TFVGeom>);
	this->set_fsh_err_est_elem_loop(	id, &T::template fsh_err_est_elem_loop<TElem, TFVGeom>);
}

////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template class DiracSourceDisc<Domain1d>;
#endif
#ifdef UG_DIM_2
template class DiracSourceDisc<Domain2d>;
#endif
#ifdef UG_DIM_3
template class DiracSourceDisc<Domain3d>;
#endif


} // namespace ug


#endif /*__H__UG__LIB_DISC__SPACIAL_DISCRETIZATION__ELEM_DISC__NEUMANN_BOUNDARY__FV1__INNER_BOUNDARY_IMPL__*/
