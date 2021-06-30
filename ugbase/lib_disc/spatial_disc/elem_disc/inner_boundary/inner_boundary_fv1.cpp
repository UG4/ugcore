/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
 * Author: Markus Breit
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

#ifndef __H__UG__LIB_DISC__SPACIAL_DISCRETIZATION__ELEM_DISC__NEUMANN_BOUNDARY__FV1__INNER_BOUNDARY_IMPL__
#define __H__UG__LIB_DISC__SPACIAL_DISCRETIZATION__ELEM_DISC__NEUMANN_BOUNDARY__FV1__INNER_BOUNDARY_IMPL__

#include "inner_boundary.h"
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"
#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"



namespace ug
{

template<typename TDomain>
void FV1InnerBoundaryElemDisc<TDomain>::
set_flux_scale(number scale)
{
	set_flux_scale(make_sp(new ConstUserNumber<dim>(scale)));
}


template<typename TDomain>
void FV1InnerBoundaryElemDisc<TDomain>::
set_flux_scale(SmartPtr<CplUserData<number, dim> > scaleFct)
{
	m_fluxScale.set_data(scaleFct);
	this->register_import(m_fluxScale);
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void FV1InnerBoundaryElemDisc<TDomain>::
set_flux_scale(const char* luaScaleFctName)
{
	set_flux_scale(LuaUserDataFactory<number,dim>::create(luaScaleFctName));
}
#endif



template<typename TDomain>
void FV1InnerBoundaryElemDisc<TDomain>::
prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid)
{
	// check that Lagrange 1st order
	for(size_t i = 0; i < vLfeID.size(); ++i)
		if(vLfeID[i].type() != LFEID::LAGRANGE || vLfeID[i].order() != 1)
			UG_THROW("FV1InnerBoundary: 1st order lagrange expected.");

	// update assemble functions
	m_bNonRegularGrid = bNonRegularGrid;
	register_all_fv1_funcs();
}

template<typename TDomain>
bool FV1InnerBoundaryElemDisc<TDomain>::
use_hanging() const
{
	return true;
}



template<typename TDomain>
template<typename TElem, typename TFVGeom>
void FV1InnerBoundaryElemDisc<TDomain>::
prep_elem_loop(const ReferenceObjectID roid, const int si)
{
	m_si = si;

	//	set local positions
	if (!TFVGeom::usesHangingNodes)
	{
		static const int refDim = TElem::dim;
		static const TFVGeom& geo = GeomProvider<TFVGeom>::get();
		const MathVector<refDim>* vBFip = geo.bf_local_ips();
		const size_t numBFip = geo.num_bf_local_ips();

		m_fluxScale.template set_local_ips<refDim>(vBFip, numBFip, false);
	}
}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void FV1InnerBoundaryElemDisc<TDomain>::
fsh_elem_loop()
{}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void FV1InnerBoundaryElemDisc<TDomain>::
prep_elem(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[])
{
#ifdef UG_PARALLEL
	DistributedGridManager& dgm = *this->approx_space()->domain()->grid()->distributed_grid_manager();
	m_bCurrElemIsHSlave = dgm.get_status(elem) & ES_H_SLAVE;
#endif

	// on horizontal interfaces: only treat hmasters
	if (m_bCurrElemIsHSlave) return;

	// update Geometry for this element
	static TFVGeom& geo = GeomProvider<TFVGeom>::get();
	try {geo.update(elem, vCornerCoords, &(this->subset_handler()));}
	UG_CATCH_THROW("FV1InnerBoundaryElemDisc::prep_elem: "
						"Cannot update Finite Volume Geometry.");

	// set local positions
	if (TFVGeom::usesHangingNodes)
	{
		const int refDim = TElem::dim;
		const MathVector<refDim>* vBFip = geo.bf_local_ips();
		const size_t numBFip = geo.num_bf_local_ips();

		m_fluxScale.template set_local_ips<refDim>(vBFip, numBFip);
	}

	//	set global positions
	const MathVector<dim>* vBFip = geo.bf_global_ips();
	const size_t numBFip = geo.num_bf_global_ips();

	m_fluxScale.set_global_ips(vBFip, numBFip);
}

// assemble stiffness part of Jacobian
template<typename TDomain>
template<typename TElem, typename TFVGeom>
void FV1InnerBoundaryElemDisc<TDomain>::
add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	// on horizontal interfaces: only treat hmasters
	if (m_bCurrElemIsHSlave) return;

	// get finite volume geometry
	const static TFVGeom& fvgeom = GeomProvider<TFVGeom>::get();

	FluxDerivCond fdc;
	size_t nFct = u.num_fct();
	std::vector<LocalVector::value_type> uAtCorner(nFct);
	for (size_t i = 0; i < fvgeom.num_bf(); ++i)
	{
		// get current BF
		const typename TFVGeom::BF& bf = fvgeom.bf(i);

		// get associated node and subset index
		const int co = bf.node_id();

		// get solution at the corner of the bf
		for (size_t fct = 0; fct < nFct; fct++)
			uAtCorner[fct] = u(fct,co);

		// get corner coordinates
		const MathVector<dim>& cc = bf.global_corner(0);

		if (!fluxDensityDerivFct(uAtCorner, elem, cc, m_si, fdc))
			UG_THROW("FV1InnerBoundaryElemDisc::add_jac_A_elem:"
							" Call to fluxDensityDerivFct resulted did not succeed.");
		
		// scale with volume of BF and flux scale (if applicable)
		number scale = bf.volume();
		if (m_fluxScale.data_given())
			scale *= m_fluxScale[i];

		for (size_t j=0; j<fdc.fluxDeriv.size(); j++)
			for (size_t k=0; k<fdc.fluxDeriv[j].size(); k++)
				fdc.fluxDeriv[j][k].second *= scale;
		
		// add to Jacobian
		for (size_t j=0; j<fdc.fluxDeriv.size(); j++)
		{
			for (size_t k=0; k<fdc.fluxDeriv[j].size(); k++)
			{
				if (fdc.from[j] != InnerBoundaryConstants::_IGNORE_)
					J(fdc.from[j], co, fdc.fluxDeriv[j][k].first, co) += fdc.fluxDeriv[j][k].second;
				if (fdc.to[j] != InnerBoundaryConstants::_IGNORE_)
					J(fdc.to[j], co, fdc.fluxDeriv[j][k].first, co)	-= fdc.fluxDeriv[j][k].second;
			}
		}	
	}
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void FV1InnerBoundaryElemDisc<TDomain>::
add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void FV1InnerBoundaryElemDisc<TDomain>::
add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	// on horizontal interfaces: only treat hmasters
	if (m_bCurrElemIsHSlave) return;

	// get finite volume geometry
	static TFVGeom& fvgeom = GeomProvider<TFVGeom>::get();

	FluxCond fc;
	size_t nFct = u.num_fct();
	std::vector<LocalVector::value_type> uAtCorner(nFct);

	// loop Boundary Faces
	for (size_t i = 0; i < fvgeom.num_bf(); ++i)
	{
		// get current BF
		const typename TFVGeom::BF& bf = fvgeom.bf(i);
		
		// get associated node and subset index
		const int co = bf.node_id();

		// get solution at the corner of the bf
		for (size_t fct = 0; fct < nFct; fct++)
			uAtCorner[fct] = u(fct,co);

		// get corner coordinates
		const MathVector<dim>& cc = bf.global_corner(0);

		// get flux densities in that node
		if (!fluxDensityFct(uAtCorner, elem, cc, m_si, fc))
		{
			UG_THROW("FV1InnerBoundaryElemDisc::add_def_A_elem:"
						" Call to fluxDensityFct did not succeed.");
		}

		// scale with volume of BF and flux scale (if applicable)
		number scale = bf.volume();
		if (m_fluxScale.data_given())
			scale *= m_fluxScale[i];

		for (size_t j=0; j<fc.flux.size(); j++)
			fc.flux[j] *= scale;

		// add to defect
		for (size_t j=0; j<fc.flux.size(); j++)
		{
			if (fc.from[j] != InnerBoundaryConstants::_IGNORE_) d(fc.from[j], co) += fc.flux[j];
			if (fc.to[j] != InnerBoundaryConstants::_IGNORE_) d(fc.to[j], co) -= fc.flux[j];
		}
	}
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void FV1InnerBoundaryElemDisc<TDomain>::
add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void FV1InnerBoundaryElemDisc<TDomain>::
add_rhs_elem(LocalVector& rhs, GridObject* elem, const MathVector<dim> vCornerCoords[])
{}


// ////////////////////////////////////////////////////////////////////////////////////////////////
//   error estimation (begin)   																 //

//	prepares the loop over all elements of one type for the computation of the error estimator
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void FV1InnerBoundaryElemDisc<TDomain>::
prep_err_est_elem_loop(const ReferenceObjectID roid, const int si)
{
	m_si = si;

	//	get the error estimator data object and check that it is of the right type
	//	we check this at this point in order to be able to dispense with this check later on
	//	(i.e. in prep_err_est_elem and compute_err_est_A_elem())
	if (this->m_spErrEstData.get() == NULL)
	{
		UG_THROW("No ErrEstData object has been given to this ElemDisc!");
	}

	err_est_type* err_est_data = dynamic_cast<err_est_type*>(this->m_spErrEstData.get());

	if (!err_est_data)
	{
		UG_THROW("Dynamic cast to MultipleSideAndElemErrEstData failed."
				<< std::endl << "Make sure you handed the correct type of ErrEstData to this discretization.");
	}

	if (!err_est_data->equal_side_order())
	{
		UG_THROW("The underlying error estimator data objects of this discretization's "
				 "error estimator do not all have the same integration orders. This case "
				 "is not supported by the implementation. If you need it, implement!");
	}

	if (err_est_data->num() < 1)
	{
		UG_THROW("No underlying error estimator data objects present. No IPs can be determined.");
	}

	//	set local positions
	if (!TFVGeom::usesHangingNodes)
	{
		static const int refDim = TElem::dim;

		// get local IPs
		size_t numSideIPs;
		const MathVector<refDim>* sideIPs;
		try
		{
			numSideIPs = err_est_data->get(0)->num_side_ips(roid);
			sideIPs = err_est_data->get(0)->template side_local_ips<refDim>(roid);
		}
		UG_CATCH_THROW("Integration points for error estimator cannot be set.");

		// set local IPs in import
		m_fluxScale.template set_local_ips<refDim>(sideIPs, numSideIPs, false);

		// store values of shape functions in local IPs
		LagrangeP1<typename reference_element_traits<TElem>::reference_element_type> trialSpace
					= Provider<LagrangeP1<typename reference_element_traits<TElem>::reference_element_type> >::get();

		m_shapeValues.resize(numSideIPs, trialSpace.num_sh());
		for (size_t ip = 0; ip < numSideIPs; ip++)
			trialSpace.shapes(m_shapeValues.shapesAtSideIP(ip), sideIPs[ip]);
	}
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void FV1InnerBoundaryElemDisc<TDomain>::
prep_err_est_elem(const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
#ifdef UG_PARALLEL
	DistributedGridManager& dgm = *this->approx_space()->domain()->grid()->distributed_grid_manager();
	m_bCurrElemIsHSlave = dgm.get_status(elem) & ES_H_SLAVE;
#endif

	// on horizontal interfaces: only treat hmasters
	if (m_bCurrElemIsHSlave) return;

	// get error estimator
	err_est_type* err_est_data = dynamic_cast<err_est_type*>(this->m_spErrEstData.get());

	// roid
	ReferenceObjectID roid = elem->reference_object_id();

	// set local positions
	if (TFVGeom::usesHangingNodes)
	{
		static const int refDim = TElem::dim;

		size_t numSideIPs;
		const MathVector<refDim>* sideIPs;
		try
		{
			numSideIPs = err_est_data->get(0)->num_side_ips(roid);
			sideIPs = err_est_data->get(0)->template side_local_ips<refDim>(roid);
		}
		UG_CATCH_THROW("Integration points for error estimator cannot be set.");

		// set local IPs in import
		m_fluxScale.template set_local_ips<refDim>(sideIPs, numSideIPs);

		// store values of shape functions in local IPs
		LagrangeP1<typename reference_element_traits<TElem>::reference_element_type> trialSpace
					= Provider<LagrangeP1<typename reference_element_traits<TElem>::reference_element_type> >::get();

		m_shapeValues.resize(numSideIPs, trialSpace.num_sh());
		for (size_t ip = 0; ip < numSideIPs; ip++)
			trialSpace.shapes(m_shapeValues.shapesAtSideIP(ip), sideIPs[ip]);
	}

	// set global positions
	size_t numSideIPs;
	MathVector<dim>* sideIPs;

	try
	{
		numSideIPs = err_est_data->get(0)->num_side_ips(roid);
		sideIPs = err_est_data->get(0)->side_global_ips(elem, vCornerCoords);
	}
	UG_CATCH_THROW("Global integration points for error estimator cannot be set.");

	// set local IPs in import
	m_fluxScale.set_global_ips(&sideIPs[0], numSideIPs);
}

//	computes the error estimator contribution for one element
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void FV1InnerBoundaryElemDisc<TDomain>::
compute_err_est_A_elem(const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[], const number& scale)
{
	// on horizontal interfaces: only treat hmasters
	if (m_bCurrElemIsHSlave) return;

	// get error estimator
	err_est_type* err_est_data = dynamic_cast<err_est_type*>(this->m_spErrEstData.get());

	// cast this elem to side_type of error estimator
	typename SideAndElemErrEstData<TDomain>::side_type* side =
		dynamic_cast<typename SideAndElemErrEstData<TDomain>::side_type*>(elem);
	if (!side)
	{
		UG_THROW("Error in DependentNeumannBoundaryFV1<TDomain>::compute_err_est_A_elem():\n"
				 "Element that error assembling routine is called for has the wrong type.");
	}

	// global IPs
	ReferenceObjectID roid = side->reference_object_id();
	size_t numSideIPs = err_est_data->get(0)->num_side_ips(roid);
	MathVector<dim>* globIPs = err_est_data->get(0)->side_global_ips(side, vCornerCoords);

	// loop IPs
	try
	{
		for (size_t sip = 0; sip < numSideIPs; sip++)
		{
			// get values of u at ip (interpolate)
			size_t nFct = u.num_fct();
			std::vector<LocalVector::value_type> uAtIP = std::vector<LocalVector::value_type>(nFct);

			for (size_t fct = 0; fct < nFct; fct++)
			{
				uAtIP[fct] = 0.0;
				for (size_t sh = 0; sh < m_shapeValues.num_sh(); sh++)
					uAtIP[fct] += u(fct,sh) * m_shapeValues.shapeAtSideIP(sh,sip);
			}

			// ip coordinates
			const MathVector<dim>& ipCoords = globIPs[sip];

			FluxCond fc;
			if (!fluxDensityFct(uAtIP, elem, ipCoords, m_si, fc))
			{
				UG_THROW("FV1InnerBoundaryElemDisc::compute_err_est_A_elem:"
							" Call to fluxDensityFct did not succeed.");
			}

			// scale fluxes if applicable
			if (m_fluxScale.data_given())
				for (size_t j = 0; j < fc.flux.size(); ++j)
					fc.flux[j] *= m_fluxScale[sip];


			// subtract from estimator values
			// sign must be opposite of that in add_def_A_elem
			// as the difference between this and the actual flux of the unknown is calculated
			for (size_t j=0; j<fc.flux.size(); j++)
			{
				if (fc.from[j] != InnerBoundaryConstants::_IGNORE_)
					(*err_est_data->get(this->m_fctGrp[fc.from[j]])) (side,sip) -= scale * fc.flux[j];
				if (fc.to[j] != InnerBoundaryConstants::_IGNORE_)
					(*err_est_data->get(this->m_fctGrp[fc.to[j]])) (side,sip) += scale * fc.flux[j];
			}
		}
	}
	UG_CATCH_THROW("Values for the error estimator could not be assembled at every IP." << std::endl
			<< "Maybe wrong type of ErrEstData object? This implementation needs: MultipleSideAndElemErrEstData.");
}

//	postprocesses the loop over all elements of one type in the computation of the error estimator
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void FV1InnerBoundaryElemDisc<TDomain>::
fsh_err_est_elem_loop()
{
	// nothing to do
}

//   error estimation (end)   																	 //
// ////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////

// register for 1D
template<typename TDomain>
void FV1InnerBoundaryElemDisc<TDomain>::
register_all_fv1_funcs()
{
//	get all grid element types in the dimension below
	typedef typename domain_traits<dim>::ManifoldElemList ElemList;

//	switch assemble functions
	boost::mpl::for_each<ElemList>(RegisterFV1(this));
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void FV1InnerBoundaryElemDisc<TDomain>::
register_fv1_func()
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
template class FV1InnerBoundaryElemDisc<Domain1d>;
#endif
#ifdef UG_DIM_2
template class FV1InnerBoundaryElemDisc<Domain2d>;
#endif
#ifdef UG_DIM_3
template class FV1InnerBoundaryElemDisc<Domain3d>;
#endif


} // namespace ug


#endif /*__H__UG__LIB_DISC__SPACIAL_DISCRETIZATION__ELEM_DISC__NEUMANN_BOUNDARY__FV1__INNER_BOUNDARY_IMPL__*/
