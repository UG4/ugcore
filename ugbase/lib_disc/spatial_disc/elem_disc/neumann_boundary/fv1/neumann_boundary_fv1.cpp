/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
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

#include "neumann_boundary_fv1.h"
#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"

namespace ug {


////////////////////////////////////////////////////////////////////////////////
//	Constructor
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
NeumannBoundaryFV1<TDomain>::NeumannBoundaryFV1(const char* function)
 :NeumannBoundaryBase<TDomain>(function)
{
	register_all_funcs(false);
}

template<typename TDomain>
void NeumannBoundaryFV1<TDomain>::
prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid)
{
	UG_COND_THROW(bNonRegularGrid && (TDomain::dim == 3),
				  "NeumannBoundaryFV1: Hanging Nodes not implemented.");

	if(vLfeID.size() != 1)
		UG_THROW("NeumannBoundary: Need exactly 1 function.");

	if(vLfeID[0].order() != 1 || vLfeID[0].type() != LFEID::LAGRANGE)
		UG_THROW("NeumannBoundary: FV Scheme only implemented for 1st order Lagrange.");
}

template<typename TDomain>
void NeumannBoundaryFV1<TDomain>::
add(SmartPtr<CplUserData<number, dim> > data, const char* BndSubsets, const char* InnerSubsets)
{
	m_vNumberData.push_back(NumberData(data, BndSubsets, InnerSubsets));
	this->add_inner_subsets(InnerSubsets);
}

template<typename TDomain>
void NeumannBoundaryFV1<TDomain>::
add(SmartPtr<CplUserData<number, dim, bool> > user, const char* BndSubsets, const char* InnerSubsets)
{
	m_vBNDNumberData.push_back(BNDNumberData(user, BndSubsets, InnerSubsets));
	this->add_inner_subsets(InnerSubsets);
}

template<typename TDomain>
void NeumannBoundaryFV1<TDomain>::
add(SmartPtr<CplUserData<MathVector<dim>, dim> > user, const char* BndSubsets, const char* InnerSubsets)
{
	m_vVectorData.push_back(VectorData(user, BndSubsets, InnerSubsets));
	this->add_inner_subsets(InnerSubsets);
}

template<typename TDomain>
void NeumannBoundaryFV1<TDomain>::update_subset_groups()
{
	for(size_t i = 0; i < m_vNumberData.size(); ++i)
		update_subset_groups(m_vNumberData[i]);
	for(size_t i = 0; i < m_vBNDNumberData.size(); ++i)
		update_subset_groups(m_vBNDNumberData[i]);
	for(size_t i = 0; i < m_vVectorData.size(); ++i)
		update_subset_groups(m_vVectorData[i]);
}

////////////////////////////////////////////////////////////////////////////////
//	assembling functions
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NeumannBoundaryFV1<TDomain>::
prep_elem_loop(ReferenceObjectID_t roid, const int si)
{
	update_subset_groups();
	m_si = si;

//	register subsetIndex at Geometry
	static TFVGeom& geo = GeomProvider<TFVGeom >::get();

//	request subset indices as boundary subset. This will force the
//	creation of boundary subsets when calling geo.update

	for(size_t i = 0; i < m_vNumberData.size(); ++i){
		if(!m_vNumberData[i].InnerSSGrp.contains(m_si)) continue;
		for(size_t s = 0; s < m_vNumberData[i].BndSSGrp.size(); ++s){
			const int si = m_vNumberData[i].BndSSGrp[s];
			geo.add_boundary_subset(si);
		}
	}
	for(size_t i = 0; i < m_vBNDNumberData.size(); ++i){
		if(!m_vBNDNumberData[i].InnerSSGrp.contains(m_si)) continue;
		for(size_t s = 0; s < m_vBNDNumberData[i].BndSSGrp.size(); ++s){
			const int si = m_vBNDNumberData[i].BndSSGrp[s];
			geo.add_boundary_subset(si);
		}
	}
	for(size_t i = 0; i < m_vVectorData.size(); ++i){
		if(!m_vVectorData[i].InnerSSGrp.contains(m_si)) continue;
		for(size_t s = 0; s < m_vVectorData[i].BndSSGrp.size(); ++s){
			const int si = m_vVectorData[i].BndSSGrp[s];
			geo.add_boundary_subset(si);
		}
	}

//	clear imports, since we will set them afterwards
	this->clear_imports();

	ReferenceObjectID_t id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;

//	set lin defect fct for imports
	for(size_t data = 0; data < m_vNumberData.size(); ++data)
	{
		if(!m_vNumberData[data].InnerSSGrp.contains(m_si)) continue;
		m_vNumberData[data].import.set_fct(id,
		                                   &m_vNumberData[data],
		                                   &NumberData::template lin_def<TElem, TFVGeom>);

		this->register_import(m_vNumberData[data].import);
		m_vNumberData[data].import.set_rhs_part();
	}
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NeumannBoundaryFV1<TDomain>::
prep_elem(const LocalVector& u, GridObject* elem, ReferenceObjectID_t roid, const MathVector<dim> vCornerCoords[])
{
//  update Geometry for this element
	static TFVGeom& geo = GeomProvider<TFVGeom >::get();
	try{
		geo.update(elem, vCornerCoords, &(this->subset_handler()));
	}
	UG_CATCH_THROW("NeumannBoundaryFV1::prep_elem: "
						"Cannot update Finite Volume Geometry.");

	for(size_t i = 0; i < m_vNumberData.size(); ++i)
		if(m_vNumberData[i].InnerSSGrp.contains(m_si))
			m_vNumberData[i].template extract_bip<TElem, TFVGeom>(geo);
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NeumannBoundaryFV1<TDomain>::
add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	static const TFVGeom& geo = GeomProvider<TFVGeom >::get();
	using BF = typename TFVGeom::BF;

//	Number Data
	for(size_t data = 0; data < m_vNumberData.size(); ++data){
		if(!m_vNumberData[data].InnerSSGrp.contains(m_si)) continue;
		size_t ip = 0;
		for(size_t s = 0; s < m_vNumberData[data].BndSSGrp.size(); ++s){
			const int si = m_vNumberData[data].BndSSGrp[s];
			const std::vector<BF>& vBF = geo.bf(si);

			for(size_t i = 0; i < vBF.size(); ++i, ++ip){
				const int co = vBF[i].node_id();
				d(_C_, co) -= m_vNumberData[data].import[ip]
				                                    * vBF[i].volume();
			}
		}
	}

//	conditional Number Data
	for(size_t data = 0; data < m_vBNDNumberData.size(); ++data){
		if(!m_vBNDNumberData[data].InnerSSGrp.contains(m_si)) continue;
		for(size_t s = 0; s < m_vBNDNumberData[data].BndSSGrp.size(); ++s)	{
			const int si = m_vBNDNumberData[data].BndSSGrp[s];
			const std::vector<BF>& vBF = geo.bf(si);

			for(size_t i = 0; i < vBF.size(); ++i){
				number val = 0.0;
				if(!(*m_vBNDNumberData[data].functor)(val, vBF[i].global_ip(), this->time(), si))
					continue;

				const int co = vBF[i].node_id();
				d(_C_, co) -= val * vBF[i].volume();
			}
		}
	}

// 	vector data
	for(size_t data = 0; data < m_vVectorData.size(); ++data){
		if(!m_vVectorData[data].InnerSSGrp.contains(m_si)) continue;
		for(size_t s = 0; s < m_vVectorData[data].BndSSGrp.size(); ++s){
			const int si = m_vVectorData[data].BndSSGrp[s];
			const std::vector<BF>& vBF = geo.bf(si);

			for(size_t i = 0; i < vBF.size(); ++i){
				MathVector<dim> val;
				(*m_vVectorData[data].functor)(val, vBF[i].global_ip(), this->time(), si);

				const int co = vBF[i].node_id();
				d(_C_, co) -= VecDot(val, vBF[i].normal());
			}
		}
	}
}

template<typename TDomain>
template<typename TElem, typename TGeom>
void NeumannBoundaryFV1<TDomain>::
fsh_elem_loop()
{
//	remove subsetIndex from Geometry
	static TGeom& geo = GeomProvider<TGeom >::get();


//	unrequest subset indices as boundary subset. This will force the
//	creation of boundary subsets when calling geo.update

	for(size_t i = 0; i < m_vNumberData.size(); ++i){
		if(!m_vNumberData[i].InnerSSGrp.contains(m_si)) continue;
		for(size_t s = 0; s < m_vNumberData[i].BndSSGrp.size(); ++s){
			const int si = m_vNumberData[i].BndSSGrp[s];
			geo.remove_boundary_subset(si);
			geo.reset_curr_elem();
		}
	}

	for(size_t i = 0; i < m_vBNDNumberData.size(); ++i){
		if(!m_vBNDNumberData[i].InnerSSGrp.contains(m_si)) continue;
		for(size_t s = 0; s < m_vBNDNumberData[i].BndSSGrp.size(); ++s){
			const int si = m_vBNDNumberData[i].BndSSGrp[s];
			geo.remove_boundary_subset(si);
			geo.reset_curr_elem();
		}
	}

	for(size_t i = 0; i < m_vVectorData.size(); ++i){
		if(!m_vVectorData[i].InnerSSGrp.contains(m_si)) continue;
		for(size_t s = 0; s < m_vVectorData[i].BndSSGrp.size(); ++s){
			const int si = m_vVectorData[i].BndSSGrp[s];
			geo.remove_boundary_subset(si);
			geo.reset_curr_elem();
		}
	}
}


///////////////////////////////////////////////////////////////////////////////////////////////////
//   error estimation (begin)  																	 //

//	prepares the loop over all elements of one type for the computation of the error estimator
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void NeumannBoundaryFV1<TDomain>::
prep_err_est_elem_loop(ReferenceObjectID_t roid, const int si)
{
	//	get the error estimator data object and check that it is of the right type
	//	we check this at this point in order to be able to dispense with this check later on
	//	(i.e. in prep_err_est_elem and compute_err_est_A_elem())
	if (this->m_spErrEstData.get() == nullptr)
	{
		UG_THROW("No ErrEstData object has been given to this ElemDisc!");
	}

	err_est_type* err_est_data = dynamic_cast<err_est_type*>(this->m_spErrEstData.get());

	if (!err_est_data)
	{
		UG_THROW("Dynamic cast to SideAndElemErrEstData failed."
				<< std::endl << "Make sure you handed the correct type of ErrEstData to this discretization.");
	}

	update_subset_groups();
	m_si = si;

//	clear imports, since we will set them now
	this->clear_imports();

	ReferenceObjectID_t id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;

//	set lin defect fct for imports
	for (size_t data = 0; data < m_vNumberData.size(); ++data)
	{
		if (!m_vNumberData[data].InnerSSGrp.contains(m_si)) continue;
		m_vNumberData[data].import.set_fct(id,
										   &m_vNumberData[data],
										   &NumberData::template lin_def<TElem, TFVGeom>);

		this->register_import(m_vNumberData[data].import);
		m_vNumberData[data].import.set_rhs_part();
	}

// get local IPs for imports
	static constexpr int refDim = TElem::dim;

	size_t numSideIPs;
	const MathVector<refDim>* sideIPs;
	try
	{
		numSideIPs = err_est_data->num_all_side_ips(roid);
		sideIPs = err_est_data->template side_local_ips<refDim>(roid);

		if (!sideIPs) return;	// is nullptr if TElem is not of the same dim as TDomain
	}
	UG_CATCH_THROW("Integration points for error estimator cannot be set.");

	for (size_t i = 0; i < m_vNumberData.size(); ++i)
	{
		if (m_vNumberData[i].InnerSSGrp.contains(m_si))
			m_vNumberData[i].template set_local_ips<refDim>(sideIPs, numSideIPs);
	}
};

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NeumannBoundaryFV1<TDomain>::
prep_err_est_elem(const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
// get global IPs for imports
	size_t numSideIPs;
	MathVector<dim>* sideIPs;

	try
	{
		err_est_type* err_est_data = dynamic_cast<err_est_type*>(this->m_spErrEstData.get());

		numSideIPs = err_est_data->num_all_side_ips(elem->reference_object_id());
		sideIPs = err_est_data->all_side_global_ips(elem, vCornerCoords);
	}
	UG_CATCH_THROW("Global integration points for error estimator cannot be set.");


	for (size_t i = 0; i < m_vNumberData.size(); ++i)
	{
		if (m_vNumberData[i].InnerSSGrp.contains(m_si))
			m_vNumberData[i].set_global_ips(&sideIPs[0], numSideIPs);
	}
}

//	computes the error estimator contribution for one element
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void NeumannBoundaryFV1<TDomain>::
compute_err_est_rhs_elem(GridObject* elem, const MathVector<dim> vCornerCoords[], const number& scale)
{
	using ref_elem_type = typename reference_element_traits<TElem>::reference_element_type;

	err_est_type* err_est_data = dynamic_cast<err_est_type*>(this->m_spErrEstData.get());

	if (err_est_data->surface_view().get() == nullptr) {UG_THROW("Error estimator has nullptr surface view.");}
	auto pErrEstGrid = (MultiGrid*) (err_est_data->surface_view()->subset_handler()->multi_grid());

//	get the sides of the element
	typename MultiGrid::traits<typename SideAndElemErrEstData<TDomain>::side_type>::secure_container side_list;
	pErrEstGrid->associated_elements_sorted(side_list, (TElem*) elem);
	if (side_list.size() != static_cast<size_t>(ref_elem_type::numSides))
		UG_THROW ("Mismatch of numbers of sides in 'ConvectionDiffusionFV1::compute_err_est_elem'");

	//	Number Data
	for (size_t data = 0; data < m_vNumberData.size(); ++data)
	{
		if (!m_vNumberData[data].InnerSSGrp.contains(m_si)) continue;

		// loop sides
		for (size_t side = 0; side < side_list.size(); side++)
		{
			for (size_t ss = 0; ss < m_vNumberData[data].BndSSGrp.size(); ss++)
			{
				// is the bnd cond subset index the same as the current side's?
				if (err_est_data->surface_view()->subset_handler()->get_subset_index(side_list[side])
					!= m_vNumberData[data].BndSSGrp[ss]) continue;

				const size_t ipIndexStart = err_est_data->first_side_ips(elem->reference_object_id(), side);

				for (size_t sip = 0; sip < err_est_data->num_side_ips(side_list[side]); sip++)
				{
					size_t ip = ipIndexStart + sip;
					// sign must be negative, since the data represents efflux
					(*err_est_data)(side_list[side],sip) -= scale * m_vNumberData[data].import[ip];
				}
			}
		}
	}

	// store global side IPs
	ReferenceObjectID_t roid = elem->reference_object_id();
	MathVector<dim>* glob_ips = err_est_data->all_side_global_ips(elem, vCornerCoords);

	//	conditional Number Data
	for (size_t data = 0; data < m_vBNDNumberData.size(); ++data)
	{
		if (!m_vBNDNumberData[data].InnerSSGrp.contains(m_si)) continue;

		// loop sides
		for (size_t side = 0; side < side_list.size(); side++)
		{
			for (size_t ss = 0; ss < m_vBNDNumberData[data].BndSSGrp.size(); ss++)
			{
				// is the bnd cond subset index the same as the current side's?
				if (err_est_data->surface_view()->subset_handler()->get_subset_index(side_list[side])
					!= m_vBNDNumberData[data].BndSSGrp[ss]) continue;

				int si = m_vBNDNumberData[data].BndSSGrp[ss];

				const size_t ipIndexStart = err_est_data->first_side_ips(roid, side);

				for (size_t sip = 0; sip < err_est_data->num_side_ips(side_list[side]); sip++)
				{
					size_t ip = ipIndexStart + sip;
					number val = 0.0;
					if (!(*m_vBNDNumberData[data].functor)(val, glob_ips[ip], this->time(), si)) continue;

					// sign must be negative, since the data represents efflux
					(*err_est_data)(side_list[side],sip) -= scale * val;
				}
			}
		}
	}

	// 	vector data
	for (size_t data = 0; data < m_vVectorData.size(); ++data)
	{
		if (!m_vVectorData[data].InnerSSGrp.contains(m_si)) continue;

		// loop sides
		for (size_t side = 0; side < side_list.size(); side++)
		{
			MathVector<dim> fd, normal;

			// normal on side
			SideNormal<ref_elem_type,dim>(normal, side, vCornerCoords);
			VecNormalize(normal, normal);

			for (size_t ss = 0; ss < m_vVectorData[data].BndSSGrp.size(); ss++)
			{
				// is the bnd cond subset index the same as the current side's?
				if (err_est_data->surface_view()->subset_handler()->get_subset_index(side_list[side])
					!= m_vVectorData[data].BndSSGrp[ss]) continue;

				int si = m_vBNDNumberData[data].BndSSGrp[ss];

				const size_t ipIndexStart = err_est_data->first_side_ips(roid, side);

				for (size_t sip = 0; sip < err_est_data->num_side_ips(side_list[side]); sip++)
				{
					size_t ip = ipIndexStart + sip;
					(*m_vVectorData[data].functor)(fd, glob_ips[ip], this->time(), si);

					// sign must be negative, since the data represents efflux
					(*err_est_data)(side_list[side],sip) -= scale * VecDot(fd, normal);
				}
			}
		}
	}
}

//	postprocesses the loop over all elements of one type in the computation of the error estimator
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void NeumannBoundaryFV1<TDomain>::
fsh_err_est_elem_loop()
{
	// finish the element loop in the same way as the actual discretization
	this->template fsh_elem_loop<TElem, TFVGeom> ();
}

//   error estimation (end)     																 //
///////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// Number Data
////////////////////////////////////////////////////////////////////////////////


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NeumannBoundaryFV1<TDomain>::NumberData::
lin_def(const LocalVector& u,
            std::vector<std::vector<number> > vvvLinDef[],
            const size_t nip)
{
//  get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();
	using BF = typename TFVGeom::BF;

	for(size_t s = 0; s < this->BndSSGrp.size(); ++s)
	{
		const int si = this->BndSSGrp[s];
		const std::vector<BF>& vBF = geo.bf(si);
		for(size_t i = 0; i < vBF.size(); ++i){
			const int co = vBF[i].node_id();
			vvvLinDef[i][_C_][co] -= vBF[i].volume();
		}
	}
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NeumannBoundaryFV1<TDomain>::NumberData::
extract_bip(const TFVGeom& geo)
{
	using BF = typename TFVGeom::BF;
	static constexpr int locDim = TElem::dim;

	std::vector<MathVector<locDim> >* vLocIP = local_ips<locDim>();

	vLocIP->clear();
	vGloIP.clear();
	for(size_t s = 0; s < this->BndSSGrp.size(); s++)
	{
		const int si = this->BndSSGrp[s];
		const std::vector<BF>& vBF = geo.bf(si);
		for(size_t i = 0; i < vBF.size(); ++i)
		{
			const BF& bf = vBF[i];
			vLocIP->push_back(bf.local_ip());
			vGloIP.push_back(bf.global_ip());
		}
	}

//	either both are empty or none is empty
	UG_ASSERT((vLocIP->empty() && vGloIP.empty()) || (!(vLocIP->empty() || vGloIP.empty())),
			  "Either vLocIP and vGloIP both have to be empty or both have to be filled!");

	if(vLocIP->empty()){
		import.template set_local_ips<locDim>(nullptr, 0);
		import.set_global_ips(nullptr, 0);
	}
	else{
		import.template set_local_ips<locDim>(&(*vLocIP)[0], vLocIP->size());
		import.set_global_ips(&vGloIP[0], vGloIP.size());
	}
}

template<typename TDomain>
template <int refDim>
std::vector<MathVector<refDim> >* NeumannBoundaryFV1<TDomain>::NumberData::local_ips()
{
	switch (refDim)
	{
		case 1: return (std::vector<MathVector<refDim> >*)(&vLocIP_dim1);
		case 2: return (std::vector<MathVector<refDim> >*)(&vLocIP_dim2);
		case 3: return (std::vector<MathVector<refDim> >*)(&vLocIP_dim3);
	}
}

////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template<>
void NeumannBoundaryFV1<Domain1d>::register_all_funcs(bool bHang)
{
	register_func<RegularEdge, FV1Geometry<RegularEdge, dim> >();
}
#endif

#ifdef UG_DIM_2
template<>
void NeumannBoundaryFV1<Domain2d>::register_all_funcs(bool bHang)
{
	register_func<RegularEdge, FV1Geometry<RegularEdge, dim> >();
	register_func<Triangle, FV1Geometry<Triangle, dim> >();
	register_func<Quadrilateral, FV1Geometry<Quadrilateral, dim> >();
}
#endif

#ifdef UG_DIM_3
template<>
void NeumannBoundaryFV1<Domain3d>::register_all_funcs(bool bHang)
{
	if(!bHang){
		register_func<RegularEdge, FV1Geometry<RegularEdge, dim> >();
		register_func<Triangle, FV1Geometry<Triangle, dim> >();
		register_func<Quadrilateral, FV1Geometry<Quadrilateral, dim> >();
		register_func<Tetrahedron, FV1Geometry<Tetrahedron, dim> >();
		register_func<Prism, FV1Geometry<Prism, dim> >();
		register_func<Pyramid, FV1Geometry<Pyramid, dim> >();
		register_func<Hexahedron, FV1Geometry<Hexahedron, dim> >();
		register_func<Octahedron, FV1Geometry<Octahedron, dim> >();
	}
	else UG_THROW("NeumannBoundaryFV1: Hanging Nodes not implemented.");
}
#endif

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NeumannBoundaryFV1<TDomain>::register_func()
{
	ReferenceObjectID_t id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	using T = this_type;

	this->clear_add_fct(id);

	this->set_prep_elem_loop_fct(id, &T::template prep_elem_loop<TElem, TFVGeom>);
	this->set_prep_elem_fct(	 id, &T::template prep_elem<TElem, TFVGeom>);
	this->set_add_rhs_elem_fct(	 id, &T::template add_rhs_elem<TElem, TFVGeom>);
	this->set_fsh_elem_loop_fct( id, &T::template fsh_elem_loop<TElem, TFVGeom>);

	this->set_add_jac_A_elem_fct(	 id, &T::template add_jac_A_elem<TElem, TFVGeom>);
	this->set_add_jac_M_elem_fct(	 id, &T::template add_jac_M_elem<TElem, TFVGeom>);
	this->set_add_def_A_elem_fct(	 id, &T::template add_def_A_elem<TElem, TFVGeom>);
	this->set_add_def_M_elem_fct(	 id, &T::template add_def_M_elem<TElem, TFVGeom>);

	// error estimator parts
	this->set_prep_err_est_elem_loop(id, &T::template prep_err_est_elem_loop<TElem, TFVGeom>);
	this->set_prep_err_est_elem(id, &T::template prep_err_est_elem<TElem, TFVGeom>);
	this->set_compute_err_est_rhs_elem(id, &T::template compute_err_est_rhs_elem<TElem, TFVGeom>);
	this->set_fsh_err_est_elem_loop(id, &T::template fsh_err_est_elem_loop<TElem, TFVGeom>);
}

////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template class NeumannBoundaryFV1<Domain1d>;
#endif
#ifdef UG_DIM_2
template class NeumannBoundaryFV1<Domain2d>;
#endif
#ifdef UG_DIM_3
template class NeumannBoundaryFV1<Domain3d>;
#endif

} // namespace ug

