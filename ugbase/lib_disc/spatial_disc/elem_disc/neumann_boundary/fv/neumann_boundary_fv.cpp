/*
 * neumann_boundary_fv.cpp
 *
 *  Created on: 14.10.2010
 *      Author: andreasvogel
 */

#include "neumann_boundary_fv.h"
#include "lib_disc/spatial_disc/disc_util/fvho_geom.h"
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
//	Constructor
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
NeumannBoundaryFV<TDomain>::NeumannBoundaryFV(const char* function)
 :NeumannBoundaryBase<TDomain>(function),
  m_order(1), m_lfeID(LFEID::LAGRANGE, m_order)
{
	register_all_funcs(m_order);
}

template<typename TDomain>
void NeumannBoundaryFV<TDomain>::
prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid)
{
	if(bNonRegularGrid)
		UG_THROW("NeumannBoundary: Hanging Nodes not implemented.");

//	check number
	if(vLfeID.size() != 1)
		UG_THROW("NeumannBoundaryFV: Need exactly 1 function.");

	if(vLfeID[0].type() != LFEID::LAGRANGE)
		UG_THROW("NeumannBoundary: FV Scheme only implemented for 1st order.");

//	check that not ADAPTIVE
	if(vLfeID[0].order() < 1)
		UG_THROW("NeumannBoundary: Adaptive order not implemented.");

//	set order
	m_lfeID = vLfeID[0];
	m_order = vLfeID[0].order();

	register_all_funcs(m_order);
}

template<typename TDomain>
void NeumannBoundaryFV<TDomain>::
add(SmartPtr<CplUserData<number, dim> > data, const char* BndSubsets, const char* InnerSubsets)
{
	m_vNumberData.push_back(NumberData(data, BndSubsets, InnerSubsets, this));
	this->add_inner_subsets(InnerSubsets);
}

template<typename TDomain>
void NeumannBoundaryFV<TDomain>::
add(SmartPtr<CplUserData<number, dim, bool> > user, const char* BndSubsets, const char* InnerSubsets)
{
	m_vBNDNumberData.push_back(BNDNumberData(user, BndSubsets, InnerSubsets));
	this->add_inner_subsets(InnerSubsets);
}

template<typename TDomain>
void NeumannBoundaryFV<TDomain>::
add(SmartPtr<CplUserData<MathVector<dim>, dim> > user, const char* BndSubsets, const char* InnerSubsets)
{
	m_vVectorData.push_back(VectorData(user, BndSubsets, InnerSubsets));
	this->add_inner_subsets(InnerSubsets);
}

template<typename TDomain>
void NeumannBoundaryFV<TDomain>::update_subset_groups()
{
	for(size_t i = 0; i < m_vNumberData.size(); ++i)
		base_type::update_subset_groups(m_vNumberData[i]);
	for(size_t i = 0; i < m_vBNDNumberData.size(); ++i)
		base_type::update_subset_groups(m_vBNDNumberData[i]);
	for(size_t i = 0; i < m_vVectorData.size(); ++i)
		base_type::update_subset_groups(m_vVectorData[i]);
}


////////////////////////////////////////////////////////////////////////////////
//	assembling functions
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NeumannBoundaryFV<TDomain>::
prep_elem_loop(const ReferenceObjectID roid, const int si)
{
	update_subset_groups();
	m_si = si;

//	register subsetIndex at Geometry
	TFVGeom& geo = GeomProvider<TFVGeom>::get(m_lfeID,m_order);
	typedef typename reference_element_traits<TElem>::reference_element_type reference_element_type;

	try{
		geo.update_local(roid, m_lfeID);
	}
	UG_CATCH_THROW("NeumannBoundaryFV::prep_elem_loop:"
						" Cannot update Finite Volume Geometry.");

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

	typedef typename NeumannBoundaryFV<TDomain>::NumberData T;
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;

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
void NeumannBoundaryFV<TDomain>::
prep_elem(TElem* elem, const LocalVector& u)
{
//  update Geometry for this element
	TFVGeom& geo = GeomProvider<TFVGeom>::get(m_lfeID,m_order);
	try{
		geo.update(elem,
	               this->template element_corners<TElem>(elem),
	               &(this->subset_handler()));
	}
	UG_CATCH_THROW("NeumannBoundaryFV::prep_elem: "
						"Cannot update Finite Volume Geometry.");

	for(size_t i = 0; i < m_vNumberData.size(); ++i)
		if(m_vNumberData[i].InnerSSGrp.contains(m_si))
			m_vNumberData[i].template extract_bip<TElem, TFVGeom>(geo);
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NeumannBoundaryFV<TDomain>::
add_rhs_elem(LocalVector& d)
{
	TFVGeom& geo = GeomProvider<TFVGeom>::get(m_lfeID,m_order);
	typedef typename TFVGeom::BF BF;

//	Number Data
	for(size_t data = 0; data < m_vNumberData.size(); ++data){
		if(!m_vNumberData[data].InnerSSGrp.contains(m_si)) continue;
		size_t ip = 0;
		for(size_t s = 0; s < m_vNumberData[data].BndSSGrp.size(); ++s){
			const int si = m_vNumberData[data].BndSSGrp[s];
			const std::vector<BF>& vBF = geo.bf(si);

			for(size_t b = 0; b < vBF.size(); ++b){
				const int co = vBF[b].node_id();

				for(size_t i = 0; i < vBF[b].num_ip(); ++i, ++ip){
					d(_C_, co) -= m_vNumberData[data].import[ip]
					             * vBF[b].volume() * vBF[b].weight(i);
				}
			}
		}
	}

//	conditional Number Data
	for(size_t data = 0; data < m_vBNDNumberData.size(); ++data){
		if(!m_vBNDNumberData[data].InnerSSGrp.contains(m_si)) continue;
		for(size_t s = 0; s < m_vBNDNumberData[data].BndSSGrp.size(); ++s)	{
			const int si = m_vBNDNumberData[data].BndSSGrp[s];
			const std::vector<BF>& vBF = geo.bf(si);

			for(size_t b = 0; b < vBF.size(); ++b){
				number val = 0.0;
				const int co = vBF[b].node_id();

				for(size_t i = 0; i < vBF[b].num_ip(); ++i){
					if(!(*m_vBNDNumberData[data].functor)(val, vBF[b].global_ip(i), this->time(), si))
						continue;

					d(_C_, co) -= val * vBF[b].volume() * vBF[b].weight(i);
				}
			}
		}
	}

// 	vector data
	for(size_t data = 0; data < m_vVectorData.size(); ++data){
		if(!m_vVectorData[data].InnerSSGrp.contains(m_si)) continue;
		for(size_t s = 0; s < m_vVectorData[data].BndSSGrp.size(); ++s){
			const int si = m_vVectorData[data].BndSSGrp[s];
			const std::vector<BF>& vBF = geo.bf(si);

			for(size_t b = 0; b < vBF.size(); ++b){
				MathVector<dim> val;
				const int co = vBF[b].node_id();

				for(size_t i = 0; i < vBF[b].num_ip(); ++i){
					(*m_vVectorData[data].functor)(val, vBF[b].global_ip(i), this->time(), si);

					d(_C_, co) -= vBF[b].weight(i) * VecDot(val, vBF[b].normal());
				}
			}
		}
	}
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NeumannBoundaryFV<TDomain>::
finish_elem_loop()
{
//	remove subsetIndex from Geometry
	TFVGeom& geo = GeomProvider<TFVGeom>::get(m_lfeID,m_order);

//	unrequest subset indices as boundary subset. This will force the
//	creation of boundary subsets when calling geo.update

	for(size_t i = 0; i < m_vNumberData.size(); ++i){
		if(!m_vNumberData[i].InnerSSGrp.contains(m_si)) continue;
		for(size_t s = 0; s < m_vNumberData[i].BndSSGrp.size(); ++s){
			const int si = m_vNumberData[i].BndSSGrp[s];
			geo.remove_boundary_subset(si);
		}
	}

	for(size_t i = 0; i < m_vBNDNumberData.size(); ++i){
		if(!m_vBNDNumberData[i].InnerSSGrp.contains(m_si)) continue;
		for(size_t s = 0; s < m_vBNDNumberData[i].BndSSGrp.size(); ++s){
			const int si = m_vBNDNumberData[i].BndSSGrp[s];
			geo.remove_boundary_subset(si);
		}
	}

	for(size_t i = 0; i < m_vVectorData.size(); ++i){
		if(!m_vVectorData[i].InnerSSGrp.contains(m_si)) continue;
		for(size_t s = 0; s < m_vVectorData[i].BndSSGrp.size(); ++s){
			const int si = m_vVectorData[i].BndSSGrp[s];
			geo.remove_boundary_subset(si);
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
// Number Data
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NeumannBoundaryFV<TDomain>::NumberData::
lin_def(const LocalVector& u,
            std::vector<std::vector<number> > vvvLinDef[],
            const size_t nip)
{
//  get finite volume geometry
	const TFVGeom& geo = GeomProvider<TFVGeom>::get(This->m_lfeID,This->m_order);
	typedef typename TFVGeom::BF BF;

	size_t ip = 0;
	for(size_t s = 0; s < this->BndSSGrp.size(); ++s)
	{
		const int si = this->BndSSGrp[s];
		const std::vector<BF>& vBF = geo.bf(si);
		for(size_t i = 0; i < vBF.size(); ++i){
			const int co = vBF[i].node_id();
			vvvLinDef[ip][_C_][co] -= vBF[i].volume();
		}
	}
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NeumannBoundaryFV<TDomain>::NumberData::
extract_bip(const TFVGeom& geo)
{
	typedef typename TFVGeom::BF BF;
	vLocIP.clear();
	vGloIP.clear();
	for(size_t s = 0; s < this->BndSSGrp.size(); s++)
	{
		const int si = this->BndSSGrp[s];
		const std::vector<BF>& vBF = geo.bf(si);
		for(size_t i = 0; i < vBF.size(); ++i)
		{
			const BF& bf = vBF[i];
			for(size_t ip = 0; ip < bf.num_ip(); ++ip){
				vLocIP.push_back(bf.local_ip(ip));
				vGloIP.push_back(bf.global_ip(ip));
			}
		}
	}

	import.set_local_ips(&vLocIP[0], vLocIP.size());
	import.set_global_ips(&vGloIP[0], vGloIP.size());
}
////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template<>
void NeumannBoundaryFV<Domain1d>::register_all_funcs(int order)
{
}
#endif

#ifdef UG_DIM_2
template<>
void NeumannBoundaryFV<Domain2d>::register_all_funcs(int order)
{
//	Triangle
	switch(order)
	{
		case 1:	{typedef FVGeometry<1, Triangle, dim> FVGeom;
				 register_func<Triangle, FVGeom >(); break;}
		case 2:	{typedef FVGeometry<2, Triangle, dim> FVGeom;
				 register_func<Triangle, FVGeom >(); break;}
		case 3:	{typedef FVGeometry<3, Triangle, dim> FVGeom;
				 register_func<Triangle, FVGeom >(); break;}
		default: {typedef DimFVGeometry<dim, 2> FVGeom;
				 register_func<Triangle, FVGeom >(); break;}
	}

//	Quadrilateral
	switch(order) {
		case 1:	{typedef FVGeometry<1, Quadrilateral, dim> FVGeom;
				 register_func<Quadrilateral, FVGeom >(); break;}
		case 2:	{typedef FVGeometry<2, Quadrilateral, dim> FVGeom;
				 register_func<Quadrilateral, FVGeom >(); break;}
		case 3:	{typedef FVGeometry<3, Quadrilateral, dim> FVGeom;
				 register_func<Quadrilateral, FVGeom >(); break;}
		default: {typedef DimFVGeometry<dim, 2> FVGeom;
				  register_func<Quadrilateral, FVGeom >(); break;}
	}
}
#endif

#ifdef UG_DIM_3
template<>
void NeumannBoundaryFV<Domain3d>::register_all_funcs(int order)
{
//	Tetrahedron
	switch(order)
	{
		case 1:	{typedef FVGeometry<1, Tetrahedron, dim> FVGeom;
				 register_func<Tetrahedron, FVGeom >(); break;}
		case 2:	{typedef FVGeometry<2, Tetrahedron, dim> FVGeom;
				 register_func<Tetrahedron, FVGeom >(); break;}
		case 3:	{typedef FVGeometry<3, Tetrahedron, dim> FVGeom;
				 register_func<Tetrahedron, FVGeom >(); break;}
		default: {typedef DimFVGeometry<dim, 3> FVGeom;
				  register_func<Tetrahedron, FVGeom >(); break;}
	}

//	Prism
	switch(order) {
		case 1:	{typedef FVGeometry<1, Prism, dim> FVGeom;
				 register_func<Prism, FVGeom >(); break;}
		default: {typedef DimFVGeometry<dim, 3> FVGeom;
				  register_func<Prism, FVGeom >(); break;}
	}

//	Hexahedron
	switch(order)
	{
		case 1:	{typedef FVGeometry<1, Hexahedron, dim> FVGeom;
				 register_func<Hexahedron, FVGeom >(); break;}
		case 2:	{typedef FVGeometry<2, Hexahedron, dim> FVGeom;
				 register_func<Hexahedron, FVGeom >(); break;}
		case 3:	{typedef FVGeometry<3, Hexahedron, dim> FVGeom;
				 register_func<Hexahedron, FVGeom >(); break;}
		default: {typedef DimFVGeometry<dim, 3> FVGeom;
				  register_func<Hexahedron, FVGeom >(); break;}
	}
}
#endif

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NeumannBoundaryFV<TDomain>::register_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef this_type T;

	this->enable_fast_add_elem(true);

	this->set_prep_elem_loop_fct(id, &T::template prep_elem_loop<TElem, TFVGeom>);
	this->set_prep_elem_fct(	 id, &T::template prep_elem<TElem, TFVGeom>);
	this->set_add_rhs_elem_fct(	 id, &T::template add_rhs_elem<TElem, TFVGeom>);
	this->set_fsh_elem_loop_fct( id, &T::template finish_elem_loop<TElem, TFVGeom>);

	this->set_add_jac_A_elem_fct(	 id, &T::template add_JA_elem<TElem, TFVGeom>);
	this->set_add_jac_M_elem_fct(	 id, &T::template add_JM_elem<TElem, TFVGeom>);
	this->set_add_def_A_elem_fct(	 id, &T::template add_dA_elem<TElem, TFVGeom>);
	this->set_add_def_M_elem_fct(	 id, &T::template add_dM_elem<TElem, TFVGeom>);
}

////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template class NeumannBoundaryFV<Domain1d>;
#endif
#ifdef UG_DIM_2
template class NeumannBoundaryFV<Domain2d>;
#endif
#ifdef UG_DIM_3
template class NeumannBoundaryFV<Domain3d>;
#endif

} // namespace ug

