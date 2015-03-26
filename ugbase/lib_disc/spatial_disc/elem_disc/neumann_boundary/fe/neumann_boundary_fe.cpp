/*
 * neumann_boundary_fe.cpp
 *
 *  Created on: 14.10.2010
 *      Author: andreasvogel
 */

#include "neumann_boundary_fe.h"
#include "lib_disc/spatial_disc/disc_util/fe_geom.h"
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
//	Constructor
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
NeumannBoundaryFE<TDomain>::NeumannBoundaryFE(const char* function)
 :NeumannBoundaryBase<TDomain>(function),
  m_order(1), m_lfeID(LFEID::LAGRANGE, TDomain::dim, m_order)
{
	this->clear_add_fct();
}

template<typename TDomain>
void NeumannBoundaryFE<TDomain>::
prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid)
{
//	check number
	if(vLfeID.size() != 1)
		UG_THROW("NeumannBoundaryFE: needs exactly 1 function.");

//	check that not ADAPTIVE
	if(vLfeID[0].order() < 1)
		UG_THROW("NeumannBoundaryFE: Adaptive order not implemented.");

//	set order
	m_lfeID = vLfeID[0];
	m_order = vLfeID[0].order();
	m_quadOrder = 2*m_order+1;

	register_all_funcs(m_order);
}

template<typename TDomain>
void NeumannBoundaryFE<TDomain>::
add(SmartPtr<CplUserData<number, dim> > data, const char* BndSubsets, const char* InnerSubsets)
{
	m_vNumberData.push_back(NumberData(data, BndSubsets, InnerSubsets, this));
	this->add_inner_subsets(InnerSubsets);
}

template<typename TDomain>
void NeumannBoundaryFE<TDomain>::
add(SmartPtr<CplUserData<number, dim, bool> > user, const char* BndSubsets, const char* InnerSubsets)
{
	m_vBNDNumberData.push_back(BNDNumberData(user, BndSubsets, InnerSubsets));
	this->add_inner_subsets(InnerSubsets);
}

template<typename TDomain>
void NeumannBoundaryFE<TDomain>::
add(SmartPtr<CplUserData<MathVector<dim>, dim> > user, const char* BndSubsets, const char* InnerSubsets)
{
	m_vVectorData.push_back(VectorData(user, BndSubsets, InnerSubsets));
	this->add_inner_subsets(InnerSubsets);
}

template<typename TDomain>
void NeumannBoundaryFE<TDomain>::update_subset_groups()
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
template<typename TElem, typename TFEGeom>
void NeumannBoundaryFE<TDomain>::
prep_elem_loop(const ReferenceObjectID roid, const int si)
{
	update_subset_groups();
	m_si = si;

//	register subsetIndex at Geometry
	TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

	try{
		geo.update_local(roid, m_lfeID, m_quadOrder);
	}
	UG_CATCH_THROW("NeumannBoundaryFE::prep_elem_loop:"
						" Cannot update Finite Element Geometry.");

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

	typedef typename NeumannBoundaryFE<TDomain>::NumberData T;
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;

//	set lin defect fct for imports
	for(size_t data = 0; data < m_vNumberData.size(); ++data)
	{
		if(!m_vNumberData[data].InnerSSGrp.contains(m_si)) continue;
		m_vNumberData[data].import.set_fct(id,
		                                   &m_vNumberData[data],
		                                   &NumberData::template lin_def<TElem, TFEGeom>);

		this->register_import(m_vNumberData[data].import);
		m_vNumberData[data].import.set_rhs_part();
	}
}

template<typename TDomain>
template<typename TElem, typename TFEGeom>
void NeumannBoundaryFE<TDomain>::
prep_elem(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[])
{
//  update Geometry for this element
	TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);
	try{
		geo.update_boundary_faces(elem, vCornerCoords,
	               m_quadOrder,
	               &(this->subset_handler()));
	}
	UG_CATCH_THROW("NeumannBoundaryFE::prep_elem: "
						"Cannot update Finite Element Geometry.");

	for(size_t i = 0; i < m_vNumberData.size(); ++i)
		if(m_vNumberData[i].InnerSSGrp.contains(m_si))
			m_vNumberData[i].template extract_bip<TElem, TFEGeom>(geo);
}

template<typename TDomain>
template<typename TElem, typename TFEGeom>
void NeumannBoundaryFE<TDomain>::
add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);
	typedef typename TFEGeom::BF BF;

//	Number Data
	for(size_t data = 0; data < m_vNumberData.size(); ++data){
		if(!m_vNumberData[data].InnerSSGrp.contains(m_si)) continue;
		for(size_t s = 0; s < m_vNumberData[data].BndSSGrp.size(); ++s){
			const int si = m_vNumberData[data].BndSSGrp[s];
			const std::vector<BF>& vBF = geo.bf(si);

			for(size_t b = 0; b < vBF.size(); ++b){
				for(size_t ip = 0; ip < vBF[b].num_ip(); ++ip){
					for(size_t sh = 0; sh < vBF[b].num_sh(); ++sh){
						d(_C_, sh) -= m_vNumberData[data].import[ip]
									 * vBF[b].shape(ip, sh) * vBF[b].weight(ip);
					}
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
				for(size_t ip = 0; ip < vBF[b].num_ip(); ++ip){
					if(!(*m_vBNDNumberData[data].functor)(val, vBF[b].global_ip(ip), this->time(), si))
						continue;

					for(size_t sh = 0; sh < vBF[b].num_sh(); ++sh)
						d(_C_, sh) -= val * vBF[b].shape(ip, sh) * vBF[b].weight(ip);
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
				for(size_t ip = 0; ip < vBF[b].num_ip(); ++ip){
					(*m_vVectorData[data].functor)(val, vBF[b].global_ip(ip), this->time(), si);

					for(size_t sh = 0; sh < vBF[b].num_sh(); ++sh)
						d(_C_, sh) -= vBF[b].shape(ip, sh) * vBF[b].weight(ip) * VecDot(val, vBF[b].normal());
				}
			}
		}
	}
}

template<typename TDomain>
template<typename TElem, typename TFEGeom>
void NeumannBoundaryFE<TDomain>::
finish_elem_loop()
{
//	remove subsetIndex from Geometry
	TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID,m_order);

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
template<typename TElem, typename TFEGeom>
void NeumannBoundaryFE<TDomain>::NumberData::
lin_def(const LocalVector& u,
            std::vector<std::vector<number> > vvvLinDef[],
            const size_t nip)
{
//  get finite volume geometry
	const TFEGeom& geo = GeomProvider<TFEGeom>::get(This->m_lfeID,This->m_order);
	typedef typename TFEGeom::BF BF;

	for(size_t s = 0; s < this->BndSSGrp.size(); ++s)
	{
		const int si = this->BndSSGrp[s];
		const std::vector<BF>& vBF = geo.bf(si);
		for(size_t b = 0; b < vBF.size(); ++b){
			for(size_t ip = 0; ip < vBF[b].num_ip(); ++ip){
				for(size_t sh = 0; sh < vBF[b].num_sh(); ++sh)
					vvvLinDef[ip][_C_][sh] -= vBF[b].weight(ip) * vBF[b].shape(ip, sh);
			}
		}
	}
}

template<typename TDomain>
template<typename TElem, typename TFEGeom>
void NeumannBoundaryFE<TDomain>::NumberData::
extract_bip(const TFEGeom& geo)
{
	typedef typename TFEGeom::BF BF;
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
void NeumannBoundaryFE<Domain1d>::register_all_funcs(int order)
{
	register_func<RegularEdge, DimFEGeometry<dim, 1> >();
}
#endif

#ifdef UG_DIM_2
template<>
void NeumannBoundaryFE<Domain2d>::register_all_funcs(int order)
{
	register_func<Triangle, DimFEGeometry<dim, 2> >();
	register_func<Quadrilateral, DimFEGeometry<dim, 2> >();
}
#endif

#ifdef UG_DIM_3
template<>
void NeumannBoundaryFE<Domain3d>::register_all_funcs(int order)
{
	register_func<Tetrahedron, DimFEGeometry<dim, 3> >();
	register_func<Prism, DimFEGeometry<dim, 3> >();
	register_func<Pyramid, DimFEGeometry<dim, 3> >();
	register_func<Hexahedron, DimFEGeometry<dim, 3> >();
	register_func<Octahedron, DimFEGeometry<dim, 3> >();
}
#endif

template<typename TDomain>
template<typename TElem, typename TFEGeom>
void NeumannBoundaryFE<TDomain>::register_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef this_type T;

	this->clear_add_fct(id);

	this->set_prep_elem_loop_fct(id, &T::template prep_elem_loop<TElem, TFEGeom>);
	this->set_prep_elem_fct(	 id, &T::template prep_elem<TElem, TFEGeom>);
	this->set_add_rhs_elem_fct(	 id, &T::template add_rhs_elem<TElem, TFEGeom>);
	this->set_fsh_elem_loop_fct( id, &T::template finish_elem_loop<TElem, TFEGeom>);

	this->set_add_jac_A_elem_fct(	 id, &T::template add_jac_A_elem<TElem, TFEGeom>);
	this->set_add_jac_M_elem_fct(	 id, &T::template add_jac_M_elem<TElem, TFEGeom>);
	this->set_add_def_A_elem_fct(	 id, &T::template add_def_A_elem<TElem, TFEGeom>);
	this->set_add_def_M_elem_fct(	 id, &T::template add_def_M_elem<TElem, TFEGeom>);
}

////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template class NeumannBoundaryFE<Domain1d>;
#endif
#ifdef UG_DIM_2
template class NeumannBoundaryFE<Domain2d>;
#endif
#ifdef UG_DIM_3
template class NeumannBoundaryFE<Domain3d>;
#endif

} // namespace ug

