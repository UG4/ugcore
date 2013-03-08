/*
 * neumann_boundary_fv1.cpp
 *
 *  Created on: 14.10.2010
 *      Author: andreasvogel
 */

#include "neumann_boundary_fv1.h"
#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"

namespace ug{


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
bool NeumannBoundaryFV1<TDomain>::
request_finite_element_id(const std::vector<LFEID>& vLfeID)
{
//	check number
	if(vLfeID.size() != 1)
	{
		UG_LOG("NeumannBoundaryFV1:"
				" Wrong number of functions given. Need exactly "<<1<<"\n");
		return false;
	}

	if(vLfeID[0].order() != 1 || vLfeID[0].type() != LFEID::LAGRANGE)
	{
		UG_LOG("NeumannBoundaryFV1:"
				" FV Scheme only implemented for 1st order.\n");
		return false;
	}

//	is supported
	return true;
}

template<typename TDomain>
bool NeumannBoundaryFV1<TDomain>::
request_non_regular_grid(bool bNonRegular)
{
	if(bNonRegular)
		UG_THROW("NeumannBoundary: Hanging Nodes not implemented.");

	return true;
}


template<typename TDomain>
void NeumannBoundaryFV1<TDomain>::
add(SmartPtr<UserData<number, dim> > data, const char* BndSubsets, const char* InnerSubsets)
{
	m_vNumberData.push_back(NumberData(data, BndSubsets, InnerSubsets));
	this->add_inner_subsets(InnerSubsets);
}

template<typename TDomain>
void NeumannBoundaryFV1<TDomain>::
add(SmartPtr<UserData<number, dim, bool> > user, const char* BndSubsets, const char* InnerSubsets)
{
	m_vBNDNumberData.push_back(BNDNumberData(user, BndSubsets, InnerSubsets));
	this->add_inner_subsets(InnerSubsets);
}

template<typename TDomain>
void NeumannBoundaryFV1<TDomain>::
add(SmartPtr<UserData<MathVector<dim>, dim> > user, const char* BndSubsets, const char* InnerSubsets)
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
prep_elem_loop(const ReferenceObjectID roid, const int si)
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

	typedef typename NeumannBoundaryFV1<TDomain>::NumberData T;
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
void NeumannBoundaryFV1<TDomain>::
prep_elem(TElem* elem, const LocalVector& u)
{
//  update Geometry for this element
	static TFVGeom& geo = GeomProvider<TFVGeom >::get();
	try{
		geo.update(elem,
	               this->template element_corners<TElem>(elem),
	               &(this->subset_handler()));
	}
	UG_CATCH_THROW("NeumannBoundary::prep_elem: "
						"Cannot update Finite Volume Geometry.");

	for(size_t i = 0; i < m_vNumberData.size(); ++i)
		if(m_vNumberData[i].InnerSSGrp.contains(m_si))
			m_vNumberData[i].template extract_bip<TElem, TFVGeom>(geo);
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NeumannBoundaryFV1<TDomain>::
add_rhs_elem(LocalVector& d)
{
	const static TFVGeom& geo = GeomProvider<TFVGeom >::get();
	typedef typename TFVGeom::BF BF;

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
finish_elem_loop()
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
void NeumannBoundaryFV1<TDomain>::NumberData::
lin_def(const LocalVector& u,
            std::vector<std::vector<number> > vvvLinDef[],
            const size_t nip)
{
//  get finite volume geometry
	const static TFVGeom& geo = GeomProvider<TFVGeom>::get();
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
void NeumannBoundaryFV1<TDomain>::NumberData::
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
			vLocIP.push_back(bf.local_ip());
			vGloIP.push_back(bf.global_ip());
		}
	}

	import.set_local_ips(&vLocIP[0], vLocIP.size());
	import.set_global_ips(&vGloIP[0], vGloIP.size());
}

////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////

template<>
void NeumannBoundaryFV1<Domain1d>::register_all_funcs(bool bHang)
{
	if(!bHang){
		register_func<Edge, FV1Geometry<Edge, dim> >();
	}
	else UG_THROW("NeumannBoundary: Hanging Nodes not implemented.");
}

template<>
void NeumannBoundaryFV1<Domain2d>::register_all_funcs(bool bHang)
{
	if(!bHang){
		register_func<Triangle, FV1Geometry<Triangle, dim> >();
		register_func<Quadrilateral, FV1Geometry<Quadrilateral, dim> >();
	}
	else UG_THROW("NeumannBoundary: Hanging Nodes not implemented.");
}

template<>
void NeumannBoundaryFV1<Domain3d>::register_all_funcs(bool bHang)
{
	if(!bHang){
		register_func<Tetrahedron, FV1Geometry<Tetrahedron, dim> >();
		register_func<Prism, FV1Geometry<Prism, dim> >();
		register_func<Pyramid, FV1Geometry<Pyramid, dim> >();
		register_func<Hexahedron, FV1Geometry<Hexahedron, dim> >();
	}
	else UG_THROW("NeumannBoundary: Hanging Nodes not implemented.");
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NeumannBoundaryFV1<TDomain>::register_func()
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
template class NeumannBoundaryFV1<Domain1d>;
#endif
#ifdef UG_DIM_2
template class NeumannBoundaryFV1<Domain2d>;
#endif
#ifdef UG_DIM_3
template class NeumannBoundaryFV1<Domain3d>;
#endif

} // namespace ug

