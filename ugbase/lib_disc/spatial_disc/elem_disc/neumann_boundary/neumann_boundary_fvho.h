/*
 * neumann_boundary_fvho.h
 *
 *  Created on: 14.10.2010
 *      Author: andreasvogel
 */

#include "neumann_boundary.h"
#include "lib_disc/spatial_disc/disc_util/fvho_geom.h"

namespace ug{

template<typename TDomain>
template<typename TElem, typename TGeomProvider>
void NeumannBoundary<TDomain>::
prep_elem_loop_fvho()
{
//	register subsetIndex at Geometry
	static typename TGeomProvider::Type& geo = TGeomProvider::get();

//	request subset indices as boundary subset. This will force the
//	creation of boundary subsets when calling geo.update

	for(size_t i = 0; i < m_vNumberData.size(); ++i){
		for(size_t s = 0; s < m_vNumberData[i].ssGrp.size(); ++s){
			const int si = m_vNumberData[i].ssGrp[s];
			geo.add_boundary_subset(si);
		}
	}
	for(size_t i = 0; i < m_vBNDNumberData.size(); ++i){
		for(size_t s = 0; s < m_vBNDNumberData[i].ssGrp.size(); ++s){
			const int si = m_vBNDNumberData[i].ssGrp[s];
			geo.add_boundary_subset(si);
		}
	}
	for(size_t i = 0; i < m_vVectorData.size(); ++i){
		for(size_t s = 0; s < m_vVectorData[i].ssGrp.size(); ++s){
			const int si = m_vVectorData[i].ssGrp[s];
			geo.add_boundary_subset(si);
		}
	}

//	clear imports, since we will set them afterwards
	this->clear_imports();

	typedef typename NeumannBoundary<TDomain>::NumberData T;
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;

//	set lin defect fct for imports
	for(size_t data = 0; data < m_vNumberData.size(); ++data)
	{
		m_vNumberData[data].import.set_fct(id,
		                                   &m_vNumberData[data],
		                                   &NumberData::template lin_def_fvho<TElem, TGeomProvider>);

		this->register_import(m_vNumberData[data].import);
		m_vNumberData[data].import.set_rhs_part();
	}
}

template<typename TDomain>
template<typename TElem, typename TGeomProvider>
void NeumannBoundary<TDomain>::
prep_elem_fvho(TElem* elem, const LocalVector& u)
{
//  update Geometry for this element
	static typename TGeomProvider::Type& geo = TGeomProvider::get();
	try{
		geo.update(elem,
	               this->template element_corners<TElem>(elem),
	               &(this->subset_handler()));
	}
	UG_CATCH_THROW("NeumannBoundary::prep_elem_fvho: "
						"Cannot update Finite Volume Geometry.");

	for(size_t i = 0; i < m_vNumberData.size(); ++i)
		m_vNumberData[i].template extract_bip_fvho<TElem, typename TGeomProvider::Type>(geo);
}

template<typename TDomain>
template<typename TElem, typename TGeomProvider>
void NeumannBoundary<TDomain>::
add_rhs_elem_fvho(LocalVector& d)
{
	static const typename TGeomProvider::Type& geo = TGeomProvider::get();
	typedef typename TGeomProvider::Type::BF BF;

//	Number Data
	for(size_t data = 0; data < m_vNumberData.size(); ++data){
		size_t ip = 0;
		for(size_t s = 0; s < m_vNumberData[data].ssGrp.size(); ++s){
			const int si = m_vNumberData[data].ssGrp[s];
			const std::vector<BF>& vBF = geo.bf(si);

			for(size_t b = 0; b < vBF.size(); ++b, ++ip){
				const int co = vBF[b].node_id();

				for(size_t i = 0; i < vBF[b].num_ip(); ++i){
					d(m_vNumberData[data].locFct, co) -= m_vNumberData[data].import[ip]
														* vBF[b].volume()
														* vBF[b].weight(ip);
				}
			}
		}
	}

//	conditional Number Data
	for(size_t data = 0; data < m_vBNDNumberData.size(); ++data){
		for(size_t s = 0; s < m_vBNDNumberData[data].ssGrp.size(); ++s)	{
			const int si = m_vBNDNumberData[data].ssGrp[s];
			const std::vector<BF>& vBF = geo.bf(si);

			for(size_t b = 0; b < vBF.size(); ++b){
				number val = 0.0;
				const int co = vBF[b].node_id();

				for(size_t i = 0; i < vBF[b].num_ip(); ++i){
					if(!(*m_vBNDNumberData[data].functor)(val, vBF[b].global_ip(i), this->time(), si))
						continue;

					d(m_vBNDNumberData[data].locFct, co) -= val
															* vBF[b].volume()
															* vBF[b].weight(i);
				}
			}
		}
	}

// 	vector data
	for(size_t data = 0; data < m_vVectorData.size(); ++data){
		for(size_t s = 0; s < m_vVectorData[data].ssGrp.size(); ++s){
			const int si = m_vVectorData[data].ssGrp[s];
			const std::vector<BF>& vBF = geo.bf(si);

			for(size_t b = 0; b < vBF.size(); ++b){
				MathVector<dim> val;
				const int co = vBF[b].node_id();

				for(size_t i = 0; i < vBF[b].num_ip(); ++i){
					(*m_vVectorData[data].functor)(val, vBF[b].global_ip(i), this->time(), si);

					d(m_vVectorData[data].locFct, co) -= vBF[b].weight(i) *
														 VecDot(val, vBF[b].normal());
				}
			}
		}
	}
}

template<typename TDomain>
template<typename TElem, typename TGeomProvider>
void NeumannBoundary<TDomain>::
finish_elem_loop_fvho()
{
//	remove subsetIndex from Geometry
	static typename TGeomProvider::Type& geo = TGeomProvider::get();

//	unrequest subset indices as boundary subset. This will force the
//	creation of boundary subsets when calling geo.update

	for(size_t i = 0; i < m_vNumberData.size(); ++i){
		for(size_t s = 0; s < m_vNumberData[i].ssGrp.size(); ++s){
			const int si = m_vNumberData[i].ssGrp[s];
			geo.remove_boundary_subset(si);
		}
	}

	for(size_t i = 0; i < m_vBNDNumberData.size(); ++i){
		for(size_t s = 0; s < m_vBNDNumberData[i].ssGrp.size(); ++s){
			const int si = m_vBNDNumberData[i].ssGrp[s];
			geo.remove_boundary_subset(si);
		}
	}

	for(size_t i = 0; i < m_vVectorData.size(); ++i){
		for(size_t s = 0; s < m_vVectorData[i].ssGrp.size(); ++s){
			const int si = m_vVectorData[i].ssGrp[s];
			geo.remove_boundary_subset(si);
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
// Number Data
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
template<typename TElem, typename TGeomProvider>
void NeumannBoundary<TDomain>::NumberData::
lin_def_fvho(const LocalVector& u,
            std::vector<std::vector<number> > vvvLinDef[],
            const size_t nip)
{
//  get finite volume geometry
	static const typename TGeomProvider::Type& geo = TGeomProvider::get();
	typedef typename TGeomProvider::Type::BF BF;

	size_t ip = 0;
	for(size_t s = 0; s < this->ssGrp.size(); ++s)
	{
		const int si = this->ssGrp[s];
		const std::vector<BF>& vBF = geo.bf(si);
		for(size_t i = 0; i < vBF.size(); ++i){
			const int co = vBF[i].node_id();
			vvvLinDef[ip][this->locFct][co] -= vBF[i].volume();
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////

template<>
void NeumannBoundary<Domain1d>::register_all_fvho_funcs(int order)
{
}

template<>
void NeumannBoundary<Domain2d>::register_all_fvho_funcs(int order)
{
//	Triangle
	switch(order)
	{
		case 1:	{typedef FVGeometry<1, Triangle, dim> FVGeom;
				 register_fvho_func<Triangle, Provider<FVGeom> >(); break;}
		case 2:	{typedef FVGeometry<2, Triangle, dim> FVGeom;
				 register_fvho_func<Triangle, Provider<FVGeom> >(); break;}
		case 3:	{typedef FVGeometry<3, Triangle, dim> FVGeom;
				 register_fvho_func<Triangle, Provider<FVGeom> >(); break;}
		default: {typedef DimFVGeometry<2, dim> FVGeom;
				 register_fvho_func<Triangle, FlexGeomProvider<FVGeom> >(); break;}
	}

//	Quadrilateral
	switch(order) {
		case 1:	{typedef FVGeometry<1, Quadrilateral, dim> FVGeom;
				 register_fvho_func<Quadrilateral, Provider<FVGeom> >(); break;}
		case 2:	{typedef FVGeometry<2, Quadrilateral, dim> FVGeom;
				 register_fvho_func<Quadrilateral, Provider<FVGeom> >(); break;}
		case 3:	{typedef FVGeometry<3, Quadrilateral, dim> FVGeom;
				 register_fvho_func<Quadrilateral, Provider<FVGeom> >(); break;}
		default: {typedef DimFVGeometry<2, dim> FVGeom;
				  register_fvho_func<Quadrilateral, FlexGeomProvider<FVGeom> >(); break;}
	}
}

template<>
void NeumannBoundary<Domain3d>::register_all_fvho_funcs(int order)
{
//	Tetrahedron
	switch(order)
	{
		case 1:	{typedef FVGeometry<1, Tetrahedron, dim> FVGeom;
				 register_fvho_func<Tetrahedron, Provider<FVGeom> >(); break;}
		case 2:	{typedef FVGeometry<2, Tetrahedron, dim> FVGeom;
				 register_fvho_func<Tetrahedron, Provider<FVGeom> >(); break;}
		case 3:	{typedef FVGeometry<3, Tetrahedron, dim> FVGeom;
				 register_fvho_func<Tetrahedron, Provider<FVGeom> >(); break;}
		default: {typedef DimFVGeometry<3, dim> FVGeom;
				  register_fvho_func<Tetrahedron, FlexGeomProvider<FVGeom> >(); break;}
	}

//	Prism
	switch(order) {
		case 1:	{typedef FVGeometry<1, Prism, dim> FVGeom;
				 register_fvho_func<Prism, Provider<FVGeom> >(); break;}
		default: {typedef DimFVGeometry<3, dim> FVGeom;
				  register_fvho_func<Prism, FlexGeomProvider<FVGeom> >(); break;}
	}

//	Hexahedron
	switch(order)
	{
		case 1:	{typedef FVGeometry<1, Hexahedron, dim> FVGeom;
				 register_fvho_func<Hexahedron, Provider<FVGeom> >(); break;}
		case 2:	{typedef FVGeometry<2, Hexahedron, dim> FVGeom;
				 register_fvho_func<Hexahedron, Provider<FVGeom> >(); break;}
		case 3:	{typedef FVGeometry<3, Hexahedron, dim> FVGeom;
				 register_fvho_func<Hexahedron, Provider<FVGeom> >(); break;}
		default: {typedef DimFVGeometry<3, dim> FVGeom;
				  register_fvho_func<Hexahedron, FlexGeomProvider<FVGeom> >(); break;}
	}
}


template<typename TDomain>
template<typename TElem, typename TGeomProvider>
void NeumannBoundary<TDomain>::register_fvho_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef this_type T;

	this->enable_fast_add_elem(true);

	this->set_prep_elem_loop_fct(id, &T::template prep_elem_loop_fvho<TElem, TGeomProvider>);
	this->set_prep_elem_fct(	 id, &T::template prep_elem_fvho<TElem, TGeomProvider>);
	this->set_add_rhs_elem_fct(	 id, &T::template add_rhs_elem_fvho<TElem, TGeomProvider>);
	this->set_fsh_elem_loop_fct( id, &T::template finish_elem_loop_fvho<TElem, TGeomProvider>);

	this->set_add_jac_A_elem_fct(	 id, &T::template add_JA_elem<TElem, TGeomProvider>);
	this->set_add_jac_M_elem_fct(	 id, &T::template add_JM_elem<TElem, TGeomProvider>);
	this->set_add_def_A_elem_fct(	 id, &T::template add_dA_elem<TElem, TGeomProvider>);
	this->set_add_def_M_elem_fct(	 id, &T::template add_dM_elem<TElem, TGeomProvider>);
}
} // namespace ug

