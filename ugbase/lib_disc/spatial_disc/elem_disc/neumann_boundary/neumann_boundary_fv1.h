/*
 * neumann_boundary_fv1.h
 *
 *  Created on: 14.10.2010
 *      Author: andreasvogel
 */

#include "neumann_boundary.h"
#include "lib_disc/spatial_disc/disc_util/finite_volume_geometry.h"

namespace ug{

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NeumannBoundary<TDomain>::
prep_elem_loop_fv1()
{
//	register subsetIndex at Geometry
	static TFVGeom& geo = Provider<TFVGeom >::get();

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
		                                   &NumberData::template lin_def_fv1<TElem, TFVGeom>);

		this->register_import(m_vNumberData[data].import);
		m_vNumberData[data].import.set_rhs_part();
	}
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NeumannBoundary<TDomain>::
prep_elem_fv1(TElem* elem, const LocalVector& u)
{
//  update Geometry for this element
	static TFVGeom& geo = Provider<TFVGeom >::get();
	if(!geo.update(elem,
	               this->template element_corners<TElem>(elem),
	               &(this->subset_handler())))
		UG_THROW("NeumannBoundary::prep_elem_fv1: "
						"Cannot update Finite Volume Geometry.");

	for(size_t i = 0; i < m_vNumberData.size(); ++i)
		m_vNumberData[i].template extract_bip<TElem, TFVGeom>(geo);
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NeumannBoundary<TDomain>::
add_rhs_elem_fv1(LocalVector& d)
{
	const static TFVGeom& geo = Provider<TFVGeom >::get();
	typedef typename TFVGeom::BF BF;

//	Number Data
	for(size_t data = 0; data < m_vNumberData.size(); ++data){
		size_t ip = 0;
		for(size_t s = 0; s < m_vNumberData[data].ssGrp.size(); ++s){
			const int si = m_vNumberData[data].ssGrp[s];
			const std::vector<BF>& vBF = geo.bf(si);

			for(size_t i = 0; i < vBF.size(); ++i, ++ip){
				const int co = vBF[i].node_id();
				d(m_vNumberData[data].locFct, co) -= m_vNumberData[data].import[ip]
				                                    * vBF[i].volume();
			}
		}
	}

//	conditional Number Data
	for(size_t data = 0; data < m_vBNDNumberData.size(); ++data){
		for(size_t s = 0; s < m_vBNDNumberData[data].ssGrp.size(); ++s)	{
			const int si = m_vBNDNumberData[data].ssGrp[s];
			const std::vector<BF>& vBF = geo.bf(si);

			for(size_t i = 0; i < vBF.size(); ++i){
				number val = 0.0;
				if(!(*m_vBNDNumberData[data].functor)(val, vBF[i].global_ip(), this->time(), si))
					continue;

				const int co = vBF[i].node_id();
				d(m_vBNDNumberData[data].locFct, co) -= val * vBF[i].volume();
			}
		}
	}

// 	vector data
	for(size_t data = 0; data < m_vVectorData.size(); ++data){
		for(size_t s = 0; s < m_vVectorData[data].ssGrp.size(); ++s){
			const int si = m_vVectorData[data].ssGrp[s];
			const std::vector<BF>& vBF = geo.bf(si);

			for(size_t i = 0; i < vBF.size(); ++i){
				MathVector<dim> val;
				(*m_vVectorData[data].functor)(val, vBF[i].global_ip(), this->time(), si);

				const int co = vBF[i].node_id();
				d(m_vVectorData[data].locFct, co) -= VecDot(val, vBF[i].normal());
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
// Number Data
////////////////////////////////////////////////////////////////////////////////


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NeumannBoundary<TDomain>::NumberData::
lin_def_fv1(const LocalVector& u,
            std::vector<std::vector<number> > vvvLinDef[],
            const size_t nip)
{
//  get finite volume geometry
	const static TFVGeom& geo = Provider<TFVGeom>::get();
	typedef typename TFVGeom::BF BF;

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
void NeumannBoundary<Domain1d>::register_all_fv1_funcs(bool bHang)
{
	if(!bHang){
		register_fv1_func<Edge, FV1Geometry<Edge, dim> >();
	}
	else UG_THROW("NeumannBoundary: Hanging Nodes not implemented.");
}

template<>
void NeumannBoundary<Domain2d>::register_all_fv1_funcs(bool bHang)
{
	if(!bHang){
		register_fv1_func<Triangle, FV1Geometry<Triangle, dim> >();
		register_fv1_func<Quadrilateral, FV1Geometry<Quadrilateral, dim> >();
	}
	else UG_THROW("NeumannBoundary: Hanging Nodes not implemented.");
}

template<>
void NeumannBoundary<Domain3d>::register_all_fv1_funcs(bool bHang)
{
	if(!bHang){
		register_fv1_func<Tetrahedron, FV1Geometry<Tetrahedron, dim> >();
		register_fv1_func<Prism, FV1Geometry<Prism, dim> >();
		register_fv1_func<Pyramid, FV1Geometry<Pyramid, dim> >();
		register_fv1_func<Hexahedron, FV1Geometry<Hexahedron, dim> >();
	}
	else UG_THROW("NeumannBoundary: Hanging Nodes not implemented.");
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void NeumannBoundary<TDomain>::register_fv1_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef this_type T;

	this->enable_fast_ass_elem(true);
	this->set_prep_elem_loop_fct(id, &T::template prep_elem_loop_fv1<TElem, TFVGeom>);
	this->set_prep_elem_fct(	 id, &T::template prep_elem_fv1<TElem, TFVGeom>);
	this->set_fsh_elem_loop_fct( id, &T::template finish_elem_loop<TElem, TFVGeom>);
	this->set_ass_JA_elem_fct(	 id, &T::template add_JA_elem_fv1<TElem, TFVGeom>);
	this->set_ass_JM_elem_fct(	 id, &T::template add_JM_elem_fv1<TElem, TFVGeom>);
	this->set_ass_dA_elem_fct(	 id, &T::template add_dA_elem_fv1<TElem, TFVGeom>);
	this->set_ass_dM_elem_fct(	 id, &T::template add_dM_elem_fv1<TElem, TFVGeom>);
	this->set_ass_rhs_elem_fct(	 id, &T::template add_rhs_elem_fv1<TElem, TFVGeom>);
}
} // namespace ug

