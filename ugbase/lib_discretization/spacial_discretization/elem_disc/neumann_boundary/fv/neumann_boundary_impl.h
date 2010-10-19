/*
 * neumann_boundary_impl.h
 *
 *  Created on: 14.10.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__NEUMANN_BOUNDARY__FV1__NEUMANN_BOUNDARY_IMPL__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__NEUMANN_BOUNDARY__FV1__NEUMANN_BOUNDARY_IMPL__

#include "neumann_boundary.h"

namespace ug{


template<template <class TElem, int TWorldDim> class TFVGeom, typename TDomain, typename TAlgebra>
template<typename TElem >
inline
bool
FVNeumannBoundaryElemDisc<TFVGeom, TDomain, TAlgebra>::
prepare_element_loop()
{
// 	resize corner coordinates
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
	m_vCornerCoords.resize(ref_elem_type::num_corners);

//  check domain
	if(m_pDomain == NULL)
	{
		UG_LOG("ERROR in 'FVNeumannBoundaryElemDisc::prepare_element_loop':"
				" Domain not set.");
		return false;
	}

// 	remember position attachement
	m_aaPos = m_pDomain->get_position_accessor();

//	register subsetIndex at Geometry
	TFVGeom<TElem, dim>& geo = FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();

	typename std::map<int, std::vector<UserDataFunction> >::const_iterator subsetIter;
	for(subsetIter = m_mBoundarySegment.begin(); subsetIter != m_mBoundarySegment.end(); ++subsetIter)
	{
		const int bndSubset = (*subsetIter).first;

		geo.add_boundary_subset(bndSubset);
	}

//	we're done
	return true;
}

template<template <class TElem, int TWorldDim> class TFVGeom, typename TDomain, typename TAlgebra>
template<typename TElem >
inline
bool
FVNeumannBoundaryElemDisc<TFVGeom, TDomain, TAlgebra>::
finish_element_loop()
{
//	remove subsetIndex from Geometry
	TFVGeom<TElem, dim>& geo = FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();

	typename std::map<int, std::vector<UserDataFunction> >::const_iterator subsetIter;
	for(subsetIter = m_mBoundarySegment.begin(); subsetIter != m_mBoundarySegment.end(); ++subsetIter)
	{
		const int bndSubset = (*subsetIter).first;

		geo.remove_boundary_subset(bndSubset);
	}

//	we're done
	return true;
}

template<template <class TElem, int TWorldDim> class TFVGeom, typename TDomain, typename TAlgebra>

template<typename TElem >
inline
bool
FVNeumannBoundaryElemDisc<TFVGeom, TDomain, TAlgebra>::
prepare_element(TElem* elem, const local_vector_type& u, const local_index_type& glob_ind)
{
	// load corners of this element
	for(size_t i = 0; i < m_vCornerCoords.size(); ++i)
	{
		VertexBase* vert = elem->vertex(i);
		m_vCornerCoords[i] = m_aaPos[vert];
	}

	// update Geometry for this element
	TFVGeom<TElem, dim>& geo = FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();
	if(!geo.update(elem, m_pDomain->get_subset_handler(), &m_vCornerCoords[0]))
		{UG_LOG("FVNeumannBoundaryElemDisc::prepare_element: Cannot update Finite Volume Geometry.\n"); return false;}

	return true;
}

template<template <class TElem, int TWorldDim> class TFVGeom, typename TDomain, typename TAlgebra>
template<typename TElem >
inline
bool
FVNeumannBoundaryElemDisc<TFVGeom, TDomain, TAlgebra>::
assemble_JA(local_matrix_type& J, const local_vector_type& u, number time)
{
	// we're done
	return true;
}


template<template <class TElem, int TWorldDim> class TFVGeom, typename TDomain, typename TAlgebra>
template<typename TElem >
inline
bool
FVNeumannBoundaryElemDisc<TFVGeom, TDomain, TAlgebra>::
assemble_JM(local_matrix_type& J, const local_vector_type& u, number time)
{
	// we're done
	return true;
}


template<template <class TElem, int TWorldDim> class TFVGeom, typename TDomain, typename TAlgebra>
template<typename TElem >
inline
bool
FVNeumannBoundaryElemDisc<TFVGeom, TDomain, TAlgebra>::
assemble_A(local_vector_type& d, const local_vector_type& u, number time)
{
	// we're done
	return true;
}


template<template <class TElem, int TWorldDim> class TFVGeom, typename TDomain, typename TAlgebra>
template<typename TElem >
inline
bool
FVNeumannBoundaryElemDisc<TFVGeom, TDomain, TAlgebra>::
assemble_M(local_vector_type& d, const local_vector_type& u, number time)
{
	// we're done
	return true;
}


template<template <class TElem, int TWorldDim> class TFVGeom, typename TDomain, typename TAlgebra>
template<typename TElem >
inline
bool
FVNeumannBoundaryElemDisc<TFVGeom, TDomain, TAlgebra>::
assemble_f(local_vector_type& d, number time)
{
	// get finite volume geometry
	TFVGeom<TElem, dim>& geo = FVGeometryProvider::get_geom<TFVGeom, TElem, dim>();

	// loop registered boundary segments
	typename std::map<int, std::vector<UserDataFunction> >::const_iterator subsetIter;
	for(subsetIter = m_mBoundarySegment.begin(); subsetIter != m_mBoundarySegment.end(); ++subsetIter)
	{
		const int bndSubset = (*subsetIter).first;
		const std::vector<UserDataFunction>& vSegmentFunction = (*subsetIter).second;

		// loop Boundary Faces
		for(size_t i = 0; i < geo.num_bf(bndSubset); ++i)
		{
		// get current BF
			const typename TFVGeom<TElem, dim>::BF& bf = geo.bf(bndSubset, i);

		//	loop functions, where neumann bnd is set for this bndSubset
			for(size_t fct = 0; fct < vSegmentFunction.size(); ++fct)
			{
				// first value
				number val = 0.0;
				vSegmentFunction[fct].functor(val, bf.global_ip(0), time);

				// other values
				for(size_t ip = 1; ip < bf.num_ip(); ++ip)
				{
					number ip_val;
					vSegmentFunction[fct].functor(ip_val, bf.global_ip(ip), time);

					// TODO: add weights for integration
					val += ip_val;
				}

				// scale with volume of BF
				val *= bf.volume();

				// get associated node
				const int co = bf.node_id();

				// Add to local matrix
				d(vSegmentFunction[fct].loc_fct, co) -= val;
			}
		}
	}

	// we're done
	return true;
}

} // namespace ug


#endif /*__H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__NEUMANN_BOUNDARY__FV1__NEUMANN_BOUNDARY_IMPL__*/
