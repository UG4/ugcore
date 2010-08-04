/*
 * p1conform_impl.h
 *
 *  Created on: 13.07.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__DOF_MANAGER__P1CONFORM_IMPL__
#define __H__LIB_DISCRETIZATION__DOF_MANAGER__P1CONFORM_IMPL__

#include <vector>

#include "./p1conform.h"

namespace ug{

///////////////////////////////////////
// P1ConformDoFDistribution
///////////////////////////////////////

///////////// LocalIndex access /////////////////

template<typename TElem>
void
P1ConformDoFDistribution::
update_indices(TElem* elem, LocalIndices& ind) const
{
	const ReferenceObjectID refID = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	const ReferenceElement& refElem = ReferenceElementFactory::get_reference_element(refID);

	for(size_t i = 0; i <  refElem.num_obj(0); ++i)
	{
		VertexBase* vrt = get_vertex(elem, i);
		int si = m_pISubsetHandler->get_subset_index(vrt);

		const size_t index = m_pStorageManager->m_vSubsetInfo[si].aaDoFVRT[vrt];
		for(size_t fct = 0; fct < ind.num_fct(); ++fct)
		{
			ind.set_index(i + fct*refElem.num_obj(0), index + ind.fct_id(fct));
		}
	}
}

template<typename TElem>
void
P1ConformDoFDistribution::
update_inner_indices(TElem* elem, LocalIndices& ind) const
{
	const ReferenceObjectID refID = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	if(refID != ROID_VERTEX) return;
	else return update_indices(elem, ind);
}

///////////// Multi index access /////////////////

template<typename TElem>
size_t
P1ConformDoFDistribution::
num_multi_indices(TElem* elem, size_t fct) const
{
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
	return ref_elem_type::num_corners;
}

template<typename TElem>
size_t
P1ConformDoFDistribution::
num_inner_multi_indices(TElem* elem, size_t fct) const
{
	if(geometry_traits<TElem>::REFERENCE_OBJECT_ID == ROID_VERTEX) return 1;
	else return 0;
}

template<typename TElem>
size_t
P1ConformDoFDistribution::
get_multi_indices(TElem* elem, size_t fct, multi_index_vector_type& ind) const
{
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;

	const size_t numDofs = ref_elem_type::num_corners;
	ind.resize(numDofs);
	for(size_t i = 0; i < numDofs; ++i)
	{
		VertexBase* vrt = get_vertex(elem, i);
		int si = m_pISubsetHandler->get_subset_index(vrt);

		ind[i][0] = m_pStorageManager->m_vSubsetInfo[si].aaDoFVRT[vrt] + fct;
		ind[i][1] = 0;
	}
	return numDofs;
}

template<typename TElem>
size_t
P1ConformDoFDistribution::
get_inner_multi_indices(TElem* elem, size_t fct, multi_index_vector_type& ind) const
{
	if(geometry_traits<TElem>::REFERENCE_OBJECT_ID == ROID_VERTEX)
	{
		VertexBase* vrt = get_vertex(elem, 0);
		int si = m_pISubsetHandler->get_subset_index(vrt);
		ind.resize(1);
		ind[0][0] = m_pStorageManager->m_vSubsetInfo[si].aaDoFVRT[vrt] + fct;
		ind[0][1] = 0;
		return 1;
	}
	else return 0;
}

///////////// Algebra index access /////////////////

template<typename TElem>
size_t
P1ConformDoFDistribution::
num_algebra_indices(TElem* elem, size_t fct) const
{
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
	return ref_elem_type::num_corners;
}

template<typename TElem>
size_t
P1ConformDoFDistribution::
num_inner_algebra_indices(TElem* elem, size_t fct) const
{
	if(geometry_traits<TElem>::REFERENCE_OBJECT_ID == ROID_VERTEX) return 1;
	else return 0;
}

template<typename TElem>
void
P1ConformDoFDistribution::
get_algebra_indices(TElem* elem, algebra_index_vector_type& ind) const
{
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;

	ind.clear();

	const int elem_si = m_pISubsetHandler->get_subset_index(elem);
	for(size_t fct = 0; fct < num_fct(elem_si); ++fct)
	{
		for(size_t i = 0; i < ref_elem_type::num_corners; ++i)
		{
			VertexBase* vrt = get_vertex(elem, i);
			int si = m_pISubsetHandler->get_subset_index(vrt);
			const size_t index = m_pStorageManager->m_vSubsetInfo[si].aaDoFVRT[vrt] + fct;
			ind.push_back(index);
		}
	}
}

template<typename TElem>
void
P1ConformDoFDistribution::
get_inner_algebra_indices(TElem* elem, algebra_index_vector_type& ind) const
{
	ind.clear();
	if(geometry_traits<TElem>::REFERENCE_OBJECT_ID == ROID_VERTEX)
	{
		VertexBase* vrt = get_vertex(elem, 0);
		int si = m_pISubsetHandler->get_subset_index(vrt);

		const size_t index =  m_pStorageManager->m_vSubsetInfo[si].aaDoFVRT[vrt];
		for(size_t fct = 0; fct < num_fct(si); ++fct)
		{
			ind.push_back(index + fct);
		}
	}
}


///////////////////////////////////////
// GroupedP1ConformDoFDistribution
///////////////////////////////////////


///////////// LocalIndex access /////////////////

template<typename TElem>
void
GroupedP1ConformDoFDistribution::
update_indices(TElem* elem, LocalIndices& ind) const
{
	const ReferenceObjectID refID = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	const ReferenceElement& refElem = ReferenceElementFactory::get_reference_element(refID);

	for(size_t i = 0; i <  refElem.num_obj(0); ++i)
	{
		VertexBase* vrt = elem->vertex(i);
		int si = m_pISubsetHandler->get_subset_index(vrt);

		const size_t index = m_pStorageManager->m_vSubsetInfo[si].aaDoFVRT[vrt];
		ind.set_index(i, index);
	}
}

template<typename TElem>
void
GroupedP1ConformDoFDistribution::
update_inner_indices(TElem* elem, LocalIndices& ind) const
{
	const ReferenceObjectID refID = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	if(refID != ROID_VERTEX) return;
	else return update_indices(elem, ind);
}


///////////// Multi Index access /////////////////

template<typename TElem>
size_t
GroupedP1ConformDoFDistribution::
num_multi_indices(TElem* obj, size_t fct) const
{
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
	return ref_elem_type::num_corners;
}

template<typename TElem>
size_t
GroupedP1ConformDoFDistribution::
num_inner_multi_indices(TElem* obj, size_t fct) const
{
	if(geometry_traits<TElem>::REFERENCE_OBJECT_ID == ROID_VERTEX) return 1;
	else return 0;
}

template<typename TElem>
size_t
GroupedP1ConformDoFDistribution::
get_multi_indices(TElem* elem, size_t fct, multi_index_vector_type& ind) const
{
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;

	const size_t numDofs = ref_elem_type::num_corners;
	ind.resize(numDofs);
	for(size_t i = 0; i < numDofs; ++i)
	{
		VertexBase* vrt = get_vertex(elem, i);
		int si = m_pISubsetHandler->get_subset_index(vrt);

		ind[i][0] = m_pStorageManager->m_vSubsetInfo[si].aaDoFVRT[vrt];
		ind[i][1] = fct;
	}
	return numDofs;
}

template<typename TElem>
size_t
GroupedP1ConformDoFDistribution::
get_inner_multi_indices(TElem* obj, size_t fct, multi_index_vector_type& ind) const
{
	if(geometry_traits<TElem>::REFERENCE_OBJECT_ID == ROID_VERTEX)
	{
		VertexBase* vrt = get_vertex(obj, 0);
		int si = m_pISubsetHandler->get_subset_index(vrt);
		ind.resize(1);
		ind[0][0] = m_pStorageManager->m_vSubsetInfo[si].aaDoFVRT[vrt];
		ind[0][1] = fct;
		return 1;
	}
	else
	{
		ind.clear();
		return 0;
	}
}


///////////// Algebra Index access /////////////////

/// number of algebra indices on element for a function (Element + Closure of Element)
template<typename TElem>
size_t
GroupedP1ConformDoFDistribution::
num_algebra_indices(TElem* elem, size_t fct) const
{
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
	return ref_elem_type::num_corners;
}

/// number of algebras indices on element for a function (only inner part of Element)
template<typename TElem>
size_t
GroupedP1ConformDoFDistribution::
num_inner_algebra_indices(TElem* elem, size_t fct) const
{
	if(geometry_traits<TElem>::REFERENCE_OBJECT_ID == ROID_VERTEX) return 1;
	else return 0;
}

template<typename TElem>
void
GroupedP1ConformDoFDistribution::
get_algebra_indices(TElem* elem, algebra_index_vector_type& ind) const
{
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;

	ind.clear();

	// if no functions, return
	const int elem_si = m_pISubsetHandler->get_subset_index(elem);
	if(num_fct(elem_si) == 0) return;

	for(size_t i = 0; i < ref_elem_type::num_corners; ++i)
	{
			VertexBase* vrt = elem->vertex(i);
			int si = m_pISubsetHandler->get_subset_index(vrt);
			const size_t index = m_pStorageManager->m_vSubsetInfo[si].aaDoFVRT[vrt];
			ind.push_back(index);
	}
}

template<typename TElem>
void
GroupedP1ConformDoFDistribution::
get_inner_algebra_indices(TElem* obj, algebra_index_vector_type& ind) const
{
	ind.clear();
	if((GeometricBaseObject)geometry_traits<TElem>::BASE_OBJECT_TYPE_ID == (GeometricBaseObject)VERTEX)
	{
		VertexBase* vrt = (VertexBase*)obj;
		int si = m_pISubsetHandler->get_subset_index(vrt);
		if(num_fct(si) == 0) return;

		const size_t index =  m_pStorageManager->m_vSubsetInfo[si].aaDoFVRT[vrt];
		ind.push_back(index);
	}
}


} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__DOF_MANAGER__P1CONFORM_IMPL__ */
