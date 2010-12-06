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
update_indices(TElem* elem, LocalIndices& ind, bool withHanging) const
{
	const ReferenceObjectID refID = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	const ReferenceElement& refElem =
			ReferenceElementFactory::get_reference_element(refID);

	if(!withHanging)
	{
		for(size_t i = 0; i <  refElem.num_obj(0); ++i)
		{
			VertexBase* vrt = GetVertex(elem, i);
			int si = m_pISubsetHandler->get_subset_index(vrt);

			const size_t index = m_pStorageManager->m_vSubsetInfo[si].aaDoFVRT[vrt];

			size_t numFct = 0;
			for(size_t fct = 0; fct < ind.num_fct(); ++fct)
			{
				if(!is_def_in_subset(ind.fct_id(fct), si)) continue;
				ind.set_index(i + numFct*refElem.num_obj(0),
				              index + m_vvOffsets[si][ind.fct_id(fct)]);
				numFct++;
			}
		}
	}
	else
	{
		ind.clear();
		size_t algDof = 0;

		// natural dofs
		for(size_t i = 0; i < refElem.num_obj(0); ++i)
		{
			VertexBase* vrt = GetVertex(elem, i);
			int si = m_pISubsetHandler->get_subset_index(vrt);
			const size_t index = m_pStorageManager->m_vSubsetInfo[si].aaDoFVRT[vrt];

			for(size_t fct = 0; fct < ind.num_fct(); ++fct)
			{
				if(!is_def_in_subset(ind.fct_id(fct), si)) continue;

				ind.set_num_indices(algDof+1);
				ind.set_index(algDof, index + m_vvOffsets[si][ind.fct_id(fct)]);

				LocalIndices::multi_index_type dof_ind;
				dof_ind[0] = algDof;
				dof_ind[1] = 0;
				ind.add_dof(fct, dof_ind);

				algDof++;
			}
		}

		// hanging dofs

		// get natural edges
		{
			std::vector<EdgeBase*> vEdges;
			EdgeBase* ed = dynamic_cast<EdgeBase*>(elem);
			if(ed != NULL)
				CollectEdgesSorted(vEdges,
				                   *(m_pISubsetHandler->get_assigned_grid()), ed);

			Face* face = dynamic_cast<Face*>(elem);
			if(face != NULL)
				CollectEdgesSorted(vEdges,
				                   *(m_pISubsetHandler->get_assigned_grid()), face);

			Volume* vol = dynamic_cast<Volume*>(elem);
			if(vol != NULL)
				CollectEdgesSorted(vEdges,
				                   *(m_pISubsetHandler->get_assigned_grid()), vol);

			for(size_t i = 0; i < vEdges.size(); ++i)
			{
				ConstrainingEdge* edge = dynamic_cast<ConstrainingEdge*>(vEdges[i]);
				if(edge == NULL) continue;
				for(VertexBaseIterator iter = edge->constrained_vertices_begin();
						iter != edge->constrained_vertices_end(); ++iter)
				{
					VertexBase* vrt = *iter;
					int si = m_pISubsetHandler->get_subset_index(vrt);
					const size_t index =
							m_pStorageManager->m_vSubsetInfo[si].aaDoFVRT[vrt];

					for(size_t fct = 0; fct < ind.num_fct(); ++fct)
					{
						if(!is_def_in_subset(ind.fct_id(fct), si)) continue;

						ind.set_num_indices(algDof+1);
						ind.set_index(algDof, index + m_vvOffsets[si][ind.fct_id(fct)]);

						LocalIndices::multi_index_type dof_ind;
						dof_ind[0] = algDof;
						dof_ind[1] = 0;
						ind.add_dof(fct, dof_ind);

						algDof++;
					}
				}
			}
		}
		// get natural faces
		{
			std::vector<Face*> vFaces; vFaces.clear();
			Face* face = dynamic_cast<Face*>(elem);
			if(face != NULL)
				CollectFacesSorted(vFaces,
				                   *(m_pISubsetHandler->get_assigned_grid()), face);
			Volume* vol = dynamic_cast<Volume*>(elem);
			if(vol != NULL)
				CollectFacesSorted(vFaces,
				                   *(m_pISubsetHandler->get_assigned_grid()), vol);

			for(size_t i = 0; i < vFaces.size(); ++i)
			{
				ConstrainingQuadrilateral* quad =
						dynamic_cast<ConstrainingQuadrilateral*>(vFaces[i]);
				if(quad == NULL) continue;
				for(VertexBaseIterator iter = quad->constrained_vertices_begin();
						iter != quad->constrained_vertices_end(); ++iter)
				{
					VertexBase* vrt = *iter;
					int si = m_pISubsetHandler->get_subset_index(vrt);
					const size_t index =
							m_pStorageManager->m_vSubsetInfo[si].aaDoFVRT[vrt];

					for(size_t fct = 0; fct < ind.num_fct(); ++fct)
					{
						if(!is_def_in_subset(ind.fct_id(fct), si)) continue;

						ind.set_num_indices(algDof+1);
						ind.set_index(algDof, index + m_vvOffsets[si][ind.fct_id(fct)]);

						LocalIndices::multi_index_type dof_ind;
						dof_ind[0] = algDof;
						dof_ind[1] = 0;
						ind.add_dof(fct, dof_ind);

						algDof++;
					}
				}
			}
		}
	}
	return;
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
		VertexBase* vrt = GetVertex(elem, i);
		int si = m_pISubsetHandler->get_subset_index(vrt);

		ind[i][0] = m_pStorageManager->m_vSubsetInfo[si].aaDoFVRT[vrt]
		                                              + m_vvOffsets[si][fct];
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
		VertexBase* vrt = GetVertex(elem, 0);
		int si = m_pISubsetHandler->get_subset_index(vrt);
		ind.resize(1);
		ind[0][0] = m_pStorageManager->m_vSubsetInfo[si].aaDoFVRT[vrt]
		                                                + m_vvOffsets[si][fct];
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
	typedef typename reference_element_traits<TElem>::reference_element_type
			ref_elem_type;
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
	typedef typename reference_element_traits<TElem>::reference_element_type
			ref_elem_type;

	ind.clear();

	const int elem_si = m_pISubsetHandler->get_subset_index(elem);
	for(size_t fct = 0; fct < num_fct(); ++fct)
	{
		for(size_t i = 0; i < ref_elem_type::num_corners; ++i)
		{
			VertexBase* vrt = GetVertex(elem, i);
			int si = m_pISubsetHandler->get_subset_index(vrt);
			if(!is_def_in_subset(fct, si)) continue;
			const size_t index =
					m_pStorageManager->m_vSubsetInfo[si].aaDoFVRT[vrt]
					                            + m_vvOffsets[si][fct];
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
		VertexBase* vrt = GetVertex(elem, 0);
		int si = m_pISubsetHandler->get_subset_index(vrt);

		const size_t index =  m_pStorageManager->m_vSubsetInfo[si].aaDoFVRT[vrt];
		for(size_t fct = 0; fct < num_fct(); ++fct)
		{
			if(!is_def_in_subset(fct, si)) continue;
			ind.push_back(index + m_vvOffsets[si][fct]);
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
update_indices(TElem* elem, LocalIndices& ind, bool withHanging) const
{
	const ReferenceObjectID refID = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	const ReferenceElement& refElem =
			ReferenceElementFactory::get_reference_element(refID);

	if(withHanging) throw(UGError("Not implemented"));

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
	typedef typename reference_element_traits<TElem>::reference_element_type
			ref_elem_type;
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
	typedef typename reference_element_traits<TElem>::reference_element_type
			ref_elem_type;

	const size_t numDofs = ref_elem_type::num_corners;
	ind.resize(numDofs);
	for(size_t i = 0; i < numDofs; ++i)
	{
		VertexBase* vrt = GetVertex(elem, i);
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
		VertexBase* vrt = GetVertex(obj, 0);
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
	typedef typename reference_element_traits<TElem>::reference_element_type
			ref_elem_type;
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
	typedef typename reference_element_traits<TElem>::reference_element_type
			ref_elem_type;

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
	if((GeometricBaseObject)geometry_traits<TElem>::BASE_OBJECT_TYPE_ID
			== (GeometricBaseObject)VERTEX)
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
