/*
 * conform_impl.h
 *
 *  Created on: 12.07.2011
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__DOF_MANAGER__CONFORM_IMPL__
#define __H__LIB_DISCRETIZATION__DOF_MANAGER__CONFORM_IMPL__

#include <vector>

#include "conform.h"
#include "lib_discretization/local_finite_element/local_dof_set.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
//	Local Indices
////////////////////////////////////////////////////////////////////////////////

template<typename TElem, typename TBaseElem>
void DoFDistribution::indices(TElem* elem, LocalIndices& ind,
                              std::vector<TBaseElem*> vElem, size_t numNatural,
                              bool bHang) const
{
//	get reference element type
	typedef typename reference_element_traits<TElem>::reference_element_type
			ref_elem_type;

//	get dimension
	static const int d = geometry_traits<TBaseElem>::BASE_OBJECT_TYPE_ID;

//	get roid type
	static const ReferenceObjectID roid = ref_elem_type::REFERENCE_OBJECT_ID;

//	add normal dofs
	for(size_t i = 0; i < numNatural; ++i)
	{
	//	get vertex
		TBaseElem* subElem = vElem[i];

	//	get subset index
		const int si = m_pISubsetHandler->get_subset_index(subElem);
		UG_ASSERT(si >= 0, "Invalid subset index " << si);

	//	read algebra index
		const size_t firstIndex = first_index(subElem);

	//	get reference object id for subselement
		const int baseRoid = subElem->reference_object_id();

	//	loop all functions
		for(size_t fct = 0; fct < num_fct(); ++fct)
		{
		//	check if function is defined on the subset
			if(!is_def_in_subset(fct, si)) continue;

		//	get local shape function id
			LFEID lsfsID = local_finite_element_id(fct);

		//	get trial space
			const ILocalDoFSet& lds = LocalDoFSetProvider::get(lsfsID, roid);

		//	get number of DoFs in this sub-geometric object
			const size_t numDoFsOnSub = lds.num_dof(d, i);

			if(!m_bGrouped)
			{
			//	compute index
				const size_t index = firstIndex + m_vvvOffsets[baseRoid][si][fct];

			//	add dof to local indices
			// \TODO: ORIENTATION !!!
				for(size_t j = 0; j < numDoFsOnSub; ++j)
					ind.push_back_index(fct, index+j);
			}
			else
			{
			//	compute index
				const size_t index = firstIndex;
				const size_t comp = m_vvvOffsets[roid][si][fct];

			//	add dof to local indices
			// \TODO: ORIENTATION !!!
				for(size_t j = 0; j < numDoFsOnSub; ++j)
					ind.push_back_multi_index(fct, index, comp+j);
			}
		}
	}

}

template<typename TElem>
void DoFDistribution::indices(TElem* elem, LocalIndices& ind, bool bHang) const
{
//	get reference element type
	typedef typename reference_element_traits<TElem>::reference_element_type
			ref_elem_type;

//	resize the number of functions
	ind.resize_fct(num_fct());
	for(size_t fct = 0; fct < num_fct(); ++fct) ind.clear_dof(fct);

// 	Grid
	Grid* grid = m_pStorageManager->get_assigned_grid();

//	get all sub-elements and add indices on those
	if(m_vMaxDoFsInDim[VERTEX] > 0)
	{
		std::vector<VertexBase*> vElem;
		CollectVertices(vElem, *grid, elem);
		const size_t numNatural = ref_elem_type::num_corners;
		indices<TElem, VertexBase>(elem, ind, vElem, numNatural, bHang);
	}
	if(m_vMaxDoFsInDim[EDGE] > 0)
	{
		std::vector<EdgeBase*> vElem;
		CollectEdgesSorted(vElem, *grid, elem);
		const size_t numNatural = ref_elem_type::num_edges;
		indices<TElem, EdgeBase>(elem, ind, vElem, numNatural, bHang);
	}
	if(m_vMaxDoFsInDim[FACE] > 0)
	{
		std::vector<Face*> vElem;
		CollectFacesSorted(vElem, *grid, elem);
		const size_t numNatural = ref_elem_type::num_faces;
		indices<TElem, Face>(elem, ind, vElem, numNatural, bHang);
	}
	if(m_vMaxDoFsInDim[VOLUME] > 0)
	{
		std::vector<Volume*> vElem;
		CollectVolumes(vElem, *grid, elem);
		const size_t numNatural = ref_elem_type::num_volumes;
		indices<TElem, Volume>(elem, ind, vElem, numNatural, bHang);
	}

//	If no hanging dofs are required, we're done
	if(!bHang) return;

	UG_LOG("ERROR in 'DoFManager::indices': Hanging DoFs are currently not "
			"supported by this DoFManager. \n");
	throw(UGFatalError("Hanging DoFs not implemented."));

//	we're done
	return;
}

////////////////////////////////////////////////////////////////////////////////
//	Multi Indices
////////////////////////////////////////////////////////////////////////////////

template<typename TElem>
size_t DoFDistribution::multi_indices(TElem* elem, size_t fct,
                                      multi_index_vector_type& ind,
                                      bool bClear) const
{
//	clear indices
	if(bClear) ind.clear();

// 	Grid
	Grid* grid = m_pStorageManager->get_assigned_grid();

//	get all sub-elements and add indices on those
	if(m_vMaxDoFsInDim[0] > 0)
	{
		std::vector<VertexBase*> vElem;
		CollectVertices(vElem, *grid, elem);
		multi_indices<VertexBase>(vElem, fct, ind);
	}
	if(m_vMaxDoFsInDim[1] > 0)
	{
		std::vector<EdgeBase*> vElem;
		CollectEdges(vElem, *grid, elem);
		multi_indices<EdgeBase>(vElem, fct, ind);
	}
	if(m_vMaxDoFsInDim[2] > 0)
	{
		std::vector<Face*> vElem;
		CollectFaces(vElem, *grid, elem);
		multi_indices<Face>(vElem, fct, ind);
	}
	if(m_vMaxDoFsInDim[3] > 0)
	{
		std::vector<Volume*> vElem;
		CollectVolumes(vElem, *grid, elem);
		multi_indices<Volume>(vElem, fct, ind);
	}

//	return number of indices
	return ind.size();
}

template<typename TBaseElem>
size_t DoFDistribution::multi_indices(std::vector<TBaseElem*> vElem, size_t fct,
                                      multi_index_vector_type& ind) const
{
//	loop passed elements
	for(size_t i = 0; i < vElem.size(); ++i)
	{
	//	get reference object type
		const ReferenceObjectID roid = static_cast<ReferenceObjectID>(vElem[i]->reference_object_id());

	//	get subset index
		const int si = m_pISubsetHandler->get_subset_index(vElem[i]);
		UG_ASSERT(si >= 0, "Invalid subset index " << si);

	//	get first index
		const size_t firstIndex = first_index(vElem[i]);

	//	forward request to inner element
		inner_multi_indices(ind, firstIndex, si, fct, roid);
	}

//	return number of indices
	return ind.size();
}

template<typename TElem>
size_t DoFDistribution::inner_multi_indices(TElem* elem, size_t fct,
                                            multi_index_vector_type& ind,
                                            bool bClear) const
{
//	get reference element type
	typedef typename reference_element_traits<TElem>::reference_element_type
			ref_elem_type;

//	get roid type
	static const ReferenceObjectID roid = ref_elem_type::REFERENCE_OBJECT_ID;

//	get subset index
	const int si = m_pISubsetHandler->get_subset_index(elem);
	UG_ASSERT(si >= 0, "Invalid subset index " << si);

//	get first index
	const size_t firstIndex = first_index(elem);

//	clear if requested
	if(bClear) ind.clear();

//	forward request to inner element
	return inner_multi_indices(ind, firstIndex, si, fct, roid);
}

////////////////////////////////////////////////////////////////////////////////
//	Algebra Indices
////////////////////////////////////////////////////////////////////////////////

template<typename TElem>
size_t DoFDistribution::algebra_indices(TElem* elem,
                                        algebra_index_vector_type& ind,
                                        bool bClear) const
{
//	clear indices
	if(bClear) ind.clear();

// 	Grid
	Grid* grid = m_pStorageManager->get_assigned_grid();

//	get all sub-elements and add indices on those
	if(m_vMaxDoFsInDim[0] > 0)
	{
		std::vector<VertexBase*> vElem;
		CollectVertices(vElem, *grid, elem);
		algebra_indices<VertexBase>(vElem, ind);
	}
	if(m_vMaxDoFsInDim[1] > 0)
	{
		std::vector<EdgeBase*> vElem;
		CollectEdges(vElem, *grid, elem);
		algebra_indices<EdgeBase>(vElem, ind);
	}
	if(m_vMaxDoFsInDim[2] > 0)
	{
		std::vector<Face*> vElem;
		CollectFaces(vElem, *grid, elem);
		algebra_indices<Face>(vElem, ind);
	}
	if(m_vMaxDoFsInDim[3] > 0)
	{
		std::vector<Volume*> vElem;
		CollectVolumes(vElem, *grid, elem);
		algebra_indices<Volume>(vElem, ind);
	}

//	return number of indices
	return ind.size();
}

template<typename TBaseElem>
size_t DoFDistribution::algebra_indices(std::vector<TBaseElem*> vElem,
                                        algebra_index_vector_type& ind) const
{
//	loop passed elements
	for(size_t i = 0; i < vElem.size(); ++i)
	{
	//	get reference object type
		const ReferenceObjectID roid = static_cast<ReferenceObjectID>(vElem[i]->reference_object_id());

	//	get subset index
		const int si = m_pISubsetHandler->get_subset_index(vElem[i]);
		UG_ASSERT(si >= 0, "Invalid subset index " << si);

	//	get first index
		const size_t firstIndex = first_index(vElem[i]);

	//	forward request to inner element
		inner_algebra_indices(ind, firstIndex, si, roid);
	}

//	return number of indices
	return ind.size();
}

template<typename TElem>
size_t DoFDistribution::inner_algebra_indices(TElem* elem,
                                              algebra_index_vector_type& ind,
                                              bool bClear) const
{
//	get reference element type
	typedef typename reference_element_traits<TElem>::reference_element_type
			ref_elem_type;

//	get roid type
	static const ReferenceObjectID roid = ref_elem_type::REFERENCE_OBJECT_ID;

//	get subset index
	const int si = m_pISubsetHandler->get_subset_index(elem);
	UG_ASSERT(si >= 0, "Invalid subset index " << si);

//	get first algebra index
	const size_t firstIndex = first_index(elem);

//	clear indices
	if(bClear) ind.clear();

//	return
	return inner_algebra_indices(ind, firstIndex, si, roid);
}

} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__DOF_MANAGER__CONFORM_IMPL__ */
