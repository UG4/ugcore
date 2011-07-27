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
#include "common/util/provider.h"

namespace ug{

template <typename TBaseElem, typename TRefElem>
struct OrientationOffset{
	static bool get(std::vector<size_t>& vOrientOffset, int p,
	                const TRefElem& rRefElem,
	                size_t nrObj, std::vector<VertexBase*> vCorner)
	{
		throw(UGFatalError("GetOrientationOffset not implemented for that type"
				"of base element."));
	}
};

template <typename TRefElem>
struct OrientationOffset<EdgeBase, TRefElem>{
	static bool get(std::vector<size_t>& vOrientOffset, int p,
	                const TRefElem& rRefElem,
	                size_t nrEdge, std::vector<VertexBase*> vCorner)
	{
	//	compare the two corners of the edge
		const size_t co0 = rRefElem.id(1, nrEdge, 0, 0);
		const size_t co1 = rRefElem.id(1, nrEdge, 0, 1);

	//	resize the orientation array
		vOrientOffset.resize(p-1);

	//	the standard orientation is from co0 -> co1.
	//	we define to store the dofs that way if co0 has smaller address than
	//	co1. If this is not the case now, we have to invert the order.
		if(vCorner[co0] < vCorner[co1])
		{
			for(int i = 0; i <= p-2; ++i)
				vOrientOffset[i] = i;
		}
		else
		{
			for(int i = 0; i <= p-2; ++i)
				vOrientOffset[i] = (p-2)-i;
		}

	//	done
		return true;
	}
};



////////////////////////////////////////////////////////////////////////////////
//	Local Indices
////////////////////////////////////////////////////////////////////////////////

template<typename TElem, typename TBaseElem>
void DoFDistribution::indices(TElem* elem, LocalIndices& ind,
                              const std::vector<TBaseElem*>& vElem, size_t numNatural,
                              bool bHang, const std::vector<VertexBase*>& vCorner) const
{
//	get reference element type
	typedef typename reference_element_traits<TElem>::reference_element_type
			ref_elem_type;

//	get dimension
	static const int d = geometry_traits<TBaseElem>::BASE_OBJECT_TYPE_ID;

//	get roid type
	static const ReferenceObjectID roid = ref_elem_type::REFERENCE_OBJECT_ID;

//	reference element
	static const ref_elem_type& rRef = Provider::get<ref_elem_type>();

//	add normal dofs
	UG_ASSERT(vElem.size() == numNatural, "Size must match");
	for(size_t i = 0; i < numNatural; ++i)
	{
	//	get subelement
		TBaseElem* subElem = vElem[i];

	//	get subset index
		const int si = m_pISubsetHandler->get_subset_index(subElem);
		UG_ASSERT(si >= 0, "Invalid subset index " << si);

	//	read algebra index
		const size_t firstIndex = first_index(subElem);

	//	get reference object id for subselement
		const int subRoid = subElem->reference_object_id();

	//	loop all functions
		for(size_t fct = 0; fct < num_fct(); ++fct)
		{
		//	check if function is defined on the subset
			if(!is_def_in_subset(fct, si)) continue;

		//	get local shape function id
			LFEID lsfsID = local_finite_element_id(fct);

		//	get trial space
			const ILocalDoFSet& lds = LocalDoFSetProvider::get(roid, lsfsID);

		//	get number of DoFs in this sub-geometric object
			const size_t numDoFsOnSub = lds.num_dof(d, i);

		//	vector storing the computed offsets. If in correct order,
		//	this would be: [0, 1, 2, ...]. But this is usually not the
		// 	case and the numbers 0 to numDoFsOnSub-1 are permuted
			std::vector<size_t> vOrientOffset(numDoFsOnSub);

		//	check if on subobject or if there are no more than one dof on
		//	the subobject. If not, just add the dofs, no orientation needed.
			if(d < m_maxDimWithDoFs && numDoFsOnSub > 1)
			{
				//	get the orientation for this
					OrientationOffset<TBaseElem, ref_elem_type>::get
						(vOrientOffset, lsfsID.order(), rRef, i, vCorner);

					UG_ASSERT(vOrientOffset.size() == numDoFsOnSub, "Number of"
					          " offsets "<<vOrientOffset.size()<<" != Number of "
					          "DoFs on subelement "<<numDoFsOnSub);
			}

			if(!m_bGrouped)
			{
			//	compute index
				const size_t index = firstIndex + m_vvvOffsets[subRoid][si][fct];

			//	check if on subobject or if there are no more than one dof on
			//	the subobject. If not, just add the dofs, no orientation needed.
				if(d >= m_maxDimWithDoFs || numDoFsOnSub <= 1)
				{
				//	add dof to local indices
					for(size_t j = 0; j < numDoFsOnSub; ++j)
						ind.push_back_index(fct, index+j);
				}
			//	if on a subobject and more than one dof on it, the orientation
			//	must be taken into account in order to ensure a continuous space
				else
				{
					for(size_t j = 0; j < numDoFsOnSub; ++j)
						ind.push_back_index(fct, vOrientOffset[j]);
				}
			}
			else
			{
			//	compute index
				const size_t index = firstIndex;
				const size_t comp = m_vvvOffsets[roid][si][fct];

			//	check if on subobject or if there are no more than one dof on
			//	the subobject. If not, just add the dofs, no orientation needed.
				if(d >= m_maxDimWithDoFs || numDoFsOnSub <= 1)
				{
				//	add dof to local indices
					for(size_t j = 0; j < numDoFsOnSub; ++j)
						ind.push_back_multi_index(fct, index, comp+j);
				}
			//	if on a subobject and more than one dof on it, the orientation
			//	must be taken into account in order to ensure a continuous space
				else
				{
					for(size_t j = 0; j < numDoFsOnSub; ++j)
						ind.push_back_multi_index(fct, index, comp+vOrientOffset[j]);
				}
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

//	Flag to indicate iff orientation is required. If this is the case, the
//	vertices must be extracted in any case, since the orientation is done via
//	the vertices of a subelement
//	\todo: this must not always be true, introduce finer decision
	bool bOrientRequired = true;

//	needed for orientation
	std::vector<VertexBase*> vCorner;

//	get all sub-elements and add indices on those
	if(m_vMaxDoFsInDim[VERTEX] > 0 || bOrientRequired == true)
		CollectVertices(vCorner, *grid, elem);

	if(m_vMaxDoFsInDim[VERTEX] > 0)
	{
		const size_t numNatural = ref_elem_type::num_corners;
		indices<TElem, VertexBase>(elem, ind, vCorner, numNatural, bHang, vCorner);
	}
	if(m_vMaxDoFsInDim[EDGE] > 0)
	{
		std::vector<EdgeBase*> vElem;
		CollectEdgesSorted(vElem, *grid, elem);
		const size_t numNatural = ref_elem_type::num_edges;
		indices<TElem, EdgeBase>(elem, ind, vElem, numNatural, bHang, vCorner);
	}
	if(m_vMaxDoFsInDim[FACE] > 0)
	{
		std::vector<Face*> vElem;
		CollectFacesSorted(vElem, *grid, elem);
		const size_t numNatural = ref_elem_type::num_faces;
		indices<TElem, Face>(elem, ind, vElem, numNatural, bHang, vCorner);
	}
	if(m_vMaxDoFsInDim[VOLUME] > 0)
	{
		std::vector<Volume*> vElem;
		CollectVolumes(vElem, *grid, elem);
		const size_t numNatural = ref_elem_type::num_volumes;
		indices<TElem, Volume>(elem, ind, vElem, numNatural, bHang, vCorner);
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

//	get reference element type
	typedef typename reference_element_traits<TElem>::reference_element_type
			ref_elem_type;

// 	Grid
	Grid* grid = m_pStorageManager->get_assigned_grid();

//	Flag to indicate iff orientation is required. If this is the case, the
//	vertices must be extracted in any case, since the orientation is done via
//	the vertices of a subelement
//	\todo: this must not always be true, introduce finer decision
	bool bOrientRequired = true;

//	needed for orientation
	std::vector<VertexBase*> vCorner;

//	get all sub-elements and add indices on those
	if(m_vMaxDoFsInDim[VERTEX] > 0 || bOrientRequired == true)
		CollectVertices(vCorner, *grid, elem);

	if(m_vMaxDoFsInDim[VERTEX] > 0)
	{
		const size_t numNatural = ref_elem_type::num_corners;
		multi_indices<TElem, VertexBase>(elem, fct, ind, vCorner, numNatural, vCorner);
	}
	if(m_vMaxDoFsInDim[EDGE] > 0)
	{
		std::vector<EdgeBase*> vElem;
		CollectEdgesSorted(vElem, *grid, elem);
		const size_t numNatural = ref_elem_type::num_edges;
		multi_indices<TElem, EdgeBase>(elem, fct, ind, vElem, numNatural, vCorner);
	}
	if(m_vMaxDoFsInDim[FACE] > 0)
	{
		std::vector<Face*> vElem;
		CollectFacesSorted(vElem, *grid, elem);
		const size_t numNatural = ref_elem_type::num_faces;
		multi_indices<TElem, Face>(elem, fct, ind, vElem, numNatural, vCorner);
	}
	if(m_vMaxDoFsInDim[VOLUME] > 0)
	{
		std::vector<Volume*> vElem;
		CollectVolumes(vElem, *grid, elem);
		const size_t numNatural = ref_elem_type::num_volumes;
		multi_indices<TElem, Volume>(elem, fct, ind, vElem, numNatural, vCorner);
	}

//	return number of indices
	return ind.size();
}

template<typename TElem, typename TBaseElem>
void DoFDistribution::multi_indices(TElem* elem, size_t fct, multi_index_vector_type& ind,
                                    const std::vector<TBaseElem*>& vElem, size_t numNatural,
                                    const std::vector<VertexBase*>& vCorner) const
{
//	get reference element type
	typedef typename reference_element_traits<TElem>::reference_element_type
			ref_elem_type;

//	get dimension
	static const int d = geometry_traits<TBaseElem>::BASE_OBJECT_TYPE_ID;

//	get roid type
	static const ReferenceObjectID roid = ref_elem_type::REFERENCE_OBJECT_ID;

//	reference element
	static const ref_elem_type& rRef = Provider::get<ref_elem_type>();

//	loop passed elements
	UG_ASSERT(vElem.size() == numNatural, "Size must match");
	for(size_t i = 0; i < numNatural; ++i)
	{
	//	get subelement
		TBaseElem* subElem = vElem[i];

	//	get subset index
		const int si = m_pISubsetHandler->get_subset_index(subElem);
		UG_ASSERT(si >= 0, "Invalid subset index " << si);

	//	read algebra index
		const size_t firstIndex = first_index(subElem);

	//	get reference object id for subselement
		const int subRoid = subElem->reference_object_id();

	//	check if function is defined on the subset
		if(!is_def_in_subset(fct, si)) continue;

	//	check if dof given
		if(m_vvNumDoFsOnROID[subRoid][si] == 0) continue;

	//	get local shape function id
		LFEID lsfsID = local_finite_element_id(fct);

	//	get trial space
		const ILocalDoFSet& lsfs = LocalDoFSetProvider::get(roid, lsfsID);

	//	get number of DoFs in this sub-geometric object
		const size_t numDoFsOnSub = lsfs.num_dof(d, i);

	//	vector storing the computed offsets. If in correct order,
	//	this would be: [0, 1, 2, ...]. But this is usually not the
	// 	case and the numbers 0 to numDoFsOnSub-1 are permuted
		std::vector<size_t> vOrientOffset;

	//	check if on subobject or if there are no more than one dof on
	//	the subobject. If not, just add the dofs, no orientation needed.
		if(d < m_maxDimWithDoFs && numDoFsOnSub > 1)
		{
			//	get the orientation for this
				OrientationOffset<TBaseElem, ref_elem_type>::get
					(vOrientOffset, lsfsID.order(), rRef, i, vCorner);

				UG_ASSERT(vOrientOffset.size() == numDoFsOnSub, "Number of"
						  " offsets "<<vOrientOffset.size()<<" != Number of "
						  "DoFs on subelement "<<numDoFsOnSub);
		}

		if(!m_bGrouped)
		{
		//	compute index
			const size_t index = firstIndex + m_vvvOffsets[subRoid][si][fct];

		//	add dof to local indices
			if(d >= m_maxDimWithDoFs || numDoFsOnSub <= 1)
			{
				for(size_t j = 0; j < numDoFsOnSub; ++j)
					ind.push_back(index_type(index+j,0));
			}
			else
			{
				for(size_t j = 0; j < numDoFsOnSub; ++j)
					ind.push_back(index_type(index+vOrientOffset[j],0));
			}
		}
		else
		{
		//	compute index
			const size_t comp = m_vvvOffsets[subRoid][si][fct];

		//	add dof to local indices
			if(d >= m_maxDimWithDoFs || numDoFsOnSub <= 1)
			{
				for(size_t j = 0; j < numDoFsOnSub; ++j)
					ind.push_back(index_type(firstIndex, comp+j));
			}
			else
			{
				for(size_t j = 0; j < numDoFsOnSub; ++j)
					ind.push_back(index_type(firstIndex, comp+vOrientOffset[j]));
			}
		}
	}

//	return number of indices
	return;
}

template<typename TElem>
size_t DoFDistribution::inner_multi_indices(TElem* elem, size_t fct,
                                            multi_index_vector_type& ind,
                                            bool bClear) const
{
//	clear if requested
	if(bClear) ind.clear();

//	get reference element type
	typedef typename reference_element_traits<TElem>::reference_element_type
			ref_elem_type;

//	get base object type
	typedef typename geometry_traits<TElem>::geometric_base_object geometric_base_object;

//	get dimension
	static const int d = geometry_traits<geometric_base_object>::BASE_OBJECT_TYPE_ID;

//	get roid type
	static const ReferenceObjectID roid = ref_elem_type::REFERENCE_OBJECT_ID;

//	reference element
	static const ref_elem_type& rRef = Provider::get<ref_elem_type>();

//	get subset index
	const int si = m_pISubsetHandler->get_subset_index(elem);
	UG_ASSERT(si >= 0, "Invalid subset index " << si);

//	read algebra index
	const size_t firstIndex = first_index(elem);

//	check if function is defined on the subset
	if(!is_def_in_subset(fct, si)) return ind.size();

//	check if dof given
	if(m_vvNumDoFsOnROID[roid][si] == 0) return ind.size();

//	get local shape function id
	LFEID lsfsID = local_finite_element_id(fct);

//	get trial space
	const ILocalDoFSet& lsfs = LocalDoFSetProvider::get(roid, lsfsID);

//	get number of DoFs in this sub-geometric object
	const size_t numDoFsOnSub = lsfs.num_dof(roid);

//	vector storing the computed offsets. If in correct order,
//	this would be: [0, 1, 2, ...]. But this is usually not the
// 	case and the numbers 0 to numDoFsOnSub-1 are permuted
	std::vector<size_t> vOrientOffset;

//	check if on subobject or if there are no more than one dof on
//	the subobject. If not, just add the dofs, no orientation needed.
	if(d < m_maxDimWithDoFs && numDoFsOnSub > 1)
	{
	// 	Grid
		Grid* grid = m_pStorageManager->get_assigned_grid();

	//	get corners
		std::vector<VertexBase*> vCorner;
		CollectVertices(vCorner, *grid, elem);

	//	get the orientation for this
		OrientationOffset<geometric_base_object, ref_elem_type>::get
			(vOrientOffset, lsfsID.order(), rRef, 0, vCorner);

		UG_ASSERT(vOrientOffset.size() == numDoFsOnSub, "Number of"
				  " offsets "<<vOrientOffset.size()<<" != Number of "
				  "DoFs on subelement "<<numDoFsOnSub);
	}

	if(!m_bGrouped)
	{
	//	compute index
		const size_t index = firstIndex + m_vvvOffsets[roid][si][fct];

	//	add dof to local indices
		if(d >= m_maxDimWithDoFs || numDoFsOnSub <= 1)
		{
			for(size_t j = 0; j < numDoFsOnSub; ++j)
				ind.push_back(index_type(index+j,0));
		}
		else
		{
			for(size_t j = 0; j < numDoFsOnSub; ++j)
				ind.push_back(index_type(index+vOrientOffset[j],0));
		}
	}
	else
	{
	//	compute index
		const size_t comp = m_vvvOffsets[roid][si][fct];

	//	add dof to local indices
		if(d >= m_maxDimWithDoFs || numDoFsOnSub <= 1)
		{
			for(size_t j = 0; j < numDoFsOnSub; ++j)
				ind.push_back(index_type(firstIndex, comp+j));
		}
		else
		{
			for(size_t j = 0; j < numDoFsOnSub; ++j)
				ind.push_back(index_type(firstIndex, comp+vOrientOffset[j]));
		}
	}

//	forward request to inner element
	return ind.size();
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
//	\todo: Do we really need them sorted ?!
	if(m_vMaxDoFsInDim[0] > 0)
	{
		std::vector<VertexBase*> vElem;
		CollectVertices(vElem, *grid, elem);
		algebra_indices<VertexBase>(vElem, ind);
	}
	if(m_vMaxDoFsInDim[1] > 0)
	{
		std::vector<EdgeBase*> vElem;
		CollectEdgesSorted(vElem, *grid, elem);
		algebra_indices<EdgeBase>(vElem, ind);
	}
	if(m_vMaxDoFsInDim[2] > 0)
	{
		std::vector<Face*> vElem;
		CollectFacesSorted(vElem, *grid, elem);
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
