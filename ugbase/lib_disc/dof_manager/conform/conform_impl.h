/*
 * conform_impl.h
 *
 *  Created on: 12.07.2011
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__DOF_MANAGER__CONFORM_IMPL__
#define __H__UG__LIB_DISC__DOF_MANAGER__CONFORM_IMPL__

#include <vector>

#include "conform.h"
#include "lib_disc/local_finite_element/local_dof_set.h"
#include "common/util/provider.h"

namespace ug{

template <typename TBaseElem, typename TRefElem>
struct OrientationOffset{
	static bool get(std::vector<size_t>& vOrientOffset, const size_t p,
	                const TRefElem& rRefElem,
	                const size_t nrObj, const size_t numDoFsOnSub,
	                const std::vector<VertexBase*>& vCorner)
	{
		throw(UGFatalError("GetOrientationOffset not implemented for that type"
				"of base element."));
	}
};

/*
 * Orientation of an edge:
 * If DoFs are assigned to a lower-dimensional edge and we have a degree higher
 * than 2 (i.e. more than one DoF on the edge) orientation is required to
 * ensure continuity of the shape functions. This means, that each element
 * that has the edge as a subelement, must number the dofs on the edge equally
 * in global numbering.
 *
 * The idea is as follows: We induce a global ordering of dofs on the edge by
 * using the storage position of the vertices of the edge. By our own definition
 * we say, that dofs are always assigned in a line from the vertices with lower
 * storage position to the vertices with the higher one. Now, in the local ordering
 * of dofs on the reference element, the edge may have been a different numbering
 * for the corners.
 * Thus, we have to distinguish two case:
 * a) corner 0 in reference element numbering has the smaller storage position: in
 * 		this case we can simply use the usual offset numbering
 * b) corner 0 in reference element numbering has a bigger storage position: in
 * 		this case we have to use the reverse order as offset numbering
 */
template <typename TRefElem>
struct OrientationOffset<EdgeBase, TRefElem>{
	static bool get(std::vector<size_t>& vOrientOffset, const size_t p,
	                const TRefElem& rRefElem,
	                const size_t nrEdge, const size_t numDoFsOnSub,
	                const std::vector<VertexBase*>& vCorner)
	{
	//	check
		UG_ASSERT(p-1 == numDoFsOnSub, "Wrong number of dofs on sub");
		UG_ASSERT(p > 2, "Orientation only needed for p > 2, but given p="<<p);

	//	compare the two corners of the edge
		const size_t co0 = rRefElem.id(1, nrEdge, 0, 0);
		const size_t co1 = rRefElem.id(1, nrEdge, 0, 1);

	//	resize the orientation array
		vOrientOffset.resize(numDoFsOnSub);

	//	the standard orientation is from co0 -> co1.
	//	we define to store the dofs that way if co0 has smaller address than
	//	co1. If this is not the case now, we have to invert the order.
		if(vCorner[co0] < vCorner[co1])
		{
			for(size_t i = 0; i < numDoFsOnSub; ++i)
				vOrientOffset[i] = i;
		}
		else
		{
			for(size_t i = 0; i < numDoFsOnSub; ++i)
				vOrientOffset[i] = (p-2)-i;
		}

	//	done
		return true;
	}
};

/*
 * Orientation of a Face:
 * If DoFs are assigned to a lower-dimensional face and we have a degree higher
 * than 2 (i.e. more than one DoF on the face) orientation is required to
 * ensure continuity of the shape functions. This means, that each element
 * that has the face as a subelement, must number the dofs on the face equally
 * in global numbering.
 *
 * The idea of the ordering is as follows:
 * We define that DoFs are always assigned to the face in a prescribed order.
 * In ordr to do so, we find the vertex of the face with the smallest storage
 * position. Then we find the vertex, that is connected to the smallest vertex
 * by an edge, with the (second) smallest storage position. Thus, we have a
 * situation like this:
 *
 *  	*						*--------*			^ j
 *  	|  \					|		 |			|
 *  	|    \				    |	     |			|
 * 		|      \				|		 |			|-----> i
 *  	*------ *				*--------*
 *  smallest   second		smallest	second
 *
 * In this picture all rotations and mirroring can appear. We define that the
 * DoFs on the face are always numbered in x direction first, continuing in the
 * next row in y direction, numbering in x again and continuing in y, etc.
 * E.g. this gives for a p = 4 element (showing only inner dofs):
 *
 *  	*						*-------*			^ j
 *  	|5 \					| 6	7 8	|			|
 *  	|3 4 \					| 3 4 5	|			|
 * 		|0 1 2 \				| 0 1 2 |			|-----> i
 *  	*-------*				*-------*
 *  smallest   second		smallest	second
 *
 * Now, the face vertices of the given element have a local numbering vCo[] in
 * the reference element. The DoFs in the reference element meaning are ordered
 * as if vCo[0] was smallest and vCo[1] was second. Thus, now in real world,
 * these must not match and we have to find the orientation of the face. Smallest
 * and second is computed and than a mapping is set up.
 */
template <typename TRefElem>
struct OrientationOffset<Face, TRefElem>{
	static bool get(std::vector<size_t>& vOrientOffset, const size_t p,
	                const TRefElem& rRefElem,
	                const size_t nrFace, const size_t numDoFsOnSub,
	                const std::vector<VertexBase*>& vCorner)
	{
	//	should only be called for p > 2, since else no orientation needed
		UG_ASSERT(p > 2, "Orientation only needed for p > 2, but given p="<<p);

	//	resize array
		vOrientOffset.resize(numDoFsOnSub);

	//	get corners of face
		const int numCo = rRefElem.num(2, nrFace, 0);
		std::vector<size_t> vCo(numCo);
		for(int i = 0; i < numCo; ++i)
			vCo[i] = rRefElem.id(2, nrFace, 0, i);

	//	find smallest
		int smallest = 0;
		for(int i = 1; i < numCo; ++i)
			if(vCorner[ vCo[i] ] < vCorner[ vCo[smallest] ]) smallest = i;

	//	find second smallest
		int second = (smallest+numCo-1)%numCo;
		const int next = (smallest+numCo+1)%numCo;
		if(vCorner[ vCo[next] ] < vCorner[ vCo[second] ]) second = next;

	//	map the i,j
		size_t map_i, map_j;

	//	in the inner, the number of dofs is as if it would be an element
	//	of degree p-2. We cache this here
		const size_t pInner = p-2;

	//	loop 'y'-direction
		size_t index = 0;
		for(size_t j = 0; j <= pInner; ++j)
		{
		//	for a quadrilateral we have a quadratic loop, but for a
		//	triangle we need to stop at the diagonal
			const size_t off = ((vCo.size() == 3) ? j : 0);

		//	loop 'x'-direction
			for(size_t i = 0; i <= pInner-off; ++i)
			{
			//	map the index-pair (i,j)
				map_face(map_i, map_j, i, j, numCo, pInner, smallest, second);

			//	linearize index and mapped index
				const size_t mappedIndex = (p-1) * map_j + map_i;

			//	check
				UG_ASSERT(mappedIndex < numDoFsOnSub, "Wrong mapped index");
				UG_ASSERT(index < numDoFsOnSub, "Wrong index");

			//	set mapping
				vOrientOffset[index++] = mappedIndex;
			}
		}

	//	case not found
		return false;
	}

	static void map_face(size_t& map_i, size_t& map_j,
	                     const size_t i, const size_t j,
	                     const int numCo, const size_t p,
	                     const int smallest, const int second)
	{
		UG_ASSERT(i <= p, "Wrong index");
		UG_ASSERT(j <= p, "Wrong index");

	//	handle rotation
		switch(numCo)
		{
			case 3:
			{
				switch(smallest)
				{
					case 0: map_i = i; map_j = j; break;
					case 1: map_i = j; map_j = p-i-j; break;
					case 2: map_i = p-i-j; map_j = i; break;
					default: throw(UGFatalError("Corner not found."));
				}
				break;
			}
			case 4:
			{
				switch(smallest)
				{
					case 0: map_i = i; map_j = j; break;
					case 1: map_i = j; map_j = p-i; break;
					case 2: map_i = p-i; map_j = p-j; break;
					case 3: map_i = p-j; map_j = i; break;
					default: throw(UGFatalError("Corner not found."));
				}
				break;
			}
			default: throw(UGFatalError("Num Corners not supported."));
		}

	//	handle mirroring
		if(second == (smallest+numCo-1)%numCo)
		{
			const size_t h = map_i; map_i = map_j; map_j = h;
		}
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
	static const ref_elem_type& rRef = Provider<ref_elem_type>::get();

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
		//	get the orientation for this subelement
			if(d < m_maxDimWithDoFs && numDoFsOnSub > 1)
					OrientationOffset<TBaseElem, ref_elem_type>::get
						(vOrientOffset, lsfsID.order(), rRef, i, numDoFsOnSub, vCorner);

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
						ind.push_back_index(fct, index+vOrientOffset[j]);
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

//	reference dimension
	static const int dim = ref_elem_type::dim;

//	resize the number of functions
	ind.resize_fct(num_fct());
	for(size_t fct = 0; fct < num_fct(); ++fct) ind.clear_dof(fct);

// 	Grid
	Grid* grid = m_pStorageManager->grid();

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
	if(dim >= EDGE && m_vMaxDoFsInDim[EDGE] > 0)
	{
		std::vector<EdgeBase*> vElem;
		CollectEdgesSorted(vElem, *grid, elem);
		const size_t numNatural = ref_elem_type::num_edges;
		indices<TElem, EdgeBase>(elem, ind, vElem, numNatural, bHang, vCorner);
	}
	if(dim >= FACE && m_vMaxDoFsInDim[FACE] > 0)
	{
		std::vector<Face*> vElem;
		CollectFacesSorted(vElem, *grid, elem);
		const size_t numNatural = ref_elem_type::num_faces;
		indices<TElem, Face>(elem, ind, vElem, numNatural, bHang, vCorner);
	}
	if(dim >= VOLUME && m_vMaxDoFsInDim[VOLUME] > 0)
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

//	reference dimension
	static const int dim = ref_elem_type::dim;

// 	Grid
	Grid* grid = m_pStorageManager->grid();

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
		static const size_t numNatural = ref_elem_type::num_corners;
		UG_ASSERT(vCorner.size() == numNatural, "Wrong number of sub-vertices of a "<<
				  ref_elem_type::REFERENCE_OBJECT_ID <<": num collected="<<vCorner.size()
				  <<", numNatural="<<numNatural);
		multi_indices<TElem, VertexBase>(elem, fct, ind, vCorner, numNatural, vCorner);
	}
	if(dim >= EDGE && m_vMaxDoFsInDim[EDGE] > 0)
	{
		std::vector<EdgeBase*> vElem;
		CollectEdgesSorted(vElem, *grid, elem);
		static const size_t numNatural = ref_elem_type::num_edges;
		UG_ASSERT(vElem.size() == numNatural, "Wrong number of sub-edges of a "<<
				  ref_elem_type::REFERENCE_OBJECT_ID <<": num collected="<<vElem.size()
				  <<", numNatural="<<numNatural);
		multi_indices<TElem, EdgeBase>(elem, fct, ind, vElem, numNatural, vCorner);
	}
	if(dim >= FACE && m_vMaxDoFsInDim[FACE] > 0)
	{
		std::vector<Face*> vElem;
		CollectFacesSorted(vElem, *grid, elem);
		static const size_t numNatural = ref_elem_type::num_faces;
		UG_ASSERT(vElem.size() == numNatural, "Wrong number of sub-faces of a "<<
				  ref_elem_type::REFERENCE_OBJECT_ID <<": num collected="<<vElem.size()
				  <<", numNatural="<<numNatural);
		multi_indices<TElem, Face>(elem, fct, ind, vElem, numNatural, vCorner);
	}
	if(dim >= VOLUME && m_vMaxDoFsInDim[VOLUME] > 0)
	{
		std::vector<Volume*> vElem;
		CollectVolumes(vElem, *grid, elem);
		static const size_t numNatural = ref_elem_type::num_volumes;
		UG_ASSERT(vElem.size() == numNatural, "Wrong number of sub-volumes of a "<<
				  ref_elem_type::REFERENCE_OBJECT_ID <<": num collected="<<vElem.size()
				  <<", numNatural="<<numNatural);
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
	static const ref_elem_type& rRef = Provider<ref_elem_type>::get();

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
	//	get the orientation for this subelement
		if(d < m_maxDimWithDoFs && numDoFsOnSub > 1)
				OrientationOffset<TBaseElem, ref_elem_type>::get
					(vOrientOffset, lsfsID.order(), rRef, i, numDoFsOnSub, vCorner);

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
	static const ref_elem_type& rRef = Provider<ref_elem_type>::get();

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
		Grid* grid = m_pStorageManager->grid();

	//	get corners
		std::vector<VertexBase*> vCorner;
		CollectVertices(vCorner, *grid, elem);

	//	get the orientation for this
		OrientationOffset<geometric_base_object, ref_elem_type>::get
			(vOrientOffset, lsfsID.order(), rRef, 0, numDoFsOnSub, vCorner);
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

//	get reference element type
	typedef typename reference_element_traits<TElem>::reference_element_type
			ref_elem_type;

//	reference dimension
	static const int dim = ref_elem_type::dim;

// 	Grid
	Grid* grid = m_pStorageManager->grid();

//	get all sub-elements and add indices on those
//	\todo: Do we really need them sorted ?!
	if(m_vMaxDoFsInDim[VERTEX] > 0)
	{
		std::vector<VertexBase*> vElem;
		CollectVertices(vElem, *grid, elem);
		algebra_indices<VertexBase>(vElem, ind);
	}
	if(dim >= EDGE && m_vMaxDoFsInDim[EDGE] > 0)
	{
		std::vector<EdgeBase*> vElem;
		CollectEdgesSorted(vElem, *grid, elem);
		algebra_indices<EdgeBase>(vElem, ind);
	}
	if(dim >= FACE && m_vMaxDoFsInDim[FACE] > 0)
	{
		std::vector<Face*> vElem;
		CollectFacesSorted(vElem, *grid, elem);
		algebra_indices<Face>(vElem, ind);
	}
	if(dim >= VOLUME && m_vMaxDoFsInDim[VOLUME] > 0)
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

#endif /* __H__UG__LIB_DISC__DOF_MANAGER__CONFORM_IMPL__ */
