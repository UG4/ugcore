/*
 * mg_dof_distribution.cpp
 *
 *  Created on: 06.12.2011
 *      Author: andreasvogel
 */

#include "common/log.h"
#include "mg_dof_distribution.h"
#include "mg_dof_distribution_impl.h"
#include "lib_disc/domain.h"
#include "lib_disc/local_finite_element/local_dof_set.h"
#include "lib_disc/reference_element/reference_element_util.h"
#include "lib_disc/common/groups_util.h"
#include "common/util/string_util.h"

using namespace std;

namespace ug{

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//	ComputeOrientationOffset
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

template <typename TBaseElem>
void ComputeOrientationOffset(std::vector<size_t>& vOrientOffset, const size_t p,
                              const ReferenceElement& rRefElem,
                              const size_t nrEdge, const size_t numDoFsOnSub,
                              const Grid::SecureVertexContainer& vCorner)
{}


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
template <>
void ComputeOrientationOffset<EdgeBase>
(std::vector<size_t>& vOrientOffset, const size_t p,
 const ReferenceElement& rRefElem,
 const size_t nrEdge, const size_t numDoFsOnSub,
 const Grid::SecureVertexContainer& vCorner)
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
};




static void MapFace(size_t& map_i, size_t& map_j,
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
				default: UG_THROW("Corner not found.");
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
				default: UG_THROW("Corner not found.");
			}
			break;
		}
		default: UG_THROW("Num Corners not supported.");
	}

//	handle mirroring
	if(second == (smallest+numCo-1)%numCo)
	{
		const size_t h = map_i; map_i = map_j; map_j = h;
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
template <>
void ComputeOrientationOffset<Face>
(std::vector<size_t>& vOrientOffset, const size_t p,
 const ReferenceElement& rRefElem,
 const size_t nrFace, const size_t numDoFsOnSub,
 const Grid::SecureVertexContainer& vCorner)
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
			MapFace(map_i, map_j, i, j, numCo, pInner, smallest, second);

		//	linearize index and mapped index
			const size_t mappedIndex = (p-1) * map_j + map_i;

		//	check
			UG_ASSERT(mappedIndex < numDoFsOnSub, "Wrong mapped index");
			UG_ASSERT(index < numDoFsOnSub, "Wrong index");

		//	set mapping
			vOrientOffset[index++] = mappedIndex;
		}
	}
};


////////////////////////////////////////////////////////////////////////////////
// MGDoFDistribution
////////////////////////////////////////////////////////////////////////////////

MGDoFDistribution::
MGDoFDistribution(SmartPtr<MultiGrid> spMG,
                  SmartPtr<MGSubsetHandler> spMGSH,
				  ConstSmartPtr<DoFDistributionInfo> spDDInfo,
                  bool bGrouped)
	: DoFDistributionInfoProvider(spDDInfo),
      m_bGrouped(bGrouped),
	  m_strictSubsetChecks(true),
	  m_parallelRedistributionMode(false),
	  m_spMG(spMG),
	  m_pMG(spMG.get()),
	  m_spMGSH(spMGSH)
{
	check_subsets();
	init_attachments();
	register_observer();
};

MGDoFDistribution::
~MGDoFDistribution()
{
	clear_attachments();
	unregister_observer();
}

void MGDoFDistribution::
begin_parallel_redistribution()
{
	UG_ASSERT(!parallel_redistribution_mode(),
			  "The dof distribution must not be in parallel distribution mode when"
			  " begin_parallel_redistribution is called.");

	m_parallelRedistributionMode = true;
	m_strictSubsetChecks = false;
}

void MGDoFDistribution::
end_parallel_redistribution()
{
	UG_ASSERT(parallel_redistribution_mode(),
			  "The dof distribution has to be in parallel distribution mode when"
			  " end_parallel_redistribution is called.");

	m_parallelRedistributionMode = false;
	m_strictSubsetChecks = true;
	parallel_redistribution_ended();
}

void MGDoFDistribution::check_subsets()
{
//	check, that all geom objects are assigned to a subset
	if(	m_spMGSH->num<VertexBase>() != multi_grid()->num<VertexBase>())
		UG_THROW("All Vertices "
			   " must be assigned to a subset. The passed subset handler "
			   " contains non-assigned elements, thus the dof distribution"
			   " is not possible, aborting.");

	if(	m_spMGSH->num<EdgeBase>() != multi_grid()->num<EdgeBase>())
		UG_THROW("All Edges "
			   " must be assigned to a subset. The passed subset handler "
			   " contains non-assigned elements, thus the dof distribution"
			   " is not possible, aborting.");

	if(	m_spMGSH->num<Face>() != multi_grid()->num<Face>())
		UG_THROW("All Faces "
			   " must be assigned to a subset. The passed subset handler "
			   " contains non-assigned elements, thus the dof distribution"
			   " is not possible, aborting.");

	if(	m_spMGSH->num<Volume>() != multi_grid()->num<Volume>())
		UG_THROW("All Volumes "
			   " must be assigned to a subset. The passed subset handler "
			   " contains non-assigned elements, thus the dof distribution"
			   " is not possible, aborting.");
}

////////////////////////////////////////////////////////////////////////////////
// MGDoFDistribution: Index Access

template <typename TBaseElem>
size_t MGDoFDistribution::
extract_inner_algebra_indices(TBaseElem* elem,
                              std::vector<size_t>& ind) const
{
//	get roid type and subset index
	const int si = m_spMGSH->get_subset_index(elem);
	static const ReferenceObjectID roid = elem->reference_object_id();

//	check if dofs present
	if(num_dofs(roid,si) > 0)
	{
	//	get index
		const size_t firstIndex = obj_index(elem);

		if(!m_bGrouped)
		{
			for(size_t fct = 0; fct < num_fct(); ++fct)
			{
			//	check that function is def on subset
				if(!is_def_in_subset(fct, si)) continue;

			//	get number of DoFs in this sub-geometric object
				const size_t numDoFsOnSub = num_dofs(fct,roid,roid);

			//	compute index
				const size_t index = firstIndex + offset(roid,si,fct);

			//	add dof to local indices
				for(size_t j = 0; j < numDoFsOnSub; ++j)
					ind.push_back(index+j);
			}
		}
		else
		{
		//	add dof to local indices
			ind.push_back(firstIndex);
		}
	}

//	return number of indices
	return ind.size();
}

template<typename TBaseElem>
void MGDoFDistribution::
extract_inner_algebra_indices(const typename Grid::traits<TBaseElem>::secure_container& vElem,
                              std::vector<size_t>& ind) const
{
//	loop passed elements
	for(size_t i = 0; i < vElem.size(); ++i)
		inner_algebra_indices(vElem[i], ind, false);
}

template<typename TBaseElem>
size_t MGDoFDistribution::inner_algebra_indices(TBaseElem* elem,
                                                std::vector<size_t>& ind,
                                                bool bClear) const
{
//	clear indices
	if(bClear) ind.clear();

//	return
	return extract_inner_algebra_indices(elem, ind);
}

template<typename TBaseElem>
size_t MGDoFDistribution::algebra_indices(TBaseElem* elem,
                                          std::vector<size_t>& ind,
                                          bool bClear) const
{
//	clear indices
	if(bClear) ind.clear();

//	reference dimension
	static const int dim = TBaseElem::dim;

//	get all sub-elements and add indices on those
	if(max_dofs(VERTEX) > 0)
	{
		//std::vector<VertexBase*> vElem;
		Grid::SecureVertexContainer vElem;
		m_pMG->associated_elements(vElem, elem);
		extract_inner_algebra_indices<VertexBase>(vElem, ind);
	}
	if(dim >= EDGE && max_dofs(EDGE) > 0)
	{
		Grid::SecureEdgeContainer vElem;
		m_pMG->associated_elements(vElem, elem);
		extract_inner_algebra_indices<EdgeBase>(vElem, ind);
	}
	if(dim >= FACE && max_dofs(FACE) > 0)
	{
		Grid::SecureFaceContainer vElem;
		m_pMG->associated_elements(vElem, elem);
		extract_inner_algebra_indices<Face>(vElem, ind);
	}
	if(dim >= VOLUME && max_dofs(VOLUME) > 0)
	{
		Grid::SecureVolumeContainer vElem;
		m_pMG->associated_elements(vElem, elem);
		extract_inner_algebra_indices<Volume>(vElem, ind);
	}

//	return number of indices
	return ind.size();
}

template<typename TBaseElem, typename TSubBaseElem>
void MGDoFDistribution::
multi_indices(TBaseElem* elem, const ReferenceObjectID roid,
              size_t fct, std::vector<multi_index_type>& ind,
              const typename Grid::traits<TSubBaseElem>::secure_container& vElem,
              const Grid::SecureVertexContainer& vCorner, bool bHang) const
{
//	get dimension of subelement
	static const int d = TSubBaseElem::dim;

//	loop passed elements
	for(size_t i = 0; i < vElem.size(); ++i)
	{
	//	get subelement
		TSubBaseElem* subElem = vElem[i];

	//	get subset index
		const int si = m_spMGSH->get_subset_index(subElem);

	//	check if function is defined on the subset
		if(!is_def_in_subset(fct, si)) continue;

	//	get reference object id for subselement
		const ReferenceObjectID subRoid = subElem->reference_object_id();

	//	check if dof given
		if(num_dofs(subRoid,si) == 0) continue;

	//	get number of DoFs in this sub-geometric object
		const size_t numDoFsOnSub = num_dofs(fct, roid, subRoid);

	//	a) Orientation required
		if(d <= max_dim_to_order_dofs(fct) && numDoFsOnSub > 1)
		{
			std::vector<size_t> vOrientOffset(numDoFsOnSub);

		//	get the orientation for this subelement
			ComputeOrientationOffset<TSubBaseElem>
				(vOrientOffset, local_finite_element_id(fct).order(),
				 ReferenceElementProvider::get(roid),
				 	 i, numDoFsOnSub, vCorner);

			if(!m_bGrouped)
			{
			//	compute index
				const size_t index = obj_index(subElem) + offset(subRoid,si,fct);

				for(size_t j = 0; j < numDoFsOnSub; ++j)
					ind.push_back(multi_index_type(index+vOrientOffset[j],0));
			}
			else
			{
			//	compute index
				const size_t comp = offset(subRoid,si,fct);
				const size_t firstIndex = obj_index(subElem);

				for(size_t j = 0; j < numDoFsOnSub; ++j)
					ind.push_back(multi_index_type(firstIndex, comp+vOrientOffset[j]));
			}
		}
	//	b) No Orientation required
		else
		{
			if(!m_bGrouped)
			{
			//	compute index
				const size_t index = obj_index(subElem) + offset(subRoid,si,fct);

				for(size_t j = 0; j < numDoFsOnSub; ++j)
					ind.push_back(multi_index_type(index+j,0));
			}
			else
			{
			//	compute index
				const size_t comp = offset(subRoid,si,fct);
				const size_t firstIndex = obj_index(subElem);

				for(size_t j = 0; j < numDoFsOnSub; ++j)
					ind.push_back(multi_index_type(firstIndex, comp+j));
			}
		} // end switch "Orientation required"
	} // end loop sub elements

//	return number of indices
	return;
}

template<typename TBaseElem>
size_t MGDoFDistribution::inner_multi_indices(TBaseElem* elem, size_t fct,
                                              std::vector<multi_index_type>& ind,
                                              bool bClear) const
{
//	clear if requested
	if(bClear) ind.clear();

//	get dimension
	static const int d = TBaseElem::dim;

//	get subset index
	const int si = m_spMGSH->get_subset_index(elem);

//	check if function is defined on the subset
	if(!is_def_in_subset(fct, si)) return ind.size();

//	get roid type
	static const ReferenceObjectID roid = elem->reference_object_id();

//	get number of DoFs in this sub-geometric object
	const size_t numDoFsOnSub = num_dofs(fct,roid,roid);

//	check if dof given
	if(numDoFsOnSub == 0) return ind.size();

//	a) Orientation required:
	if(d <= max_dim_to_order_dofs(fct) && numDoFsOnSub > 1)
	{
	//	get corners
		Grid::SecureVertexContainer vCorner;
		m_pMG->associated_elements(vCorner, elem);
	//	get the orientation for this
		std::vector<size_t> vOrientOffset(numDoFsOnSub);
		ComputeOrientationOffset<TBaseElem>
			(vOrientOffset, local_finite_element_id(fct).order(),
			 ReferenceElementProvider::get(roid),
			 0, numDoFsOnSub, vCorner);

		if(!m_bGrouped)
		{
		//	compute index
			const size_t index = obj_index(elem) + offset(roid,si,fct);

			for(size_t j = 0; j < numDoFsOnSub; ++j)
				ind.push_back(multi_index_type(index+vOrientOffset[j],0));
		}
		else
		{
		//	compute index
			const size_t comp = offset(roid,si,fct);
			const size_t firstIndex = obj_index(elem);

			for(size_t j = 0; j < numDoFsOnSub; ++j)
				ind.push_back(multi_index_type(firstIndex, comp+vOrientOffset[j]));
		}
	}
//	b) No orientation needed
	else
	{
		if(!m_bGrouped)
		{
		//	compute index
			const size_t index = obj_index(elem) + offset(roid,si,fct);

			for(size_t j = 0; j < numDoFsOnSub; ++j)
				ind.push_back(multi_index_type(index+j,0));
		}
		else
		{
		//	compute index
			const size_t comp = offset(roid,si,fct);
			const size_t firstIndex = obj_index(elem);

			for(size_t j = 0; j < numDoFsOnSub; ++j)
				ind.push_back(multi_index_type(firstIndex, comp+j));
		}
	}

//	done
	return ind.size();
}

template<typename TBaseElem>
size_t MGDoFDistribution::multi_indices(TBaseElem* elem, size_t fct,
                                        std::vector<multi_index_type>& ind,
                                        bool bHang, bool bClear) const
{
//	clear indices
	if(bClear) ind.clear();

//	reference dimension
	static const int dim = TBaseElem::dim;

//	reference object id
	const ReferenceObjectID roid = elem->reference_object_id();

//	get all sub-elements and add indices on those
	Grid::SecureVertexContainer vCorner;
	m_pMG->associated_elements_sorted(vCorner, elem);
	if(max_dofs(VERTEX) > 0)
	{
		multi_indices<TBaseElem, VertexBase>(elem, roid, fct, ind, vCorner, vCorner, bHang);
	}
	if(dim >= EDGE && max_dofs(EDGE) > 0)
	{
		Grid::SecureEdgeContainer vElem;
		m_pMG->associated_elements_sorted(vElem, elem);
		multi_indices<TBaseElem, EdgeBase>(elem, roid, fct, ind, vElem, vCorner, bHang);
	}
	if(dim >= FACE && max_dofs(FACE) > 0)
	{
		Grid::SecureFaceContainer vElem;
		m_pMG->associated_elements_sorted(vElem, elem);
		multi_indices<TBaseElem, Face>(elem, roid, fct, ind, vElem, vCorner, bHang);
	}
	if(dim >= VOLUME && max_dofs(VOLUME) > 0)
	{
		Grid::SecureVolumeContainer vElem;
		m_pMG->associated_elements_sorted(vElem, elem);
		multi_indices<TBaseElem, Volume>(elem, roid, fct, ind, vElem, vCorner, bHang);
	}

//	todo: add hanging nodes
//	If no hanging dofs are required, we're done
	if(!bHang) return ind.size();

	UG_THROW("Hanging DoFs are currently not supported by this DoFManager.");

//	return number of indices
	return ind.size();
}

template<typename TBaseElem>
void MGDoFDistribution::indices_on_vertex(TBaseElem* elem, const ReferenceObjectID roid,
                                          LocalIndices& ind,
                                          const Grid::SecureVertexContainer& vElem) const
{
//	get reference object id for subelement
	static const ReferenceObjectID subRoid = ROID_VERTEX;

//	add normal dofs
	for(size_t i = 0; i < vElem.size(); ++i)
	{
	//	get subset index
		const int si = m_spMGSH->get_subset_index(vElem[i]);

	//	loop all functions
		for(size_t fct = 0; fct < num_fct(); ++fct)
		{
		//	check if function is defined on the subset
			if(!is_def_in_subset(fct, si)) continue;

		//	get number of DoFs in this sub-geometric object
			const size_t numDoFsOnSub = num_dofs(fct,roid,subRoid);

		//	Always no orientation needed
			if(!m_bGrouped)
			{
			//	compute index
				const size_t index = obj_index(vElem[i]) + offset(subRoid,si,fct);

			//	add dof to local indices
				for(size_t j = 0; j < numDoFsOnSub; ++j)
					ind.push_back_index(fct, index+j);
			}
			else
			{
			//	compute index
				const size_t index = obj_index(vElem[i]);
				const size_t comp = offset(subRoid,si,fct);

			//	add dof to local indices
				for(size_t j = 0; j < numDoFsOnSub; ++j)
					ind.push_back_multi_index(fct, index, comp+j);
			}
		} // end loop functions
	} // end loop subelement

}

template<typename TBaseElem, typename TSubBaseElem>
void MGDoFDistribution::indices(TBaseElem* elem, const ReferenceObjectID roid,
                                LocalIndices& ind,
                                const typename Grid::traits<TSubBaseElem>::secure_container& vElem,
                                const Grid::SecureVertexContainer& vCorner) const
{
//	dimension of subelement
	static const int d = TSubBaseElem::dim;

//	add normal dofs
	for(size_t i = 0; i < vElem.size(); ++i)
	{
	//	get subelement
		TSubBaseElem* subElem = vElem[i];

	//	get subset index
		const int si = m_spMGSH->get_subset_index(subElem);

	//	get reference object id for subselement
		const ReferenceObjectID subRoid = subElem->reference_object_id();

	//	loop all functions
		for(size_t fct = 0; fct < num_fct(); ++fct)
		{
		//	check if function is defined on the subset
			if(!is_def_in_subset(fct, si)) continue;

		//	get number of DoFs in this sub-geometric object
			const size_t numDoFsOnSub = num_dofs(fct,roid,subRoid);

		//	a)	Orientation is required: Thus, we compute the offsets, that are
		//		no longer in the usual order [0, 1, 2, ...]. Orientation is
		//		required if there are more than 1 dof on a subelement of a
		//		finite element and thus, when gluing two elements together,
		//		also the dofs on the subelements have to fit in order to
		//		guarantee continuity. This is not needed for Vertices, since there
		//		no distinction can be made when all dofs are at the same position.
		//		This is also not needed for the highest dimension of a finite
		//		element, since the dofs on this geometric object must not be
		//		identified with other dofs.
			if(d <= max_dim_to_order_dofs(fct) && numDoFsOnSub > 1)
			{
			//	vector storing the computed offsets. If in correct order,
			//	this would be: [0, 1, 2, ...]. But this is usually not the
			// 	case and the numbers 0 to numDoFsOnSub-1 are permuted
				std::vector<size_t> vOrientOffset(numDoFsOnSub);

				ComputeOrientationOffset<TSubBaseElem>
					(vOrientOffset, local_finite_element_id(fct).order(),
					 ReferenceElementProvider::get(roid),
					 i, numDoFsOnSub, vCorner);

				if(!m_bGrouped)
				{
					const size_t index = obj_index(subElem) + offset(subRoid,si,fct);
					for(size_t j = 0; j < numDoFsOnSub; ++j)
						ind.push_back_index(fct, index+vOrientOffset[j]);
				}
				else
				{
				//	compute index
					const size_t index = obj_index(subElem);
					const size_t comp = offset(subRoid,si,fct);

					for(size_t j = 0; j < numDoFsOnSub; ++j)
						ind.push_back_multi_index(fct, index, comp+vOrientOffset[j]);
				}
			}
		//	b)	No orientation needed
			else
			{
				if(!m_bGrouped)
				{
				//	compute index
					const size_t index = obj_index(subElem) + offset(subRoid,si,fct);

				//	add dof to local indices
					for(size_t j = 0; j < numDoFsOnSub; ++j)
						ind.push_back_index(fct, index+j);
				}
				else
				{
				//	compute index
					const size_t index = obj_index(subElem);
					const size_t comp = offset(subRoid,si,fct);

				//	add dof to local indices
					for(size_t j = 0; j < numDoFsOnSub; ++j)
						ind.push_back_multi_index(fct, index, comp+j);
				}
			} // end switch "Orientation needed"
		} // end loop functions
	} // end loop subelement

}

template <typename TConstraining, typename TConstrained, typename TBaseElem>
void MGDoFDistribution::
constrained_indices(LocalIndices& ind,
                    const typename Grid::traits<TBaseElem>::secure_container& vSubElem) const
{
//	loop all edges
	for(size_t i = 0; i < vSubElem.size(); ++i)
	{
	//	only constraining objects are of interest
		TConstraining* constrainingObj = dynamic_cast<TConstraining*>(vSubElem[i]);
		if(constrainingObj == NULL) continue;

	//	loop constraining vertices
		for(size_t i = 0; i != constrainingObj->num_constrained_vertices(); ++i)
		{
		//	get vertex
			TConstrained* vrt = constrainingObj->constrained_vertex(i);

		//	get roid
			const ReferenceObjectID subRoid = vrt->reference_object_id();

		//	get subset index
			int si = m_spMGSH->get_subset_index(vrt);

		//	loop functions
			for(size_t fct = 0; fct < num_fct(); ++fct)
			{
			//	check that function is defined on subset
				if(!is_def_in_subset(fct, si)) continue;

				if(!m_bGrouped)
				{
				//	compute index
					const size_t index = obj_index(vrt) + offset(subRoid,si,fct);

				//	add dof to local indices
					ind.push_back_index(fct, index);
				}
				else
				{
				//	compute index
					const size_t index = obj_index(vrt);
					const size_t comp = offset(subRoid,si,fct);

				//	add dof to local indices
					ind.push_back_multi_index(fct, index, comp);
				}
			}
		}
	}
}

void MGDoFDistribution::local_finite_element_ids(LocalIndices& ind) const
{
	ind.resize_fct(num_fct());
	for(size_t fct = 0; fct < num_fct(); ++fct)
		ind.set_lfeID(fct, local_finite_element_id(fct));
}

template<typename TBaseElem>
void MGDoFDistribution::indices(TBaseElem* elem, LocalIndices& ind, bool bHang) const
{
//	reference dimension
	static const int dim = TBaseElem::dim;

//	resize the number of functions
	ind.resize_fct(num_fct());
	for(size_t fct = 0; fct < num_fct(); ++fct) ind.clear_dof(fct);

//	get all sub-elements and add indices on those
	Grid::SecureVertexContainer vCorner;
	m_pMG->associated_elements_sorted(vCorner, elem);

//	storage for (maybe needed) subelements
	Grid::SecureEdgeContainer vEdge;
	Grid::SecureFaceContainer vFace;
	Grid::SecureVolumeContainer vVol;

//	collect elements, if needed
	if(dim >= EDGE)
		if(max_dofs(EDGE) > 0 || bHang) m_pMG->associated_elements_sorted(vEdge, elem);
	if(dim >= FACE)
		if(max_dofs(FACE) > 0 || bHang) m_pMG->associated_elements_sorted(vFace, elem);
	if(dim >= VOLUME)
		if(max_dofs(VOLUME) > 0 || bHang) m_pMG->associated_elements_sorted(vVol, elem);

//	get reference object id
	const ReferenceObjectID roid = elem->reference_object_id();

//	get regular dofs on all subelements and the element itself
//	use specialized function for vertices (since only one position and one reference object)
	if(max_dofs(VERTEX) > 0) 				  indices_on_vertex<TBaseElem>(elem, roid, ind, vCorner);
	if(dim >= EDGE && max_dofs(EDGE) > 0) 	  indices<TBaseElem, EdgeBase>(elem, roid, ind, vEdge, vCorner);
	if(dim >= FACE && max_dofs(FACE) > 0) 	  indices<TBaseElem, Face>(elem, roid, ind, vFace, vCorner);
	if(dim >= VOLUME && max_dofs(VOLUME) > 0) indices<TBaseElem, Volume>(elem, roid, ind, vVol, vCorner);

//	If no hanging dofs are required, we're done
	if(!bHang) return;

//	get dofs on hanging vertices
	if(max_dofs(VERTEX > 0))
	{
		if(dim >= EDGE) constrained_indices<ConstrainingEdge, VertexBase, EdgeBase>(ind, vEdge);
		if(dim >= FACE) constrained_indices<ConstrainingQuadrilateral, VertexBase, Face>(ind, vFace);
	}

//	todo: allow also constrained dofs on other elements
	if(max_dofs(EDGE) || max_dofs(FACE) || max_dofs(VOLUME))
		UG_THROW("Hanging DoFs are only implemented for P1 by this DoFManager.");

//	we're done
	return;
}



template <typename TBaseElem>
void MGDoFDistribution::
changable_indices(std::vector<size_t>& vIndex,
                  const std::vector<TBaseElem*>& vElem) const
{
//	Get connected indices
	for(size_t i = 0; i < vElem.size(); ++i)
	{
	//	Get Vertices of adjacent edges
		TBaseElem* elem = vElem[i];
		UG_ASSERT(m_spMGSH->get_subset_index(elem) >= 0, "Must have subset");

	//	get adjacent index
		const size_t adjInd = obj_index(elem);

	//	add to index list
		vIndex.push_back(adjInd);
	}
}


////////////////////////////////////////////////////////////////////////////////
// MGDoFDistribution: DoF Handling


void MGDoFDistribution::init_attachments()
{
//	attach DoFs to vertices
	if(max_dofs(VERTEX)) {
		multi_grid()->attach_to<VertexBase>(m_aIndex);
		m_aaIndexVRT.access(*m_pMG, m_aIndex);
	}
	if(max_dofs(EDGE)) {
		multi_grid()->attach_to<EdgeBase>(m_aIndex);
		m_aaIndexEDGE.access(*m_pMG, m_aIndex);
	}
	if(max_dofs(FACE)) {
		multi_grid()->attach_to<Face>(m_aIndex);
		m_aaIndexFACE.access(*m_pMG, m_aIndex);
	}
	if(max_dofs(VOLUME)) {
		multi_grid()->attach_to<Volume>(m_aIndex);
		m_aaIndexVOL.access(*m_pMG, m_aIndex);
	}
}

void MGDoFDistribution::clear_attachments()
{
//	detach DoFs
	if(m_aaIndexVRT.valid()) multi_grid()->detach_from<VertexBase>(m_aIndex);
	if(m_aaIndexEDGE.valid()) multi_grid()->detach_from<EdgeBase>(m_aIndex);
	if(m_aaIndexFACE.valid()) multi_grid()->detach_from<Face>(m_aIndex);
	if(m_aaIndexVOL.valid()) multi_grid()->detach_from<Volume>(m_aIndex);

	m_aaIndexVRT.invalidate();
	m_aaIndexEDGE.invalidate();
	m_aaIndexFACE.invalidate();
	m_aaIndexVOL.invalidate();
}

void MGDoFDistribution::indices(GeometricObject* elem, LocalIndices& ind, bool bHang) const
{
	switch(elem->base_object_id())
	{
		case VERTEX: return indices(static_cast<VertexBase*>(elem), ind, bHang);
		case EDGE: return indices(static_cast<EdgeBase*>(elem), ind, bHang);
		case FACE: return indices(static_cast<Face*>(elem), ind, bHang);
		case VOLUME: return indices(static_cast<Volume*>(elem), ind, bHang);
		default: UG_THROW("Geometric Base element not found.");
	}
}

size_t MGDoFDistribution::multi_indices(GeometricObject* elem, size_t fct,
                                        std::vector<multi_index_type>& ind,
                                        bool bHang, bool bClear) const
{
	switch(elem->base_object_id())
	{
		case VERTEX: return multi_indices(static_cast<VertexBase*>(elem), fct, ind, bHang, bClear);
		case EDGE: return multi_indices(static_cast<EdgeBase*>(elem), fct, ind, bHang, bClear);
		case FACE: return multi_indices(static_cast<Face*>(elem), fct, ind, bHang, bClear);
		case VOLUME: return multi_indices(static_cast<Volume*>(elem), fct, ind, bHang, bClear);
		default: UG_THROW("Geometric Base element not found.");
	}
}

size_t MGDoFDistribution::inner_multi_indices(GeometricObject* elem, size_t fct,
                                              std::vector<multi_index_type>& ind,
                                              bool bClear) const
{
	switch(elem->base_object_id())
	{
		case VERTEX: return inner_multi_indices(static_cast<VertexBase*>(elem), fct, ind, bClear);
		case EDGE: return inner_multi_indices(static_cast<EdgeBase*>(elem), fct, ind, bClear);
		case FACE: return inner_multi_indices(static_cast<Face*>(elem), fct, ind, bClear);
		case VOLUME: return inner_multi_indices(static_cast<Volume*>(elem), fct, ind, bClear);
		default: UG_THROW("Geometric Base element not found.");
	}
}

size_t MGDoFDistribution::algebra_indices(GeometricObject* elem,	std::vector<size_t>& ind,
                                          bool bClear) const
{
	switch(elem->base_object_id())
	{
		case VERTEX: return algebra_indices(static_cast<VertexBase*>(elem), ind, bClear);
		case EDGE: return algebra_indices(static_cast<EdgeBase*>(elem), ind, bClear);
		case FACE: return algebra_indices(static_cast<Face*>(elem), ind, bClear);
		case VOLUME: return algebra_indices(static_cast<Volume*>(elem), ind, bClear);
		default: UG_THROW("Geometric Base element not found.");
	}
}

size_t MGDoFDistribution::inner_algebra_indices(GeometricObject* elem, std::vector<size_t>& ind,
                                                bool bClear) const
{
	switch(elem->base_object_id())
	{
		case VERTEX: return inner_algebra_indices(static_cast<VertexBase*>(elem), ind, bClear);
		case EDGE: return inner_algebra_indices(static_cast<EdgeBase*>(elem), ind, bClear);
		case FACE: return inner_algebra_indices(static_cast<Face*>(elem), ind, bClear);
		case VOLUME: return inner_algebra_indices(static_cast<Volume*>(elem), ind, bClear);
		default: UG_THROW("Geometric Base element not found.");
	}
}

void MGDoFDistribution::register_observer()
{
	int type = OT_GRID_OBSERVER;

	if(max_dofs(VERTEX)) type |= OT_VERTEX_OBSERVER;
	if(max_dofs(EDGE)) type |= OT_EDGE_OBSERVER;
	if(max_dofs(FACE)) type |= OT_FACE_OBSERVER;
	if(max_dofs(VOLUME)) type |= OT_VOLUME_OBSERVER;

#ifdef UG_PARALLEL
	if(pcl::GetNumProcesses() > 1){
	//	to correctly support parallel coarsening we have to register the dof-distribution
	//	as a full observer.
		type = OT_FULL_OBSERVER;
	}
#endif

	multi_grid()->register_observer(this, type);
}

void MGDoFDistribution::unregister_observer()
{
	multi_grid()->unregister_observer(this);
}

///////////////////////////////////////////////////////////////////////////////
// template instantiations
///////////////////////////////////////////////////////////////////////////////
template void MGDoFDistribution::indices<VertexBase>(VertexBase*, LocalIndices&, bool) const;
template void MGDoFDistribution::indices<EdgeBase>(EdgeBase*, LocalIndices&, bool) const;
template void MGDoFDistribution::indices<Face>(Face*, LocalIndices&, bool) const;
template void MGDoFDistribution::indices<Volume>(Volume*, LocalIndices&, bool) const;

template size_t MGDoFDistribution::multi_indices<VertexBase>(VertexBase*, size_t, std::vector<multi_index_type>&, bool, bool) const;
template size_t MGDoFDistribution::multi_indices<EdgeBase>(EdgeBase*, size_t, std::vector<multi_index_type>&, bool, bool) const;
template size_t MGDoFDistribution::multi_indices<Face>(Face*, size_t, std::vector<multi_index_type>&, bool, bool) const;
template size_t MGDoFDistribution::multi_indices<Volume>(Volume*, size_t, std::vector<multi_index_type>&, bool, bool) const;

template size_t MGDoFDistribution::inner_multi_indices<VertexBase>(VertexBase*, size_t,std::vector<multi_index_type>&, bool) const;
template size_t MGDoFDistribution::inner_multi_indices<EdgeBase>(EdgeBase*, size_t,std::vector<multi_index_type>&, bool) const;
template size_t MGDoFDistribution::inner_multi_indices<Face>(Face*, size_t,std::vector<multi_index_type>&, bool) const;
template size_t MGDoFDistribution::inner_multi_indices<Volume>(Volume*, size_t,std::vector<multi_index_type>&, bool) const;

template size_t MGDoFDistribution::algebra_indices<VertexBase>(VertexBase*,	std::vector<size_t>&, bool) const;
template size_t MGDoFDistribution::algebra_indices<EdgeBase>(EdgeBase*,	std::vector<size_t>&, bool) const;
template size_t MGDoFDistribution::algebra_indices<Face>(Face*,	std::vector<size_t>&, bool) const;
template size_t MGDoFDistribution::algebra_indices<Volume>(Volume*,	std::vector<size_t>&, bool) const;

template size_t MGDoFDistribution::inner_algebra_indices<VertexBase>(VertexBase*, std::vector<size_t>& , bool) const;
template size_t MGDoFDistribution::inner_algebra_indices<EdgeBase>(EdgeBase*, std::vector<size_t>& , bool) const;
template size_t MGDoFDistribution::inner_algebra_indices<Face>(Face*, std::vector<size_t>& , bool) const;
template size_t MGDoFDistribution::inner_algebra_indices<Volume>(Volume*, std::vector<size_t>& , bool) const;

template void MGDoFDistribution::changable_indices<VertexBase>(std::vector<size_t>& vIndex, const std::vector<VertexBase*>& vElem) const;
template void MGDoFDistribution::changable_indices<EdgeBase>(std::vector<size_t>& vIndex, const std::vector<EdgeBase*>& vElem) const;
template void MGDoFDistribution::changable_indices<Face>(std::vector<size_t>& vIndex, const std::vector<Face*>& vElem) const;
template void MGDoFDistribution::changable_indices<Volume>(std::vector<size_t>& vIndex, const std::vector<Volume*>& vElem) const;

} // end namespace ug
