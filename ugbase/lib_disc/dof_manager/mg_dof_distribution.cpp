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
#include "lib_disc/reference_element/reference_element.h"
#include "lib_disc/reference_element/reference_element_traits.h"
#include "lib_disc/reference_element/reference_mapping.h"
#include "lib_disc/reference_element/reference_mapping_provider.h"
#include "lib_disc/common/groups_util.h"
#include "common/util/string_util.h"
#include "orientation.h"

using namespace std;

namespace ug{


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
	  m_spMG(spMG),
	  m_pMG(spMG.get()),
	  m_spMGSH(spMGSH)
{
	check_subsets();
	init_attachments();
};

MGDoFDistribution::
~MGDoFDistribution()
{
	clear_attachments();
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
	const ReferenceObjectID roid = elem->reference_object_id();

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
				const size_t numDoFsOnSub = num_fct_dofs(fct,roid,si);

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
		Grid::SecureVertexContainer vVrt;
		m_pMG->associated_elements(vVrt, elem);
		extract_inner_algebra_indices<VertexBase>(vVrt, ind);
	}
	if(dim >= EDGE && max_dofs(EDGE) > 0)
	{
		Grid::SecureEdgeContainer vEdge;
		m_pMG->associated_elements(vEdge, elem);
		extract_inner_algebra_indices<EdgeBase>(vEdge, ind);
	}
	if(dim >= FACE && max_dofs(FACE) > 0)
	{
		Grid::SecureFaceContainer vFace;
		m_pMG->associated_elements(vFace, elem);
		extract_inner_algebra_indices<Face>(vFace, ind);
	}
	if(dim >= VOLUME && max_dofs(VOLUME) > 0)
	{
		Grid::SecureVolumeContainer vVol;
		m_pMG->associated_elements(vVol, elem);
		extract_inner_algebra_indices<Volume>(vVol, ind);
	}

//	return number of indices
	return ind.size();
}

template<typename TBaseElem, typename TSubBaseElem>
void MGDoFDistribution::
dof_indices(TBaseElem* elem, const ReferenceObjectID roid,
              size_t fct, std::vector<DoFIndex>& ind,
              const typename Grid::traits<TSubBaseElem>::secure_container& vElem) const
{
//	storage for offsets
	std::vector<size_t> vOrientOffset;

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
		const size_t numDoFsOnSub = num_fct_dofs(fct, subRoid, si);

	//	get the orientation for this subelement
		ComputeOrientationOffset(vOrientOffset, elem, subElem, i, lfeid(fct));

		UG_ASSERT(vOrientOffset.size() == numDoFsOnSub ||
		          vOrientOffset.empty(), "Something wrong with orientation");

		if(!m_bGrouped)
		{
		//	compute index
			const size_t index = obj_index(subElem) + offset(subRoid,si,fct);

			if(vOrientOffset.empty()){
				for(size_t j = 0; j < numDoFsOnSub; ++j)
					ind.push_back(DoFIndex(index + j,0));
			}
			else{
				for(size_t j = 0; j < numDoFsOnSub; ++j)
					ind.push_back(DoFIndex(index + vOrientOffset[j],0));
			}
		}
		else
		{
		//	compute index
			const size_t comp = offset(subRoid,si,fct);
			const size_t firstIndex = obj_index(subElem);

			if(vOrientOffset.empty()){
				for(size_t j = 0; j < numDoFsOnSub; ++j)
					ind.push_back(DoFIndex(firstIndex, comp + j));
			}
			else{
				for(size_t j = 0; j < numDoFsOnSub; ++j)
					ind.push_back(DoFIndex(firstIndex, comp + vOrientOffset[j]));
			}
		}

	} // end loop sub elements

//	return number of indices
	return;
}

template<typename TBaseElem>
size_t MGDoFDistribution::inner_dof_indices(TBaseElem* elem, size_t fct,
                                              std::vector<DoFIndex>& ind,
                                              bool bClear) const
{
//	clear if requested
	if(bClear) ind.clear();

//	get subset index
	const int si = m_spMGSH->get_subset_index(elem);

//	check if function is defined on the subset
	if(!is_def_in_subset(fct, si)) return ind.size();

//	get roid type
	const ReferenceObjectID roid = elem->reference_object_id();

//	get number of DoFs in this sub-geometric object
	const size_t numDoFsOnSub = num_fct_dofs(fct,roid,si);

//	check if dof given
	if(numDoFsOnSub == 0) return ind.size();

//	Note: No orientation needed
	if(!m_bGrouped)
	{
	//	compute index
		const size_t index = obj_index(elem) + offset(roid,si,fct);

		for(size_t j = 0; j < numDoFsOnSub; ++j)
			ind.push_back(DoFIndex(index+j,0));
	}
	else
	{
	//	compute index
		const size_t comp = offset(roid,si,fct);
		const size_t firstIndex = obj_index(elem);

		for(size_t j = 0; j < numDoFsOnSub; ++j)
			ind.push_back(DoFIndex(firstIndex, comp+j));
	}

//	done
	return ind.size();
}

template <typename TConstraining, typename TConstrained, typename TBaseElem>
void MGDoFDistribution::
constrained_vertex_dof_indices(size_t fct,std::vector<DoFIndex>& ind,
                    const typename Grid::traits<TBaseElem>::secure_container& vSubElem) const
{
	//	loop all edges
	for(size_t i = 0; i < vSubElem.size(); ++i)
	{
	//	only constraining objects are of interest
		TConstraining* constrainingObj = dynamic_cast<TConstraining*>(vSubElem[i]);
		if(constrainingObj == NULL) continue;
		
		//	get subset index
		const int si = m_spMGSH->get_subset_index(vSubElem[i]);

	//	loop constraining vertices
		for(size_t j = 0; j != constrainingObj->num_constrained_vertices(); ++j)
		{
		//	get vertex
			TConstrained* vrt = constrainingObj->constrained_vertex(j);

		//	get roid
			const ReferenceObjectID subRoid = vrt->reference_object_id();
			
		//	check if dof given
			if(num_dofs(subRoid,si) == 0) continue;

		//	get subset index
			int si = m_spMGSH->get_subset_index(vrt);
		
		//	check that function is defined on subset
			if(!is_def_in_subset(fct, si)) continue;
			
		//	get number of DoFs in this sub-geometric object
			const size_t numDoFsOnSub = num_fct_dofs(fct, subRoid, si);

			if(!m_bGrouped)
			{
			//	compute index
				const size_t index = obj_index(vrt) + offset(subRoid,si,fct);

				for (size_t k=0;k<numDoFsOnSub;k++)
					ind.push_back(DoFIndex(index + k,0));
			}
			else
			{
			//	compute index
				const size_t index = obj_index(vrt);
				const size_t comp = offset(subRoid,si,fct);

			//	add dof to local indices
				for(size_t k = 0; k < numDoFsOnSub; ++k)
					ind.push_back(DoFIndex(index, comp + k));
			}
		}
	}
}

template <typename TBaseElem,typename TConstraining, typename TConstrained, typename TSubElem>
void MGDoFDistribution::
constrained_edge_dof_indices(TBaseElem* elem,size_t fct,std::vector<DoFIndex>& ind,
                    const typename Grid::traits<TSubElem>::secure_container& vSubElem) const
{
	//	storage for offsets
	std::vector<size_t> vOrientOffset;

	//	loop all edges
	for(size_t i = 0; i < vSubElem.size(); ++i)
	{
	//	only constraining objects are of interest
		TConstraining* constrainingObj = dynamic_cast<TConstraining*>(vSubElem[i]);
		if(constrainingObj == NULL) continue;

		std::vector<size_t> sortedInd;
		sort_constrained_edges<TBaseElem,TConstraining,TConstrained>(sortedInd,elem,constrainingObj,i);
		
	//	get the orientation for this subelement
		ComputeOrientationOffset(vOrientOffset, elem, constrainingObj, i, lfeid(fct));

	//  loop constraining edges
		for(size_t j = 0; j != constrainingObj->num_constrained_edges(); ++j)
		{
			//	get edge
			TConstrained* edg = constrainingObj->constrained_edge(sortedInd[j]);

			//	get roid
			const ReferenceObjectID subRoid = edg->reference_object_id();

			//	get subset index
			int si = m_spMGSH->get_subset_index(edg);

			
			//	check that function is defined on subset
			if(!is_def_in_subset(fct, si)) continue;
			
			//	check if dof given
			if(num_dofs(subRoid,si) == 0) continue;
			
			const size_t numDoFsOnSub=num_fct_dofs(fct, subRoid, si);

			if(!m_bGrouped)
			{
			//	compute index
				const size_t index = obj_index(edg) + offset(subRoid,si,fct);

				if(vOrientOffset.empty()){
					for(size_t k = 0; k < numDoFsOnSub; ++k)
						ind.push_back(DoFIndex(index + k,0));
				}
				else{
					for(size_t k = 0; k < numDoFsOnSub; ++k)
						ind.push_back(DoFIndex(index + vOrientOffset[k],0));
				}
			}
			else
			{
				//	compute index
				const size_t index = obj_index(edg);
				const size_t comp = offset(subRoid,si,fct);
				
				if(vOrientOffset.empty()){
				for(size_t k = 0; k < numDoFsOnSub; ++k)
					ind.push_back(DoFIndex(index, comp + k));
				}
				else{
					for(size_t k = 0; k < numDoFsOnSub; ++k)
						ind.push_back(DoFIndex(index, comp + vOrientOffset[k]));
				}
			}
		}
	}
}

template <typename TBaseElem,typename TConstraining, typename TConstrained, typename TSubElem>
void MGDoFDistribution::
constrained_face_dof_indices(TBaseElem* elem,size_t fct,std::vector<DoFIndex>& ind,
                    const typename Grid::traits<TSubElem>::secure_container& vSubElem) const
{
	//	storage for offsets
	std::vector<size_t> vOrientOffset;

	//	loop all edges
	for(size_t i = 0; i < vSubElem.size(); ++i)
	{
	//	only constraining objects are of interest
		TConstraining* constrainingObj = dynamic_cast<TConstraining*>(vSubElem[i]);
		if(constrainingObj == NULL) continue;

		std::vector<size_t> sortedInd;
		sort_constrained_faces<TBaseElem,TConstraining,TConstrained>(sortedInd,elem,constrainingObj,i);
		
	//	get the orientation for this subelement
		ComputeOrientationOffset(vOrientOffset, elem, constrainingObj, i, lfeid(fct));

	//  loop constraining edges
		for(size_t j = 0; j != constrainingObj->num_constrained_faces(); ++j)
		{
			//	get face
			TConstrained* face = constrainingObj->constrained_face(sortedInd[j]);

			//	get roid
			const ReferenceObjectID subRoid = face->reference_object_id();

			//	get subset index
			int si = m_spMGSH->get_subset_index(face);

			
			//	check that function is defined on subset
			if(!is_def_in_subset(fct, si)) continue;
			
			//	check if dof given
			if(num_dofs(subRoid,si) == 0) continue;
			
			const size_t numDoFsOnSub=num_fct_dofs(fct, subRoid, si);

			if(!m_bGrouped)
			{
			//	compute index
				const size_t index = obj_index(face) + offset(subRoid,si,fct);

				if(vOrientOffset.empty()){
					for(size_t k = 0; k < numDoFsOnSub; ++k)
						ind.push_back(DoFIndex(index + k,0));
				}
				else{
					for(size_t k = 0; k < numDoFsOnSub; ++k)
						ind.push_back(DoFIndex(index + vOrientOffset[k],0));
				}
			}
			else
			{
				//	compute index
				const size_t index = obj_index(face);
				const size_t comp = offset(subRoid,si,fct);
				
				if(vOrientOffset.empty()){
				for(size_t k = 0; k < numDoFsOnSub; ++k)
					ind.push_back(DoFIndex(index, comp + k));
				}
				else{
					for(size_t k = 0; k < numDoFsOnSub; ++k)
						ind.push_back(DoFIndex(index, comp + vOrientOffset[k]));
				}
			}
		}
	}
}


template<typename TBaseElem>
size_t MGDoFDistribution::dof_indices(TBaseElem* elem, size_t fct,
                                        std::vector<DoFIndex>& ind,
                                        bool bHang, bool bClear) const
{
//	clear indices
	if(bClear) ind.clear();

//	reference dimension
	static const int dim = TBaseElem::dim;

//	reference object id
	const ReferenceObjectID roid = elem->reference_object_id();
	
//	storage for (maybe needed) subelements
	Grid::SecureVertexContainer vCorner;
	Grid::SecureEdgeContainer vEdge;
	Grid::SecureFaceContainer vFace;
	Grid::SecureVolumeContainer vVol;

//	collect elements, if needed
	if(dim >= VERTEX)
		if(max_dofs(VERTEX) > 0) m_pMG->associated_elements_sorted(vCorner, elem);
	if(dim >= EDGE)
		if(max_dofs(EDGE) > 0 || bHang) m_pMG->associated_elements_sorted(vEdge, elem);
	if(dim >= FACE)
		if(max_dofs(FACE) > 0 || bHang) m_pMG->associated_elements_sorted(vFace, elem);
	if(dim >= VOLUME)
		if(max_dofs(VOLUME) > 0) m_pMG->associated_elements_sorted(vVol, elem);

//	get regular dofs on all subelements and the element itself
//	use specialized function for vertices (since only one position and one reference object)
	if(dim >= VERTEX && max_dofs(VERTEX) > 0) dof_indices<TBaseElem, VertexBase>(elem, roid, fct, ind, vCorner);
	if(dim >= EDGE && max_dofs(EDGE) > 0) 	  dof_indices<TBaseElem, EdgeBase>(elem, roid, fct, ind, vEdge);
	if(dim >= FACE && max_dofs(FACE) > 0) 	  dof_indices<TBaseElem, Face>(elem, roid, fct, ind, vFace);
	if(dim >= VOLUME && max_dofs(VOLUME) > 0) dof_indices<TBaseElem, Volume>(elem, roid, fct, ind, vVol);

//	If no hanging dofs are required, we're done
	if(!bHang) return ind.size();

	//	get dofs on hanging vertices
	if(max_dofs(VERTEX > 0))
	{
		if(dim >= EDGE) constrained_vertex_dof_indices<ConstrainingEdge, VertexBase, EdgeBase>(fct,ind,vEdge);
		if(dim >= FACE) constrained_vertex_dof_indices<ConstrainingQuadrilateral, VertexBase, Face>(fct,ind,vFace);
	}

//	get dofs on hanging edges
	if (max_dofs(EDGE) > 0){
		if(dim >= EDGE) constrained_edge_dof_indices<TBaseElem,ConstrainingEdge, EdgeBase, EdgeBase>(elem,fct,ind, vEdge);
		if(dim >= FACE) constrained_edge_dof_indices<TBaseElem,ConstrainingTriangle, EdgeBase, Face>(elem,fct,ind, vFace);
		if(dim >= FACE) constrained_edge_dof_indices<TBaseElem,ConstrainingQuadrilateral, EdgeBase, Face>(elem,fct,ind, vFace);
	}

//  get dofs on hanging faces
	if (max_dofs(FACE) > 0){
		if(dim >= FACE) constrained_face_dof_indices<TBaseElem,ConstrainingTriangle, Face, Face>(elem,fct,ind,vFace);
		if(dim >= FACE) constrained_face_dof_indices<TBaseElem,ConstrainingQuadrilateral, Face, Face>(elem,fct,ind,vFace);
	}
		
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
			const size_t numDoFsOnSub = num_fct_dofs(fct,subRoid,si);

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
                                const typename Grid::traits<TSubBaseElem>::secure_container& vElem) const
{
//	storage for offsets
	std::vector<size_t> vOrientOffset;

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
			const size_t numDoFsOnSub = num_fct_dofs(fct,subRoid,si);

		//	Orientation is required: Thus, we compute the offsets, that are
		//	no longer in the usual order [0, 1, 2, ...]. Orientation is
		//	required if there are more than 1 dof on a subelement of a
		//	finite element and thus, when gluing two elements together,
		//	also the dofs on the subelements have to fit in order to
		//	guarantee continuity. This is not needed for Vertices, since there
		//	no distinction can be made when all dofs are at the same position.
		//	This is also not needed for the highest dimension of a finite
		//	element, since the dofs on this geometric object must not be
		//	identified with other dofs.
			ComputeOrientationOffset(vOrientOffset, elem, subElem, i, lfeid(fct));

			UG_ASSERT(vOrientOffset.size() == numDoFsOnSub ||
			          vOrientOffset.empty(), "Something wrong with orientation");

			if(!m_bGrouped)
			{
				const size_t index = obj_index(subElem) + offset(subRoid,si,fct);

				if(vOrientOffset.empty()){
					for(size_t j = 0; j < numDoFsOnSub; ++j)
						ind.push_back_index(fct, index + j);
				}else {
					for(size_t j = 0; j < numDoFsOnSub; ++j)
						ind.push_back_index(fct, index + vOrientOffset[j]);
				}
			}
			else
			{
			//	compute index
				const size_t index = obj_index(subElem);
				const size_t comp = offset(subRoid,si,fct);

				if(vOrientOffset.empty()){
					for(size_t j = 0; j < numDoFsOnSub; ++j)
						ind.push_back_multi_index(fct, index, comp + j);
				}else{
					for(size_t j = 0; j < numDoFsOnSub; ++j)
						ind.push_back_multi_index(fct, index, comp + vOrientOffset[j]);
				}
			}
		} // end loop functions
	} // end loop subelement

}

template <typename TConstraining, typename TConstrained, typename TBaseElem>
void MGDoFDistribution::
constrained_vertex_indices(LocalIndices& ind,
                    const typename Grid::traits<TBaseElem>::secure_container& vSubElem) const
{
//	loop all edges
	for(size_t i = 0; i < vSubElem.size(); ++i)
	{
	//	only constraining objects are of interest
		TConstraining* constrainingObj = dynamic_cast<TConstraining*>(vSubElem[i]);
		if(constrainingObj == NULL) continue;

	//	loop constraining vertices
		for(size_t j = 0; j != constrainingObj->num_constrained_vertices(); ++j)
		{
		//	get vertex
			TConstrained* vrt = constrainingObj->constrained_vertex(j);

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

// sort constrained edges by association to reference object vertices of coarse edge
template <typename TBaseElem,typename TConstraining, typename TConstrained>
void MGDoFDistribution::
sort_constrained_edges(std::vector<size_t>& sortedInd,TBaseElem* elem,TConstraining* constrainingObj,size_t objIndex) const
{
	static const int dim = TBaseElem::dim;
	ReferenceObjectID roid = elem->reference_object_id();
	const DimReferenceElement<dim>& refElem
		= ReferenceElementProvider::get<dim>(roid);
	// get edge belonging to reference id vertex 0 on edge
	const size_t vertexIndex = refElem.id(1,objIndex,0,0);
	sortedInd.resize(2);
	VertexBase* vertex0 = NULL;
	// get child of vertex
	if (dim==2){
		Face* baseElem = dynamic_cast<Face*>(elem);
		vertex0 = multi_grid()->template get_child<VertexBase,VertexBase>(baseElem->vertex(vertexIndex),0);
	}
	if (dim==3){
		Volume* baseElem = dynamic_cast<Volume*>(elem);
		vertex0 = multi_grid()->template get_child<VertexBase,VertexBase>(baseElem->vertex(vertexIndex),0);
	}
	TConstrained* edg = constrainingObj->constrained_edge(0);
	bool found = false;
	for (size_t k=0;k<2;k++){
		VertexBase* vrt = edg->vertex(k);
		if (vrt==vertex0){
			found = true;
			break;
		}
	}
	if (found==true){
		sortedInd[0]=0;
		sortedInd[1]=1;
	} else {
		sortedInd[0]=1;
		sortedInd[1]=0;
		// check
		bool found2 = false;
		TConstrained* otherEdge = constrainingObj->constrained_edge(1);
		for (size_t k=0;k<2;k++){
			VertexBase* vrt = otherEdge->vertex(k);
			if (vrt==vertex0){
				found2 = true;
				break;
			}
		}
		if (found2==false) UG_THROW("no edge found belonging to vertex 0\n");
	}
}

// sort constrained faces as given by association to reference object vertices of coarse faces
template <typename TBaseElem,typename TConstraining, typename TConstrained>
void MGDoFDistribution::
sort_constrained_faces(std::vector<size_t>& sortedInd,TBaseElem* elem,TConstraining* constrainingObj,size_t objIndex) const
{
	static const int dim = TBaseElem::dim;
	ReferenceObjectID roid = elem->reference_object_id();
	const DimReferenceElement<dim>& refElem
			= ReferenceElementProvider::get<dim>(roid);
	const size_t numVrt = constrainingObj->num_vertices();
	sortedInd.resize(4);
	VertexBase* vrt = NULL;
	Volume* baseElem = dynamic_cast<Volume*>(elem);
	for (size_t i=0;i<numVrt;i++){
		const size_t vertexIndex = refElem.id(2,objIndex,0,i);
		vrt = multi_grid()->template get_child<VertexBase,VertexBase>(baseElem->vertex(vertexIndex),0);
		// loop constrained faces to find face corresponding to vertex
		bool found = false;
		for (size_t j=0;j<numVrt;j++){
			TConstrained* face = constrainingObj->constrained_face(j);
			for (size_t k=0;k<face->num_vertices();k++){
				if (face->vertex(k)==vrt){
					found = true;
					sortedInd[i] = j;
					break;
				}
			}
		}
		if (found==false) UG_THROW("corresponding constrained object vertex not found");
	}
	// for triangle refinement inner face is still missing, it should be constrained_face(3)
	if (numVrt==3){ 
		// check if it is not face 3
		for (size_t i=0;i<3;i++) if (sortedInd[i]==3) {
			bool found = false;
			for (size_t j=0;j<3;j++){
				for (size_t k=0;k<3;k++){
					if (sortedInd[k]==j){ 
						found = true;
						break;
					}
				}
				if (found==false){ 
					sortedInd[3]=j;
					return;
				}
			}
		}
		sortedInd[3]=3;
	}
}

template <typename TBaseElem,typename TConstraining, typename TConstrained, typename TSubElem>
void MGDoFDistribution::
constrained_edge_indices(TBaseElem* elem,LocalIndices& ind,
                    const typename Grid::traits<TSubElem>::secure_container& vSubElem) const
{
	//	loop all edges
	for(size_t i = 0; i < vSubElem.size(); ++i)
	{
	//	only constraining objects are of interest
		TConstraining* constrainingObj = dynamic_cast<TConstraining*>(vSubElem[i]);
		if(constrainingObj == NULL) continue;

		std::vector<size_t> sortedInd;
		sort_constrained_edges<TBaseElem,TConstraining,TConstrained>(sortedInd,elem,constrainingObj,i);

	//  loop constraining edges
		for(size_t j = 0; j != constrainingObj->num_constrained_edges(); ++j)
		{
			//	get edge
			TConstrained* edg = constrainingObj->constrained_edge(sortedInd[j]);

			//	get roid
			const ReferenceObjectID subRoid = edg->reference_object_id();

			//	get subset index
			int si = m_spMGSH->get_subset_index(edg);

			//	loop functions
			for(size_t fct = 0; fct < num_fct(); ++fct)
			{
			//	check that function is defined on subset
				if(!is_def_in_subset(fct, si)) continue;

				if(!m_bGrouped)
				{
				//	compute index
					const size_t index = obj_index(edg) + offset(subRoid,si,fct);

				//	add dof to local indices
					ind.push_back_index(fct, index);
				}
				else
				{
				//	compute index
					const size_t index = obj_index(edg);
					const size_t comp = offset(subRoid,si,fct);

				//	add dof to local indices
					ind.push_back_multi_index(fct, index, comp);
				}
			}
		}
	}
}

template <typename TBaseElem,typename TConstraining, typename TConstrained, typename TSubElem>
void MGDoFDistribution::
constrained_face_indices(TBaseElem* elem,LocalIndices& ind,
                    const typename Grid::traits<TSubElem>::secure_container& vSubElem) const
{
	//	loop all faces
	for(size_t i = 0; i < vSubElem.size(); ++i)
	{
	//	only constraining objects are of interest
		TConstraining* constrainingObj = dynamic_cast<TConstraining*>(vSubElem[i]);
		
		if(constrainingObj == NULL) continue;
		
		std::vector<size_t> sortedInd;
		sort_constrained_faces<TBaseElem,TConstraining,TConstrained>(sortedInd,elem,constrainingObj,i);

	//  loop constraining faces
		for(size_t j = 0; j != constrainingObj->num_constrained_faces(); ++j)
		{
			//	get face
			TConstrained* face = constrainingObj->constrained_face(sortedInd[j]);

			//	get roid
			const ReferenceObjectID subRoid = face->reference_object_id();

			//	get subset index
			int si = m_spMGSH->get_subset_index(face);

			//	loop functions
			for(size_t fct = 0; fct < num_fct(); ++fct)
			{
			//	check that function is defined on subset
				if(!is_def_in_subset(fct, si)) continue;

				if(!m_bGrouped)
				{
				//	compute index
					const size_t index = obj_index(face) + offset(subRoid,si,fct);

				//	add dof to local indices
					ind.push_back_index(fct, index);
				}
				else
				{
				//	compute index
					const size_t index = obj_index(face);
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

//	storage for (maybe needed) subelements
	Grid::SecureVertexContainer vCorner;
	Grid::SecureEdgeContainer vEdge;
	Grid::SecureFaceContainer vFace;
	Grid::SecureVolumeContainer vVol;

//	collect elements, if needed
	if(dim >= VERTEX)
		if(max_dofs(VERTEX) > 0) m_pMG->associated_elements_sorted(vCorner, elem);
	if(dim >= EDGE)
		if(max_dofs(EDGE) > 0 || bHang) m_pMG->associated_elements_sorted(vEdge, elem);
	if(dim >= FACE)
		if(max_dofs(FACE) > 0 || bHang) m_pMG->associated_elements_sorted(vFace, elem);
	if(dim >= VOLUME)
		if(max_dofs(VOLUME) > 0) m_pMG->associated_elements_sorted(vVol, elem);

//	get reference object id
	const ReferenceObjectID roid = elem->reference_object_id();

//	get regular dofs on all subelements and the element itself
//	use specialized function for vertices (since only one position and one reference object)
	if(dim >= VERTEX && max_dofs(VERTEX) > 0) indices_on_vertex<TBaseElem>(elem, roid, ind, vCorner);
	if(dim >= EDGE && max_dofs(EDGE) > 0) 	  indices<TBaseElem, EdgeBase>(elem, roid, ind, vEdge);
	if(dim >= FACE && max_dofs(FACE) > 0) 	  indices<TBaseElem, Face>(elem, roid, ind, vFace);
	if(dim >= VOLUME && max_dofs(VOLUME) > 0) indices<TBaseElem, Volume>(elem, roid, ind, vVol);

//	If no hanging dofs are required, we're done
	if(!bHang) return;

//	get dofs on hanging vertices
	if(max_dofs(VERTEX > 0))
	{
		if(dim >= EDGE) constrained_vertex_indices<ConstrainingEdge, VertexBase, EdgeBase>(ind, vEdge);
		if(dim >= FACE) constrained_vertex_indices<ConstrainingQuadrilateral, VertexBase, Face>(ind, vFace);
	}

//	get dofs on hanging edges
	if (max_dofs(EDGE) > 0){
		if(dim >= EDGE) constrained_edge_indices<TBaseElem,ConstrainingEdge, EdgeBase, EdgeBase>(elem,ind, vEdge);
		if(dim >= FACE) constrained_edge_indices<TBaseElem,ConstrainingTriangle, EdgeBase, Face>(elem,ind, vFace);
		if(dim >= FACE) constrained_edge_indices<TBaseElem,ConstrainingQuadrilateral, EdgeBase, Face>(elem,ind, vFace);
	}

//  get dofs on hanging faces
	if (max_dofs(FACE) > 0){
		if(dim >= FACE) constrained_face_indices<TBaseElem,ConstrainingTriangle, Face, Face>(elem,ind, vFace);
		if(dim >= FACE) constrained_face_indices<TBaseElem,ConstrainingQuadrilateral, Face, Face>(elem,ind, vFace);
	}

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

size_t MGDoFDistribution::dof_indices(GeometricObject* elem, size_t fct,
                                        std::vector<DoFIndex>& ind,
                                        bool bHang, bool bClear) const
{
	switch(elem->base_object_id())
	{
		case VERTEX: return dof_indices(static_cast<VertexBase*>(elem), fct, ind, bHang, bClear);
		case EDGE: return dof_indices(static_cast<EdgeBase*>(elem), fct, ind, bHang, bClear);
		case FACE: return dof_indices(static_cast<Face*>(elem), fct, ind, bHang, bClear);
		case VOLUME: return dof_indices(static_cast<Volume*>(elem), fct, ind, bHang, bClear);
		default: UG_THROW("Geometric Base element not found.");
	}
}

size_t MGDoFDistribution::inner_dof_indices(GeometricObject* elem, size_t fct,
                                              std::vector<DoFIndex>& ind,
                                              bool bClear) const
{
	switch(elem->base_object_id())
	{
		case VERTEX: return inner_dof_indices(static_cast<VertexBase*>(elem), fct, ind, bClear);
		case EDGE: return inner_dof_indices(static_cast<EdgeBase*>(elem), fct, ind, bClear);
		case FACE: return inner_dof_indices(static_cast<Face*>(elem), fct, ind, bClear);
		case VOLUME: return inner_dof_indices(static_cast<Volume*>(elem), fct, ind, bClear);
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

bool MGDoFDistribution::add(GeometricObject* elem, const ReferenceObjectID roid,
                            const int si, LevInfo& li)
{
	switch(elem->base_object_id())
	{
		case VERTEX: return add(static_cast<VertexBase*>(elem), roid, si, li);
		case EDGE: return add(static_cast<EdgeBase*>(elem), roid, si, li);
		case FACE: return add(static_cast<Face*>(elem), roid, si, li);
		case VOLUME: return add(static_cast<Volume*>(elem), roid, si, li);
		default: UG_THROW("Geometric Base element not found.");
	}
	return false;
}

///////////////////////////////////////////////////////////////////////////////
// template instantiations
///////////////////////////////////////////////////////////////////////////////
template void MGDoFDistribution::indices<VertexBase>(VertexBase*, LocalIndices&, bool) const;
template void MGDoFDistribution::indices<EdgeBase>(EdgeBase*, LocalIndices&, bool) const;
template void MGDoFDistribution::indices<Face>(Face*, LocalIndices&, bool) const;
template void MGDoFDistribution::indices<Volume>(Volume*, LocalIndices&, bool) const;

template size_t MGDoFDistribution::dof_indices<VertexBase>(VertexBase*, size_t, std::vector<DoFIndex>&, bool, bool) const;
template size_t MGDoFDistribution::dof_indices<EdgeBase>(EdgeBase*, size_t, std::vector<DoFIndex>&, bool, bool) const;
template size_t MGDoFDistribution::dof_indices<Face>(Face*, size_t, std::vector<DoFIndex>&, bool, bool) const;
template size_t MGDoFDistribution::dof_indices<Volume>(Volume*, size_t, std::vector<DoFIndex>&, bool, bool) const;

template size_t MGDoFDistribution::inner_dof_indices<VertexBase>(VertexBase*, size_t,std::vector<DoFIndex>&, bool) const;
template size_t MGDoFDistribution::inner_dof_indices<EdgeBase>(EdgeBase*, size_t,std::vector<DoFIndex>&, bool) const;
template size_t MGDoFDistribution::inner_dof_indices<Face>(Face*, size_t,std::vector<DoFIndex>&, bool) const;
template size_t MGDoFDistribution::inner_dof_indices<Volume>(Volume*, size_t,std::vector<DoFIndex>&, bool) const;

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
