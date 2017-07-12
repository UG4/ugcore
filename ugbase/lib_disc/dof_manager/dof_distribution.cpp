/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#include "dof_distribution.h"
#include "lib_disc/function_spaces/grid_function.h"

#include "common/log.h"
#include "lib_disc/domain.h"
#include "lib_disc/local_finite_element/local_dof_set.h"
#include "lib_disc/reference_element/reference_element.h"
#include "lib_disc/reference_element/reference_element_traits.h"
#include "lib_disc/reference_element/reference_mapping.h"
#include "lib_disc/reference_element/reference_mapping_provider.h"
#include "lib_disc/common/groups_util.h"
#include "common/util/string_util.h"
#include "orientation.h"
#include "lib_grid/tools/periodic_boundary_manager.h"

#include "lib_grid/file_io/file_io.h"
#include "lib_grid/algorithms/debug_util.h"

using namespace std;

namespace ug{

////////////////////////////////////////////////////////////////////////////////
// DoFDistribution
////////////////////////////////////////////////////////////////////////////////

DoFDistribution::
DoFDistribution(SmartPtr<MultiGrid> spMG,
                SmartPtr<MGSubsetHandler> spMGSH,
                ConstSmartPtr<DoFDistributionInfo> spDDInfo,
                SmartPtr<SurfaceView> spSurfView,
                const GridLevel& level, bool bGrouped,
                SmartPtr<DoFIndexStorage> spDoFIndexStorage)
	: DoFDistributionInfoProvider(spDDInfo),
      m_bGrouped(bGrouped),
	  m_spMG(spMG),
	  m_pMG(m_spMG.get()),
	  m_spMGSH(spMGSH),
	  m_spSurfView(spSurfView),
	  m_gridLevel(level),
	  m_spDoFIndexStorage(spDoFIndexStorage),
	  m_numIndex(0)
{
	if(m_spDoFIndexStorage.invalid())
		m_spDoFIndexStorage = SmartPtr<DoFIndexStorage>(new DoFIndexStorage(spMG, spDDInfo));

	check_subsets();
	m_vNumIndexOnSubset.resize(num_subsets(), 0);

#ifdef UG_PARALLEL
	m_spAlgebraLayouts = SmartPtr<AlgebraLayouts>(new AlgebraLayouts);
#endif

	reinit();
}


DoFDistribution::
~DoFDistribution() {}


void DoFDistribution::check_subsets()
{
//	check, that all geom objects are assigned to a subset
	if(	m_spMGSH->num<Vertex>() != multi_grid()->num<Vertex>())
		UG_THROW("All Vertices "
			   " must be assigned to a subset. The passed subset handler "
			   " contains non-assigned elements, thus the dof distribution"
			   " is not possible, aborting.");

	if(	m_spMGSH->num<Edge>() != multi_grid()->num<Edge>())
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

SurfaceView::SurfaceConstants DoFDistribution::defaultValidSurfState() const{
	if(m_gridLevel.is_level()) return SurfaceView::ALL;
	else if(m_gridLevel.is_surface()) return SurfaceView::ALL_BUT_SHADOW_COPY;
	else UG_THROW("DoFDistribution: GridLevel::type not valid.")
}

////////////////////////////////////////////////////////////////////////////////
// DoFDistribution: Index Access
////////////////////////////////////////////////////////////////////////////////

template <typename TBaseObject>
void DoFDistribution::
add(TBaseObject* obj, const ReferenceObjectID roid, const int si)
{
	UG_ASSERT(si >= 0, "Invalid subset index passed");

//	if no dofs on this subset for the roid, do nothing
	if(num_dofs(roid,si) == 0) return;

//	check for periodicity
	bool master = false;
	if(m_spMG->has_periodic_boundaries()){
		PeriodicBoundaryManager& pbm = *m_spMG->periodic_boundary_manager();
		if(pbm.is_slave(obj)) return; // ignore slaves
		if(pbm.is_master(obj)) master = true;
	}

//	compute the number of indices needed on the Geometric object
	size_t numNewIndex = 1;
	if(!m_bGrouped) numNewIndex = num_dofs(roid,si);

// 	set first available index to the object. The first available index is the
//	first managed index plus the size of the index set. (If holes are in the
//	index set, this is not treated here, holes remain)
	obj_index(obj) = m_numIndex;

//	number of managed indices and the number of managed indices on the subset has
//	changed. Thus, increase the counters.
	m_numIndex += numNewIndex;
	m_vNumIndexOnSubset[si] += numNewIndex;

// 	if obj is a master, assign all its slaves
	if(master) {
		typedef typename PeriodicBoundaryManager::Group<TBaseObject>::SlaveContainer SlaveContainer;
		typedef typename PeriodicBoundaryManager::Group<TBaseObject>::SlaveIterator SlaveIterator;
		SlaveContainer& slaves = *m_spMG->periodic_boundary_manager()->slaves(obj);
		size_t master_index = obj_index(obj);
		for(SlaveIterator iter = slaves.begin(); iter != slaves.end(); ++iter){
			obj_index(*iter) = master_index;
		}
	}
}

template <typename TBaseElem>
size_t DoFDistribution::
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
void DoFDistribution::
extract_inner_algebra_indices(const typename Grid::traits<TBaseElem>::secure_container& vElem,
                              std::vector<size_t>& ind) const
{
//	loop passed elements
	for(size_t i = 0; i < vElem.size(); ++i)
		inner_algebra_indices(vElem[i], ind, false);
}

template<typename TBaseElem>
size_t DoFDistribution::_inner_algebra_indices(TBaseElem* elem,
                                                std::vector<size_t>& ind,
                                                bool bClear) const
{
//	clear indices
	if(bClear) ind.clear();

//	return
	return extract_inner_algebra_indices(elem, ind);
}

// this function could be merged with inner_algebra_indices with additional
// default paramter e.g. selectedFct==-1 -> all functions.
size_t DoFDistribution::inner_algebra_indices_for_fct(GridObject* elem,
                                                std::vector<size_t>& ind,
                                                bool bClear, int fct) const
{
//	clear indices
	if(bClear) ind.clear();

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
		//	check that function is def on subset
			if(!is_def_in_subset(fct, si)) return ind.size();

		//	get number of DoFs in this sub-geometric object
			const size_t numDoFsOnSub = num_fct_dofs(fct,roid,si);

		//	compute index
			const size_t index = firstIndex + offset(roid,si,fct);

		//	add dof to local indices
			for(size_t j = 0; j < numDoFsOnSub; ++j)
				ind.push_back(index+j);
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
size_t DoFDistribution::_algebra_indices(TBaseElem* elem,
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
		extract_inner_algebra_indices<Vertex>(vVrt, ind);
	}
	if(dim >= EDGE && max_dofs(EDGE) > 0)
	{
		Grid::SecureEdgeContainer vEdge;
		m_pMG->associated_elements(vEdge, elem);
		extract_inner_algebra_indices<Edge>(vEdge, ind);
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
void DoFDistribution::
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
size_t DoFDistribution::_inner_dof_indices(TBaseElem* elem, size_t fct,
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
void DoFDistribution::
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
void DoFDistribution::
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
void DoFDistribution::
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
size_t DoFDistribution::_dof_indices(TBaseElem* elem, size_t fct,
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
	if(dim >= VERTEX && max_dofs(VERTEX) > 0) dof_indices<TBaseElem, Vertex>(elem, roid, fct, ind, vCorner);
	if(dim >= EDGE && max_dofs(EDGE) > 0) 	  dof_indices<TBaseElem, Edge>(elem, roid, fct, ind, vEdge);
	if(dim >= FACE && max_dofs(FACE) > 0) 	  dof_indices<TBaseElem, Face>(elem, roid, fct, ind, vFace);
	if(dim >= VOLUME && max_dofs(VOLUME) > 0) dof_indices<TBaseElem, Volume>(elem, roid, fct, ind, vVol);

//	If no hanging dofs are required, we're done
	if(!bHang) return ind.size();

	//	get dofs on hanging vertices
	if(max_dofs(VERTEX > 0))
	{
		if(dim >= EDGE) constrained_vertex_dof_indices<ConstrainingEdge, Vertex, Edge>(fct,ind,vEdge);
		if(dim >= FACE) constrained_vertex_dof_indices<ConstrainingQuadrilateral, Vertex, Face>(fct,ind,vFace);
	}

//	get dofs on hanging edges
	if (max_dofs(EDGE) > 0){
		if(dim >= EDGE) constrained_edge_dof_indices<TBaseElem,ConstrainingEdge, Edge, Edge>(elem,fct,ind, vEdge);
		if(dim >= FACE) constrained_edge_dof_indices<TBaseElem,ConstrainingTriangle, Edge, Face>(elem,fct,ind, vFace);
		if(dim >= FACE) constrained_edge_dof_indices<TBaseElem,ConstrainingQuadrilateral, Edge, Face>(elem,fct,ind, vFace);
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
void DoFDistribution::indices_on_vertex(TBaseElem* elem, const ReferenceObjectID roid,
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
void DoFDistribution::indices(TBaseElem* elem, const ReferenceObjectID roid,
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
void DoFDistribution::
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
void DoFDistribution::
sort_constrained_edges(std::vector<size_t>& sortedInd,TBaseElem* elem,TConstraining* constrainingObj,size_t objIndex) const
{
	static const int dim = TBaseElem::dim;
	ReferenceObjectID roid = elem->reference_object_id();
	const DimReferenceElement<dim>& refElem
		= ReferenceElementProvider::get<dim>(roid);
	// get edge belonging to reference id vertex 0 on edge
	const size_t vertexIndex = refElem.id(1,objIndex,0,0);
	sortedInd.resize(2);
	Vertex* vertex0 = NULL;
	// get child of vertex
	if (dim==2){
		Face* baseElem = dynamic_cast<Face*>(elem);
		vertex0 = multi_grid()->template get_child<Vertex,Vertex>(baseElem->vertex(vertexIndex),0);
	}
	if (dim==3){
		Volume* baseElem = dynamic_cast<Volume*>(elem);
		vertex0 = multi_grid()->template get_child<Vertex,Vertex>(baseElem->vertex(vertexIndex),0);
	}
	TConstrained* edg = constrainingObj->constrained_edge(0);
	bool found = false;
	for (size_t k=0;k<2;k++){
		Vertex* vrt = edg->vertex(k);
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
			Vertex* vrt = otherEdge->vertex(k);
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
void DoFDistribution::
sort_constrained_faces(std::vector<size_t>& sortedInd,TBaseElem* elem,TConstraining* constrainingObj,size_t objIndex) const
{
	static const int dim = TBaseElem::dim;
	ReferenceObjectID roid = elem->reference_object_id();
	const DimReferenceElement<dim>& refElem
			= ReferenceElementProvider::get<dim>(roid);
	const size_t numVrt = constrainingObj->num_vertices();
	sortedInd.resize(4);
	Vertex* vrt = NULL;
	Volume* baseElem = dynamic_cast<Volume*>(elem);
	for (size_t i=0;i<numVrt;i++){
		const size_t vertexIndex = refElem.id(2,objIndex,0,i);
		vrt = multi_grid()->template get_child<Vertex,Vertex>(baseElem->vertex(vertexIndex),0);
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
void DoFDistribution::
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
void DoFDistribution::
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

template<typename TBaseElem>
void DoFDistribution::_indices(TBaseElem* elem, LocalIndices& ind, bool bHang) const
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
	if(dim >= EDGE && max_dofs(EDGE) > 0) 	  indices<TBaseElem, Edge>(elem, roid, ind, vEdge);
	if(dim >= FACE && max_dofs(FACE) > 0) 	  indices<TBaseElem, Face>(elem, roid, ind, vFace);
	if(dim >= VOLUME && max_dofs(VOLUME) > 0) indices<TBaseElem, Volume>(elem, roid, ind, vVol);

//	If no hanging dofs are required, we're done
	if(!bHang) return;

//	get dofs on hanging vertices
	if (max_dofs(VERTEX) > 0)
	{
		if(dim >= EDGE) constrained_vertex_indices<ConstrainingEdge, Vertex, Edge>(ind, vEdge);
		if(dim >= FACE) constrained_vertex_indices<ConstrainingQuadrilateral, Vertex, Face>(ind, vFace);
	}

//	get dofs on hanging edges
	if (max_dofs(EDGE) > 0){
		if(dim >= EDGE) constrained_edge_indices<TBaseElem,ConstrainingEdge, Edge, Edge>(elem,ind, vEdge);
		if(dim >= FACE) constrained_edge_indices<TBaseElem,ConstrainingTriangle, Edge, Face>(elem,ind, vFace);
		if(dim >= FACE) constrained_edge_indices<TBaseElem,ConstrainingQuadrilateral, Edge, Face>(elem,ind, vFace);
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
void DoFDistribution::
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

		UG_ASSERT(adjInd < (size_t)1e10, "adjInd = " << adjInd); // <-- I get an adjInd of (size_t) (-1) here!

	//	add to index list
		vIndex.push_back(adjInd);
	}
}

///////////////////////////////////////////////////////////////////////////////
// forwarding fcts
///////////////////////////////////////////////////////////////////////////////

void DoFDistribution::indices(Vertex* elem, LocalIndices& ind, bool bHang) const {_indices<Vertex>(elem, ind, bHang);}
void DoFDistribution::indices(Edge* elem, LocalIndices& ind, bool bHang) const {_indices<Edge>(elem, ind, bHang);}
void DoFDistribution::indices(Face* elem, LocalIndices& ind, bool bHang) const {_indices<Face>(elem, ind, bHang);}
void DoFDistribution::indices(Volume* elem, LocalIndices& ind, bool bHang) const {_indices<Volume>(elem, ind, bHang);}
void DoFDistribution::indices(GridObject* elem, LocalIndices& ind, bool bHang) const
{
	switch(elem->base_object_id())
	{
		case VERTEX: return indices(static_cast<Vertex*>(elem), ind, bHang);
		case EDGE: return indices(static_cast<Edge*>(elem), ind, bHang);
		case FACE: return indices(static_cast<Face*>(elem), ind, bHang);
		case VOLUME: return indices(static_cast<Volume*>(elem), ind, bHang);
		default: UG_THROW("Geometric Base element not found.");
	}
}

size_t DoFDistribution::dof_indices(Vertex* elem, size_t fct, std::vector<DoFIndex>& ind, bool bHang, bool bClear) const {return _dof_indices<Vertex>(elem, fct, ind, bHang, bClear);}
size_t DoFDistribution::dof_indices(Edge* elem, size_t fct, std::vector<DoFIndex>& ind, bool bHang, bool bClear) const {return _dof_indices<Edge>(elem, fct, ind, bHang, bClear);}
size_t DoFDistribution::dof_indices(Face* elem, size_t fct, std::vector<DoFIndex>& ind, bool bHang, bool bClear) const {return _dof_indices<Face>(elem, fct, ind, bHang, bClear);}
size_t DoFDistribution::dof_indices(Volume* elem, size_t fct, std::vector<DoFIndex>& ind, bool bHang, bool bClear) const {return _dof_indices<Volume>(elem, fct, ind, bHang, bClear);}
size_t DoFDistribution::dof_indices(GridObject* elem, size_t fct, std::vector<DoFIndex>& ind, bool bHang, bool bClear) const
{
	switch(elem->base_object_id())
	{
		case VERTEX: return dof_indices(static_cast<Vertex*>(elem), fct, ind, bHang, bClear);
		case EDGE: return dof_indices(static_cast<Edge*>(elem), fct, ind, bHang, bClear);
		case FACE: return dof_indices(static_cast<Face*>(elem), fct, ind, bHang, bClear);
		case VOLUME: return dof_indices(static_cast<Volume*>(elem), fct, ind, bHang, bClear);
		default: UG_THROW("Geometric Base element not found.");
	}
}

size_t DoFDistribution::inner_dof_indices(Vertex* elem, size_t fct, std::vector<DoFIndex>& ind, bool bHang) const {return _inner_dof_indices<Vertex>(elem, fct, ind, bHang);}
size_t DoFDistribution::inner_dof_indices(Edge* elem, size_t fct, std::vector<DoFIndex>& ind, bool bHang) const {return _inner_dof_indices<Edge>(elem, fct, ind, bHang);}
size_t DoFDistribution::inner_dof_indices(Face* elem, size_t fct, std::vector<DoFIndex>& ind, bool bHang) const {return _inner_dof_indices<Face>(elem, fct, ind, bHang);}
size_t DoFDistribution::inner_dof_indices(Volume* elem, size_t fct, std::vector<DoFIndex>& ind, bool bHang) const {return _inner_dof_indices<Volume>(elem, fct, ind, bHang);}
size_t DoFDistribution::inner_dof_indices(GridObject* elem, size_t fct, std::vector<DoFIndex>& ind, bool bClear) const
{
	switch(elem->base_object_id())
	{
		case VERTEX: return inner_dof_indices(static_cast<Vertex*>(elem), fct, ind, bClear);
		case EDGE: return inner_dof_indices(static_cast<Edge*>(elem), fct, ind, bClear);
		case FACE: return inner_dof_indices(static_cast<Face*>(elem), fct, ind, bClear);
		case VOLUME: return inner_dof_indices(static_cast<Volume*>(elem), fct, ind, bClear);
		default: UG_THROW("Geometric Base element not found.");
	}
}


size_t DoFDistribution::algebra_indices(Vertex* elem, std::vector<size_t>& ind, bool bClear) const {return _algebra_indices<Vertex>(elem, ind, bClear);}
size_t DoFDistribution::algebra_indices(Edge* elem, std::vector<size_t>& ind, bool bClear) const {return _algebra_indices<Edge>(elem, ind, bClear);}
size_t DoFDistribution::algebra_indices(Face* elem, std::vector<size_t>& ind, bool bClear) const {return _algebra_indices<Face>(elem, ind, bClear);}
size_t DoFDistribution::algebra_indices(Volume* elem, std::vector<size_t>& ind, bool bClear) const {return _algebra_indices<Volume>(elem, ind, bClear);}
size_t DoFDistribution::algebra_indices(GridObject* elem,	std::vector<size_t>& ind, bool bClear) const
{
	switch(elem->base_object_id())
	{
		case VERTEX: return algebra_indices(static_cast<Vertex*>(elem), ind, bClear);
		case EDGE: return algebra_indices(static_cast<Edge*>(elem), ind, bClear);
		case FACE: return algebra_indices(static_cast<Face*>(elem), ind, bClear);
		case VOLUME: return algebra_indices(static_cast<Volume*>(elem), ind, bClear);
		default: UG_THROW("Geometric Base element not found.");
	}
}

size_t DoFDistribution::inner_algebra_indices(Vertex* elem, std::vector<size_t>& ind, bool bClear) const {return _inner_algebra_indices<Vertex>(elem, ind, bClear);}
size_t DoFDistribution::inner_algebra_indices(Edge* elem, std::vector<size_t>& ind, bool bClear) const {return _inner_algebra_indices<Edge>(elem, ind, bClear);}
size_t DoFDistribution::inner_algebra_indices(Face* elem, std::vector<size_t>& ind, bool bClear) const {return _inner_algebra_indices<Face>(elem, ind, bClear);}
size_t DoFDistribution::inner_algebra_indices(Volume* elem, std::vector<size_t>& ind, bool bClear) const {return _inner_algebra_indices<Volume>(elem, ind, bClear);}
size_t DoFDistribution::inner_algebra_indices(GridObject* elem, std::vector<size_t>& ind, bool bClear) const
{
	switch(elem->base_object_id())
	{
		case VERTEX: return inner_algebra_indices(static_cast<Vertex*>(elem), ind, bClear);
		case EDGE: return inner_algebra_indices(static_cast<Edge*>(elem), ind, bClear);
		case FACE: return inner_algebra_indices(static_cast<Face*>(elem), ind, bClear);
		case VOLUME: return inner_algebra_indices(static_cast<Volume*>(elem), ind, bClear);
		default: UG_THROW("Geometric Base element not found.");
	}
}


////////////////////////////////////////////////////////////////////////////////
// Resize management
////////////////////////////////////////////////////////////////////////////////

void DoFDistribution::manage_grid_function(IGridFunction& gridFct)
{
//	if already registered, we're done
	if(std::find(m_vpGridFunction.begin(), m_vpGridFunction.end(), &gridFct)
		!= m_vpGridFunction.end())
		return;

//	add to managed functions
	m_vpGridFunction.push_back(&gridFct);
}

void DoFDistribution::unmanage_grid_function(IGridFunction& gridFct)
{
	m_vpGridFunction.erase(
	std::remove(m_vpGridFunction.begin(), m_vpGridFunction.end(), &gridFct)
	, m_vpGridFunction.end());
}

void DoFDistribution::permute_values(const std::vector<size_t>& vIndNew)
{
//	swap values of handled grid functions
	for(size_t i = 0; i < m_vpGridFunction.size(); ++i)
		m_vpGridFunction[i]->permute_values(vIndNew);
}

void DoFDistribution::copy_values(const std::vector<std::pair<size_t, size_t> >& vIndexMap,
                                         bool bDisjunct)
{
//	swap values of handled grid functions
	for(size_t i = 0; i < m_vpGridFunction.size(); ++i)
		m_vpGridFunction[i]->copy_values(vIndexMap, bDisjunct);
}

void DoFDistribution::resize_values(size_t newSize)
{
//	swap values of handled grid functions
	for(size_t i = 0; i < m_vpGridFunction.size(); ++i)
		m_vpGridFunction[i]->resize_values(newSize);
}

////////////////////////////////////////////////////////////////////////////////
// Init DoFs
////////////////////////////////////////////////////////////////////////////////


template <typename TBaseElem>
void DoFDistribution::reinit()
{
	typedef typename traits<TBaseElem>::iterator iterator;
	static const int dim = TBaseElem::dim;

//	check if indices in the dimension
	if(max_dofs(dim) == 0) return;

//	SURFACE
	if(grid_level().type() == GridLevel::SURFACE){

	//	in order to also cater for some seldom occurring parallel cases, we have
	//	to perform a slightly cumbersome iteration here. The basic idea is the following:
	//	Each PURE_SURFACE element requires a dof. Furthermore all SHADOWING elements
	//	which are not SHADOWED require a dof. SHADOWED_NONCOPY elements also require
	//	a dof and SHADOWED_COPY elements have to copy their dofs from their SHADOWING
	//	element.
	//	In parallel, however, a problem occurs. The associated SHADOWING element of
	//	a SHADOW_COPY element does not necessarily exist on the same process
	//	(should occur on vertices only).
	//	We thus can't simply iterate over PURE_SURFACE | SHADOWING and pass values on
	//	to parents. Instead we have to iterate over all candidates and also give a
	//	dof to SHADOW_COPY elements which don't have children.

		const SurfaceView& sv = *m_spSurfView;
		MultiGrid& mg = *m_spMG;

		for(int si = 0; si < num_subsets(); ++si)
		{
		// 	skip if no dofs shall exist on the given subset
			if(max_dofs(dim, si) == 0) continue;

		//	iterate over all surface elements (including shadows)
			iterator iter = begin<TBaseElem>(si, SurfaceView::ALL);
			iterator iterEnd = end<TBaseElem>(si, SurfaceView::ALL);

			for(; iter != iterEnd; ++iter){
				TBaseElem* elem = *iter;
				if(sv.is_contained(elem, grid_level(), SurfaceView::SHADOW_RIM_COPY)){
					if(mg.num_children<TBaseElem>(elem) > 0){
						TBaseElem* child = mg.get_child<TBaseElem>(elem, 0);
						if(sv.is_contained(child, grid_level(), SurfaceView::SURFACE_RIM))
							continue;
					}
				}

			//	create a dof and copy it down to SHADOW_COPY parents
				const ReferenceObjectID roid = elem->reference_object_id();
				add(elem, roid, si);

				TBaseElem* p = dynamic_cast<TBaseElem*>(mg.get_parent(elem));
				while(p && sv.is_contained(p, grid_level(), SurfaceView::SHADOW_RIM_COPY)){
					obj_index(p) = obj_index(elem);
					p = dynamic_cast<TBaseElem*>(mg.get_parent(p));
				}
			}
		} // end subset
	}

	// LEVEL
	else if(grid_level().type() == GridLevel::LEVEL){

		for(int si = 0; si < num_subsets(); ++si)
		{
		// 	skip if no dofs to be distributed
			if(max_dofs(dim, si) == 0) continue;

		//	get iterators of elems
			iterator iter = begin<TBaseElem>(si);
			iterator iterEnd = end<TBaseElem>(si);

		// 	loop elems
			for(; iter != iterEnd; ++iter)
			{
			// 	get vertex
				TBaseElem* elem = *iter;

			//	get roid
				const ReferenceObjectID roid = elem->reference_object_id();

			//	add element
				add(elem, roid, si);
			}
		}

	}
	else{
		UG_THROW("DoFDistribution: GridLevel-Type"<<grid_level().type()<<" not supported");
	}
}

void DoFDistribution::reinit()
{
	m_numIndex = 0;
	m_vNumIndexOnSubset.resize(0);
	m_vNumIndexOnSubset.resize(num_subsets(), 0);

	if(max_dofs(VERTEX)) reinit<Vertex>();
	if(max_dofs(EDGE))   reinit<Edge>();
	if(max_dofs(FACE))   reinit<Face>();
	if(max_dofs(VOLUME)) reinit<Volume>();

#ifdef UG_PARALLEL
	reinit_layouts_and_communicator();
#endif
}


#ifdef UG_PARALLEL
void DoFDistribution::reinit_layouts_and_communicator()
{
	pcl::ProcessCommunicator commWorld;

//  -----------------------------------
//	CREATE PROCESS COMMUNICATOR
//  -----------------------------------
//	The idea  of local processes is to exclude processes from
//	e.g. norm computation that does not have a grid on a given
//	level. If no DoFs exist on the level, that level is excluded from
//	norm computations. In those cases the process votes false for
//	the subcommunicator.

// 	choose if this process participates
	bool participate = !commWorld.empty() && (num_indices() > 0);

//	create process communicator for interprocess layouts
	layouts()->proc_comm() = commWorld.create_sub_communicator(participate);

//  -----------------------------------
//	CREATE INDEX LAYOUTS ON LEVEL
//  -----------------------------------

	reinit_index_layout(layouts()->master(), INT_H_MASTER);
	reinit_index_layout(layouts()->slave(), INT_H_SLAVE);
	reinit_index_layout(layouts()->vertical_slave(), INT_V_SLAVE);

//	vertical layouts for ghosts
	if(grid_level().ghosts()){
		reinit_index_layout(layouts()->vertical_master(), INT_V_MASTER);
	}else{
		layouts()->vertical_master().clear();
	}
}

void DoFDistribution::reinit_index_layout(IndexLayout& layout, int keyType)
{
//	clear layout
	layout.clear();

//	add the index from grid layouts
	if(max_dofs(VERTEX)) add_indices_from_layouts<Vertex>(layout, keyType);
	if(max_dofs(EDGE))   add_indices_from_layouts<Edge>(layout, keyType);
	if(max_dofs(FACE))   add_indices_from_layouts<Face>(layout, keyType);
	if(max_dofs(VOLUME)) add_indices_from_layouts<Volume>(layout, keyType);

//	touching an interface means creation. Thus we remove the empty interfaces
//	to avoid storage, communication (should not happen any longer) etc...
	pcl::RemoveEmptyInterfaces(layout);
}

template <typename TBaseElem>
void DoFDistribution::
add_indices_from_layouts(IndexLayout& indexLayout,int keyType)
{
//	get the grid layout map
	GridLayoutMap& layoutMap = m_spMG->distributed_grid_manager()->grid_layout_map();

//	check if layout present
	if(!layoutMap.has_layout<TBaseElem>(keyType)) return;

//	types
	typedef typename GridLayoutMap::Types<TBaseElem>::Layout::LevelLayout TLayout;
	typedef typename TLayout::iterator InterfaceIterator;
	typedef typename TLayout::Interface ElemInterface;
	typedef typename ElemInterface::iterator ElemIterator;
	typedef IndexLayout::Interface IndexInterface;

//	vector for algebra indices
	std::vector<size_t> vIndex;

//	find levels to loop
	int lvlTo = (grid_level().top() ? (multi_grid()->num_levels()-1) : grid_level().level());
	int lvlFrom = (grid_level().is_surface() ? 0 : lvlTo);

//	loop all level
	for(int lvl = lvlFrom; lvl <= lvlTo; ++lvl)
	{
	//	get element layout
		TLayout& elemLayout = layoutMap.get_layout<TBaseElem>(keyType).layout_on_level(lvl);

	//	iterate over all grid element interfaces
		for(InterfaceIterator iIter = elemLayout.begin();
			iIter != elemLayout.end(); ++iIter)
		{
		//	get a grid element interface
			ElemInterface& elemInterface = elemLayout.interface(iIter);

		//	get a corresponding index interface
			IndexInterface& indexInterface = indexLayout.interface(elemLayout.proc_id(iIter));

		//	iterate over entries in the grid element interface
			for(ElemIterator eIter = elemInterface.begin();
				eIter != elemInterface.end(); ++eIter)
			{
			//	get the grid element
				typename ElemInterface::Element elem = elemInterface.get_element(eIter);

			//	check if element is a surface element
				// Using SURFACE will also admit SHADOW elements from a level distribution.
				// SHADOW_RIM elements must not be included; otherwise, any communication
				// between hMaster and hSlaves would be done twice (also by the corresponding
				// shadowing element, which references the same DoF), which is why we use
				// SURFACE here instead of, e.g., ALL.
				if(!m_spSurfView->is_contained(elem, grid_level(), SurfaceView::SURFACE))
					continue;

			//	add the indices to the interface
				inner_algebra_indices(elem, vIndex);
				for(size_t i = 0; i < vIndex.size(); ++i)
					indexInterface.push_back(vIndex[i]);
			}
		}
	}
}
#endif

////////////////////////////////////////////////////////////////////////////////
// DoF Statistic
////////////////////////////////////////////////////////////////////////////////


template <typename TBaseElem>
void DoFDistribution::sum_dof_count(DoFCount& cnt) const
{
	typedef typename traits<TBaseElem>::const_iterator iterator;

	const SurfaceView& sv = *m_spSurfView;
	DistributedGridManager* spDstGrdMgr = sv.subset_handler()->grid()->distributed_grid_manager();

//	get iterators of elems
	iterator iter = begin<TBaseElem>();
	iterator iterEnd = end<TBaseElem>();

// 	loop elems
	for(; iter != iterEnd; ++iter){
		TBaseElem* elem = *iter;
		const ReferenceObjectID roid = elem->reference_object_id();
		const int si = sv.subset_handler()->get_subset_index(elem);

		// Surface State
		SurfaceView::SurfaceState SurfaceState = sv.surface_state(elem, grid_level());

		// skip shadow-rim-copy: those are already numbered by their SURFACE_RIM
		if(SurfaceState.contains(SurfaceView::MG_SHADOW_RIM_COPY))
			continue;

		// Interface State
		byte InterfaceState = ES_NONE;
		if(spDstGrdMgr) InterfaceState = spDstGrdMgr->get_status(elem);

		// loop all functions
		for(size_t fct = 0; fct < num_fct(); ++fct){

			const size_t numDoF = num_fct_dofs(fct,roid,si);

			cnt.add(fct, si, SurfaceState, InterfaceState, numDoF);
		}
	}
}

DoFCount DoFDistribution::dof_count() const
{
	PROFILE_FUNC();
	DoFCount cnt(grid_level(), dof_distribution_info());

	if(max_dofs(VERTEX)) sum_dof_count<Vertex>(cnt);
	if(max_dofs(EDGE)) sum_dof_count<Edge>(cnt);
	if(max_dofs(FACE)) sum_dof_count<Face>(cnt);
	if(max_dofs(VOLUME)) sum_dof_count<Volume>(cnt);

	return cnt;
}


////////////////////////////////////////////////////////////////////////////////
// Permutation of indices
////////////////////////////////////////////////////////////////////////////////

template <typename TBaseElem>
void DoFDistribution::
get_connections(std::vector<std::vector<size_t> >& vvConnection) const
{
//	dimension of Base Elem
	static const int dim = TBaseElem::dim;

//	Adjacent geometric objects
	std::vector<Vertex*> vVrts;
	std::vector<Edge*> vEdges;
	std::vector<Face*> vFaces;
	std::vector<Volume*> vVols;

// 	Iterators
	typedef typename traits<TBaseElem>::const_iterator const_iterator;
	const_iterator iterEnd = end<TBaseElem>();
	const_iterator iter = begin<TBaseElem>();

//	loop elem
	for(; iter != iterEnd; ++iter)
	{
	// 	Get elem
		TBaseElem* elem = *iter;
		std::vector<size_t> vIndex;

	//	Get connected elements
		if(dim >= VERTEX && max_dofs(VERTEX) > 0) {
			collect_associated(vVrts, elem);
			changable_indices<Vertex>(vIndex, vVrts);
		}
		if(dim >= EDGE   && max_dofs(EDGE) > 0)	{
			collect_associated(vEdges, elem);
			changable_indices<Edge>(vIndex, vEdges);
		}
		if(dim >= FACE   && max_dofs(FACE) > 0)	{
			collect_associated(vFaces, elem);
			changable_indices<Face>(vIndex, vFaces);
		}
		if(dim >= VOLUME && max_dofs(VOLUME) > 0) {
			collect_associated(vVols, elem);
			changable_indices<Volume>(vIndex, vVols);
		}

	//	remove doubles
		std::sort(vIndex.begin(), vIndex.end());
		vIndex.erase(std::unique(vIndex.begin(), vIndex.end()), vIndex.end());

	//	add coupling to adjacency graph
		std::vector<size_t>::iterator it;

		for(size_t i = 0; i < vIndex.size(); ++i){
			for(size_t j = i; j < vIndex.size(); ++j)
			{
				const size_t iIndex = vIndex[i];
				const size_t jIndex = vIndex[j];

			//	add connection (iIndex->jIndex)
				it = find(vvConnection[iIndex].begin(), vvConnection[iIndex].end(), jIndex);
				if(it == vvConnection[iIndex].end()) vvConnection[iIndex].push_back(jIndex);


			//	add the opposite direction (adjInd -> index)
				it = find(vvConnection[jIndex].begin(), vvConnection[jIndex].end(), iIndex);
				if(it == vvConnection[jIndex].end()) vvConnection[jIndex].push_back(iIndex);
			}
		}
	}
}

void DoFDistribution::
get_connections(std::vector<std::vector<size_t> >& vvConnection) const
{
	vvConnection.clear();

//	if no subset given, we're done
	if(num_subsets() == 0) return;

//	check that in all subsets same number of functions and at least one
//	if this is not the case for non-grouped DoFs, we cannot allow reordering
	if(!grouped())
	{
		size_t numDoFs = 0;
		for(int si = 0; si < num_subsets(); ++si)
			for(int i = 0; i < NUM_REFERENCE_OBJECTS; ++i)
			{
				const ReferenceObjectID roid = (ReferenceObjectID) i;
				if(num_dofs(roid,si) == 0) continue;
				if(numDoFs == 0) {numDoFs = num_dofs(roid,si); continue;}
				if(num_dofs(roid,si) != numDoFs)
					UG_THROW("DoFDistribution::get_connections: "
							"Currently only implemented iff same number of DoFs on"
							" all geometric objects in all subsets: \n"
							"num_dofs("<<roid<<","<<si<<")="<<num_dofs(roid,si)<<
							", but previously found "<<numDoFs);
			}
	}

//	clear neighbors
	vvConnection.resize(num_indices());

	get_connections<Vertex>(vvConnection);
	get_connections<Edge>(vvConnection);
	get_connections<Face>(vvConnection);
	get_connections<Volume>(vvConnection);
}

template <typename TBaseElem>
void DoFDistribution::permute_indices(const std::vector<size_t>& vNewInd)
{
// 	loop Vertices
	typedef typename traits<TBaseElem>::const_iterator const_iterator;

	const_iterator iterEnd = end<TBaseElem>(SurfaceView::ALL);
	const_iterator iter = begin<TBaseElem>(SurfaceView::ALL);
	for(; iter != iterEnd; ++iter)
	{
	// 	get vertex
		TBaseElem* elem = *iter;

	// 	get current (old) index
		const size_t oldIndex = obj_index(elem);

	//	replace old index by new one
		obj_index(elem) = vNewInd[oldIndex];
	}
}

void DoFDistribution::permute_indices(const std::vector<size_t>& vNewInd)
{
	if(max_dofs(VERTEX)) permute_indices<Vertex>(vNewInd);
	if(max_dofs(EDGE))   permute_indices<Edge>(vNewInd);
	if(max_dofs(FACE))   permute_indices<Face>(vNewInd);
	if(max_dofs(VOLUME)) permute_indices<Volume>(vNewInd);

#ifdef UG_PARALLEL
	reinit_layouts_and_communicator();
#endif

//	permute indices in associated vectors
	permute_values(vNewInd);
}

} // end namespace ug
