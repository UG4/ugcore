/*
 * parallelization.h
 *
 *  Created on: 21.5.2010
 *      Author: A. Vogel, S.Reiter
 */

#ifndef __H__UG__LIB_DISC__PARALLELIZATION__PARALLELIZATION_UTIL__
#define __H__UG__LIB_DISC__PARALLELIZATION__PARALLELIZATION_UTIL__

#include "pcl/pcl_util.h"
#include "lib_grid/parallelization/parallelization.h"
#include "lib_grid/parallelization/util/compol_interface_status.h"
#include "lib_algebra/parallelization/parallel_index_layout.h"
#include "lib_algebra/parallelization/communication_policies.h"
#include "lib_algebra/parallelization/parallel_vector.h"
#include "lib_algebra/parallelization/parallel_matrix.h"

namespace ug
{

///	Adds dof-indices of elements in elemLayout to the specified IndexLayout.
/**
 * Make sure that TLayout holds elements of type VertexBase*, EdgeBase*,
 * Face* or Volume*.
 *
 * \param pIgnoreMap	If specified (defaul: NULL), this map will be used to
 * 						check whether an element shall not be added to the new
 * 						interface. For each interface in the given layout, the
 * 						map has to hold a vector<bool> of the size of the interface.
 * 						If an entry is true, then the corresponding interface-
 * 						entry will be ignored.
 *
 *  \todo: replace IndexLayout with TDoFManager::IndexLayout.
 */
template <class TDD, class TLayout>
bool AddEntriesToLevelIndexLayout(IndexLayout& indexLayoutOut,
					  TDD& dofDistr, TLayout& elemLayout,
					  const std::map<int, std::vector<bool> >* pIgnoreMap = NULL)
{
//	iterator for grid element interfaces
	typedef typename TLayout::iterator InterfaceIterator;

//	type of grid element interfaces
	typedef typename TLayout::Interface ElemInterface;

//	iterator for grid elements
	typedef typename ElemInterface::iterator ElemIterator;

//	type of index interfaces
	typedef IndexLayout::Interface IndexInterface;

//	iterate over all grid element interfaces
	for(InterfaceIterator iIter = elemLayout.begin();
		iIter != elemLayout.end(); ++iIter)
	{
	//	get a grid element interface
		ElemInterface& elemInterface = elemLayout.interface(iIter);

	//	get a corresponding index interface
		IndexInterface& indexInterface = indexLayoutOut.interface(
											elemLayout.proc_id(iIter));

	//	if some elements shall be ignored, then we'll perform a special loop
		if(pIgnoreMap){
			std::map<int, std::vector<bool> >::const_iterator
				findIter = pIgnoreMap->find(elemInterface.get_target_proc());

			UG_ASSERT(findIter != pIgnoreMap->end(), "The vector has to exist");
			const std::vector<bool>& vec = findIter->second;

			UG_ASSERT(vec.size() == elemInterface.size(), "Sizes have to match!");

		//	iterate over entries in the grid element interface
			int counter = 0;
			for(ElemIterator eIter = elemInterface.begin();
				eIter != elemInterface.end(); ++eIter, ++counter)
			{
			//	if the corresponding vec-entry is true, then we'll ignore the elem.
				if(vec[counter])
					continue;

			//	get the grid element
				typename ElemInterface::Element elem = elemInterface.get_element(eIter);

			//	get the algebraic indices on the grid element
				std::vector<size_t> indices;
				dofDistr.inner_algebra_indices(elem, indices);

			//	add the indices to the interface
				for(size_t i = 0; i < indices.size(); ++i)
					indexInterface.push_back(indices[i]);
			}
		}
		else{
		//	iterate over entries in the grid element interface
			for(ElemIterator eIter = elemInterface.begin();
				eIter != elemInterface.end(); ++eIter)
			{
			//	get the grid element
				typename ElemInterface::Element elem = elemInterface.get_element(eIter);

			//	get the algebraic indices on the grid element
				std::vector<size_t> indices;
				dofDistr.inner_algebra_indices(elem, indices);

			//	add the indices to the interface
				for(size_t i = 0; i < indices.size(); ++i)
					indexInterface.push_back(indices[i]);
			}
		}
	}

//	touching an interface means creation. Thus we remove the empty interfaces
//	to avoid storage, communication (should not happen any longer) etc...
	pcl::RemoveEmptyInterfaces(elemLayout);

//	we're done
	return true;
}

/// creates the index layout for a level given a GridLayoutMap
/**
 * This function creates the Index layout based on a GridLayout. All elements
 * of the GridLayoutMap are loop on grid level and the indices attached to the
 * grid elements are added to the interface. Since the ordering of the grid
 * elements in the interfaces is assumed to be correct, also the ordering in
 * the index layouts are correct.
 *
 * \param[out]		layoutOut		the created index layout
 * \param[in]		dofDistr		the DoF Distribution
 * \param[in]		layoutMap		the grid layout map
 * \param[in]		keyType			key type (e.g. slave or master)
 * \param[in]		level			level, where layouts should be build
 *
 */
template <class TDD>
bool CreateLevelIndexLayout(	IndexLayout& layoutOut,
                            	TDD& dofDistr,
                            	GridLayoutMap& layoutMap,
                            	int keyType, int level)
{
//	clear the layout
	layoutOut.clear();

//	success flag
	bool bRetVal = true;

// 	add dofs on elements
	if(dofDistr.has_indices_on(VERTEX))
		if(layoutMap.has_layout<VertexBase>(keyType))
		{
			bRetVal &= AddEntriesToLevelIndexLayout(layoutOut, dofDistr,
									layoutMap.get_layout<VertexBase>(keyType).layout_on_level(level));
		}

	if(dofDistr.has_indices_on(EDGE))
		if(layoutMap.has_layout<EdgeBase>(keyType))
		{
			bRetVal &= AddEntriesToLevelIndexLayout(layoutOut, dofDistr,
									layoutMap.get_layout<EdgeBase>(keyType).layout_on_level(level));
		}

	if(dofDistr.has_indices_on(FACE))
		if(layoutMap.has_layout<Face>(keyType))
		{
			bRetVal &= AddEntriesToLevelIndexLayout(layoutOut, dofDistr,
									layoutMap.get_layout<Face>(keyType).layout_on_level(level));
		}

	if(dofDistr.has_indices_on(VOLUME))
		if(layoutMap.has_layout<Volume>(keyType))
		{
			bRetVal &= AddEntriesToLevelIndexLayout(layoutOut, dofDistr,
									layoutMap.get_layout<Volume>(keyType).layout_on_level(level));
		}

//	we're done
	return bRetVal;
}


///	Adds dof-indices of elements in elemLayout to the specified IndexLayout.
/**
 * Make sure that TLayout holds elements of type VertexBase*, EdgeBase*,
 * Face* or Volume*.
 *
 *  \todo: replace IndexLayout with TDoFManager::IndexLayout.
 *
 * \param[in]		mg				underlying MultiGrid
 */
template <class TDD, class TLayout>
bool AddEntriesToSurfaceIndexLayout(IndexLayout& indexLayoutOut,
                                    TDD& dofDistr,
                                    TLayout& elemLayout,
                                    MultiGrid& mg,
                                    DistributedGridManager& dGrMgr)
{
//	iterator for grid element interfaces
	typedef typename TLayout::iterator InterfaceIterator;

//	type of grid element interfaces
	typedef typename TLayout::Interface ElemInterface;

//	iterator for grid elements
	typedef typename ElemInterface::iterator ElemIterator;

//	type of index interfaces
	typedef IndexLayout::Interface IndexInterface;

//	iterate over all grid element interfaces
	for(InterfaceIterator iIter = elemLayout.begin();
		iIter != elemLayout.end(); ++iIter)
	{
	//	get a grid element interface
		ElemInterface& elemInterface = elemLayout.interface(iIter);

	//	get a corresponding index interface
		IndexInterface& indexInterface = indexLayoutOut.interface(
											elemLayout.proc_id(iIter));

	//	iterate over entries in the grid element interface
		for(ElemIterator eIter = elemInterface.begin();
			eIter != elemInterface.end(); ++eIter)
		{
		//	get the grid element
			typename ElemInterface::Element elem = elemInterface.get_element(eIter);

		//	check if element is on surface (i.e. has no children). Shadows are
		//	not taken into account here, since their indices are already added
		//	to the interface by the shadowing objects
			if(mg.has_children(elem)) {continue;}

		//	check if element is a ghost element, i.e. it is a surface element
		//	but only due to a hierarchical cut of the grid in order to
		//	refine it further on another process. These cuts lead to so called
		//	vertical interfaces.
			if(dGrMgr.is_ghost(elem)) {continue;}

		//	get the algebraic indices on the grid element
			std::vector<size_t> indices;
			dofDistr.inner_algebra_indices(elem, indices);

		//	add the indices to the interface
			for(size_t i = 0; i < indices.size(); ++i)
			{
				indexInterface.push_back(indices[i]);
			}
		}
	}

//	touching an interface means creation. Thus we remove the empty interfaces
//	to avoid storage, communication (should not happen any longer) etc...
	pcl::RemoveEmptyInterfaces(elemLayout);

//	we're done
	return true;
}

/// creates the index layout for a level given a GridLayoutMap
/**
 * This function creates the Index layout based on a GridLayout. All elements
 * of the GridLayoutMap are loop level by level and the indices attached to the
 * grid elements are added to the interface, if an element does not have a children.
 * Since the ordering of the grid elements in the interfaces is assumed to be
 * correct, also the ordering in the index layouts are correct.
 *
 * \param[out]		layoutOut		the created index layout
 * \param[in]		dofDistr		the DoF Distribution
 * \param[in]		layoutMap		the grid layout map
 * \param[in]		keyType			key type (e.g. slave or master)
 * \param[in]		mg				underlying MultiGrid
 * \param[in]		dGrMgr			distributed Grid Manager
 */
template <class TDD>
bool CreateSurfaceIndexLayout(	IndexLayout& layoutOut,
                            	TDD& dofDistr,
                            	GridLayoutMap& layoutMap,
                            	int keyType,
                            	MultiGrid& mg, DistributedGridManager& dGrMgr)
{
//	success flag
	bool bRetVal = true;

// 	add dofs on elements
	if(dofDistr.has_indices_on(VERTEX))
		for(size_t level = 0; level < layoutMap.get_layout<VertexBase>(keyType).num_levels(); ++level)
			if(layoutMap.has_layout<VertexBase>(keyType))
			{
				bRetVal &= AddEntriesToSurfaceIndexLayout(layoutOut, dofDistr,
										layoutMap.get_layout<VertexBase>(keyType).layout_on_level(level), mg, dGrMgr);
			}

	if(dofDistr.has_indices_on(EDGE))
		for(size_t level = 0; level < layoutMap.get_layout<EdgeBase>(keyType).num_levels(); ++level)
			if(layoutMap.has_layout<EdgeBase>(keyType))
			{
				bRetVal &= AddEntriesToSurfaceIndexLayout(layoutOut, dofDistr,
										layoutMap.get_layout<EdgeBase>(keyType).layout_on_level(level), mg, dGrMgr);
			}

	if(dofDistr.has_indices_on(FACE))
		for(size_t level = 0; level < layoutMap.get_layout<Face>(keyType).num_levels(); ++level)
			if(layoutMap.has_layout<Face>(keyType))
			{
				bRetVal &= AddEntriesToSurfaceIndexLayout(layoutOut, dofDistr,
										layoutMap.get_layout<Face>(keyType).layout_on_level(level), mg, dGrMgr);
			}

	if(dofDistr.has_indices_on(VOLUME))
		for(size_t level = 0; level < layoutMap.get_layout<Volume>(keyType).num_levels(); ++level)
			if(layoutMap.has_layout<Volume>(keyType))
			{
				bRetVal &= AddEntriesToSurfaceIndexLayout(layoutOut, dofDistr,
										layoutMap.get_layout<Volume>(keyType).layout_on_level(level), mg, dGrMgr);
			}

//	we're done
	return bRetVal;
}


/// copies all needed parallel informations into a parallel matrix
/**
 * This function copies the layouts and the communicators of a dof distribution
 * into a matrix.
 *
 * \param[in]	dd		DoFDistribution to copy the infos from
 * \param[out]	mat		Matrix filled with infos
 *
 * \tparam	TMatrix 	Sequential Matrix type
 */
template <typename TMatrix, typename TDD>
void CopyLayoutsAndCommunicatorIntoMatrix(ParallelMatrix<TMatrix>& mat,
                                          TDD& dd)
{
	mat.set_layouts(dd.master_layout(), dd.slave_layout());

	mat.set_communicator(dd.communicator());
	mat.set_process_communicator(dd.process_communicator());
}

/// copies all needed parallel informations into a parallel vector
/**
 * This function copies the layouts and the communicators of a dof distribution
 * into a vector.
 *
 * \param[in]	dd		DoFDistribution to copy the infos from
 * \param[out]	vec		Vector filled with infos
 *
 * \tparam 	TVector		Sequential vector type
 */
template <typename TVector, typename TDD>
void CopyLayoutsAndCommunicatorIntoVector(ParallelVector<TVector>& vec,
                                          TDD& dd)
{
	//	copy all horizontal layouts (for all domain decomps)
		vec.set_layouts(dd.master_layout(), dd.slave_layout());

	//	copy vertical layouts
		vec.set_vertical_layouts(dd.vertical_master_layout(),
		                         dd.vertical_slave_layout());

	//	copy communicator
		vec.set_communicator(dd.communicator());
		vec.set_process_communicator(dd.process_communicator());
}


/**
 *
 */
template <class TDD, class TLayout>
bool AddEntriesToIndexLayout_DomainDecomposition(
							IndexLayout& processLayoutOut,
							IndexLayout& subdomainLayoutOut,
							TDD& dofDistr,
							TLayout& elemLayout,
							pcl::IDomainDecompositionInfo* ddInfoIn)
{
	typedef typename TLayout::iterator InterfaceIterator;
	typedef typename TLayout::Interface ElemInterface;
	typedef typename ElemInterface::iterator ElemIterator;

	typedef IndexLayout::Interface IndexInterface;

	int localProc = pcl::GetProcRank();
	int localSubdom = ddInfoIn->map_proc_id_to_subdomain_id(localProc);

//	iterate over all interfaces
	for(InterfaceIterator iIter = elemLayout.begin();
		iIter != elemLayout.end(); ++iIter)
	{
		ElemInterface& elemInterface = elemLayout.interface(iIter);
		int targetProc = elemLayout.proc_id(iIter);
		int targetSubdom = ddInfoIn->map_proc_id_to_subdomain_id(targetProc);

		if(targetSubdom == localSubdom){
		//	create a process interface
			IndexInterface& indexInterface = processLayoutOut.interface(targetProc);

		//	iterate over entries in the elemInterface and add associated
		//	dofs to the indexInterface
			for(ElemIterator eIter = elemInterface.begin();
				eIter != elemInterface.end(); ++eIter)
			{
				typename ElemInterface::Element elem = elemInterface.get_element(eIter);
				std::vector<size_t> indices;
				dofDistr.inner_algebra_indices(elem, indices);
				for(size_t i = 0; i < indices.size(); ++i)
				{
					indexInterface.push_back(indices[i]);
				}
			}
		}
		else{
		//	create a subdomain interface
			IndexInterface& indexInterface = subdomainLayoutOut.interface(targetProc);

		//	iterate over entries in the elemInterface and add associated
		//	dofs to the indexInterface
			for(ElemIterator eIter = elemInterface.begin();
				eIter != elemInterface.end(); ++eIter)
			{
				typename ElemInterface::Element elem = elemInterface.get_element(eIter);
				std::vector<size_t> indices;
				dofDistr.inner_algebra_indices(elem, indices);
				for(size_t i = 0; i < indices.size(); ++i)
				{
					indexInterface.push_back(indices[i]);
				}
			}
		}
	}
	return true;
}


template <class TDD>
bool CreateIndexLayouts_DomainDecomposition(
						IndexLayout& processLayoutOut,
						IndexLayout& subdomainLayoutOut,
						TDD& dofDistr,
						GridLayoutMap& layoutMap,
						int keyType, int level,
						pcl::IDomainDecompositionInfo* ddInfoIn)
{
//TODO: clear the layout!
	bool bRetVal = true;
	if(layoutMap.has_layout<VertexBase>(keyType)){
		bRetVal &= AddEntriesToIndexLayout_DomainDecomposition(
								processLayoutOut,
								subdomainLayoutOut,
								dofDistr,
								layoutMap.get_layout<VertexBase>(keyType).
								layout_on_level(level),
								ddInfoIn); /*(cb_ProcIDToSubdomID)*/
	}
/*
	if(layoutMap.has_layout<EdgeBase>(keyType)){
		bRetVal &= AddEntriesToIndexLayout(layoutOut, dofManager,
								layoutMap.get_layout<EdgeBase>(keyType).layout_on_level(level));
	}
	if(layoutMap.has_layout<Face>(keyType)){
		bRetVal &= AddEntriesToIndexLayout(layoutOut, dofManager,
								layoutMap.get_layout<Face>(keyType).layout_on_level(level));
	}
	if(layoutMap.has_layout<Volume>(keyType)){
		bRetVal &= AddEntriesToIndexLayout(layoutOut, dofManager,
								layoutMap.get_layout<Volume>(keyType).layout_on_level(level));
	}
*/
	return bRetVal;
}

// returns in a vector all appearencies of an index in a layout
void FindPositionInInterfaces(std::vector<std::pair<int, size_t> >& vIndexInterface,
                                     IndexLayout& layout, size_t index);

bool AddExtraProcessEntriesToSubdomainLayout(
								size_t numIDs,
								IndexLayout& processMasterLayoutIn,
								IndexLayout& processSlaveLayoutIn,
								IndexLayout& subdomainMasterLayoutInOut,
								IndexLayout& subdomainSlaveLayoutInOut);

/// permutes an IndexLayout for the permutation of indices
/**
 * This Function changes the indices in the layout according to a given
 * permutation of the indices. (The order of the DoFs in the interfaces remains
 * the same, but the DoFs are "renamed")
 * The vector vIndNew must return the new index for each old index,
 * i.e. newIndex = vIndNew[oldIndex].
 *
 * \param[in]	vIndNew		mapping for each index
 * \returns 	success flag
 */
void PermuteIndicesInIndexLayout(	IndexLayout& layout,
									const std::vector<size_t>& vIndNew);

}//	end of namespace

#endif
