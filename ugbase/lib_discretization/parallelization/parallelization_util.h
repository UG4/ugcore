/*
 * parallelization.h
 *
 *  Created on: 21.5.2010
 *      Author: A. Vogel, S.Reiter
 */

#ifndef __H__LIB_DISCRETIZATION__PARALLELIZATION__PARALLELIZATION_UTIL__
#define __H__LIB_DISCRETIZATION__PARALLELIZATION__PARALLELIZATION_UTIL__


#include "lib_grid/parallelization/parallelization.h"

#include "lib_algebra/parallelization/parallel_index_layout.h"
#include "lib_algebra/parallelization/communication_policies.h"

namespace ug
{

///	Adds dof-indices of elements in elemLayout to the specified IndexLayout.
/**
 * Make sure that TLayout holds elements of type VertexBase*, EdgeBase*,
 * Face* or Volume*.
 *
 * \param ignoreElemsWithStatus: an or combination of ug::ElementStatus.
 * 			Elements sharing at least on status with this param are not
 * 			added to the layout.
 * 			This param only has effect if dGrMgr is specified.
 *
 *  \todo: replace IndexLayout with TDoFManager::IndexLayout.
 */
template <class TDoFDistr, class TLayout>
bool AddEntriesToLevelIndexLayout(IndexLayout& indexLayoutOut,
                                  TDoFDistr& dofDistr,
                                  TLayout& elemLayout,
                                  DistributedGridManager* dGrMgr = NULL,
                                  byte ignoreElemsWithStatus = ES_NONE)
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

		//	check whether the element also lies in a vertical interface.
		//	If so, no horizontal interfaces should be created.
			if(dGrMgr){
				if(dGrMgr->get_status(elem) & ignoreElemsWithStatus){
					continue;
				}
			}

		//	get the algebraic indices on the grid element
			typename TDoFDistr::algebra_index_vector_type indices;
			dofDistr.get_inner_algebra_indices(elem, indices);

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
 * \param ignoreElemsWithStatus: an or combination of ug::ElementStatus.
 * 			Elements sharing at least on status with this param are not
 * 			added to the layout (default ES_NONE).
 * 			This param only has effect if dGrMgr is specified.
 */
template <class TDoFDistribution>
bool CreateLevelIndexLayout(	IndexLayout& layoutOut,
                            	TDoFDistribution& dofDistr,
                            	GridLayoutMap& layoutMap,
                            	int keyType, int level,
                            	DistributedGridManager* dGrMgr = NULL,
                            	byte ignoreElemsWithStatus = ES_NONE)
{
//	clear the layout
	layoutOut.clear();

//	success flag
	bool bRetVal = true;

// 	add dofs on elements
	if(dofDistr.template has_dofs_on<VertexBase>())
		if(layoutMap.has_layout<VertexBase>(keyType))
		{
			bRetVal &= AddEntriesToLevelIndexLayout(layoutOut, dofDistr,
									layoutMap.get_layout<VertexBase>(keyType).layout_on_level(level),
									dGrMgr, ignoreElemsWithStatus);
		}

	if(dofDistr.template has_dofs_on<EdgeBase>())
		if(layoutMap.has_layout<EdgeBase>(keyType))
		{
			bRetVal &= AddEntriesToLevelIndexLayout(layoutOut, dofDistr,
									layoutMap.get_layout<EdgeBase>(keyType).layout_on_level(level),
									dGrMgr, ignoreElemsWithStatus);
		}

	if(dofDistr.template has_dofs_on<Face>())
		if(layoutMap.has_layout<Face>(keyType))
		{
			bRetVal &= AddEntriesToLevelIndexLayout(layoutOut, dofDistr,
									layoutMap.get_layout<Face>(keyType).layout_on_level(level),
									dGrMgr, ignoreElemsWithStatus);
		}

	if(dofDistr.template has_dofs_on<Volume>())
		if(layoutMap.has_layout<Volume>(keyType))
		{
			bRetVal &= AddEntriesToLevelIndexLayout(layoutOut, dofDistr,
									layoutMap.get_layout<Volume>(keyType).layout_on_level(level),
									dGrMgr, ignoreElemsWithStatus);
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
template <class TDoFDistr, class TLayout>
bool AddEntriesToSurfaceIndexLayout(IndexLayout& indexLayoutOut,
                                    TDoFDistr& dofDistr,
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
			typename TDoFDistr::algebra_index_vector_type indices;
			dofDistr.get_inner_algebra_indices(elem, indices);

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
template <class TDoFDistribution>
bool CreateSurfaceIndexLayout(	IndexLayout& layoutOut,
                            	TDoFDistribution& dofDistr,
                            	GridLayoutMap& layoutMap,
                            	int keyType,
                            	MultiGrid& mg, DistributedGridManager& dGrMgr)
{
//	success flag
	bool bRetVal = true;

// 	add dofs on elements
	if(dofDistr.template has_dofs_on<VertexBase>())
		for(size_t level = 0; level < layoutMap.get_layout<VertexBase>(keyType).num_levels(); ++level)
			if(layoutMap.has_layout<VertexBase>(keyType))
			{
				bRetVal &= AddEntriesToSurfaceIndexLayout(layoutOut, dofDistr,
										layoutMap.get_layout<VertexBase>(keyType).layout_on_level(level), mg, dGrMgr);
			}

	if(dofDistr.template has_dofs_on<EdgeBase>())
		for(size_t level = 0; level < layoutMap.get_layout<EdgeBase>(keyType).num_levels(); ++level)
			if(layoutMap.has_layout<EdgeBase>(keyType))
			{
				bRetVal &= AddEntriesToSurfaceIndexLayout(layoutOut, dofDistr,
										layoutMap.get_layout<EdgeBase>(keyType).layout_on_level(level), mg, dGrMgr);
			}

	if(dofDistr.template has_dofs_on<Face>())
		for(size_t level = 0; level < layoutMap.get_layout<Face>(keyType).num_levels(); ++level)
			if(layoutMap.has_layout<Face>(keyType))
			{
				bRetVal &= AddEntriesToSurfaceIndexLayout(layoutOut, dofDistr,
										layoutMap.get_layout<Face>(keyType).layout_on_level(level), mg, dGrMgr);
			}

	if(dofDistr.template has_dofs_on<Volume>())
		for(size_t level = 0; level < layoutMap.get_layout<Volume>(keyType).num_levels(); ++level)
			if(layoutMap.has_layout<Volume>(keyType))
			{
				bRetVal &= AddEntriesToSurfaceIndexLayout(layoutOut, dofDistr,
										layoutMap.get_layout<Volume>(keyType).layout_on_level(level), mg, dGrMgr);
			}

//	we're done
	return bRetVal;
}



template <typename TMatrix, typename TDoFDistr>
void CopyLayoutsAndCommunicatorIntoMatrix(TMatrix& mat, TDoFDistr& dofDistr)
{
	mat.set_layouts(dofDistr.get_master_layout(), dofDistr.get_slave_layout());

	mat.set_communicator(dofDistr.get_communicator());
	mat.set_process_communicator(dofDistr.get_process_communicator());
}


/**
 *
 */
template <class TDoFDistr, class TLayout>
bool AddEntriesToIndexLayout_DomainDecomposition(
							IndexLayout& processLayoutOut,
							IndexLayout& subdomainLayoutOut,
							TDoFDistr& dofDistr,
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
				typename TDoFDistr::algebra_index_vector_type indices;
				dofDistr.get_inner_algebra_indices(elem, indices);
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
				typename TDoFDistr::algebra_index_vector_type indices;
				dofDistr.get_inner_algebra_indices(elem, indices);
				for(size_t i = 0; i < indices.size(); ++i)
				{
					indexInterface.push_back(indices[i]);
				}
			}
		}
	}
	return true;
}


template <class TDoFDistribution>
bool CreateIndexLayouts_DomainDecomposition(
						IndexLayout& processLayoutOut,
						IndexLayout& subdomainLayoutOut,
						TDoFDistribution& dofDistr,
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
