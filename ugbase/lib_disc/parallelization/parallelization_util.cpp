/*
 * parallelization.cpp
 *
 *  Created on: 24.03.2011
 *      Author: andreasvogel, S.Reiter
 */

#include "parallelization_util.h"

namespace ug{


//	Adds dof-indices of elements in elemLayout to the specified IndexLayout.
/*
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
template <class TLayout>
bool AddEntriesToLevelIndexLayout(IndexLayout& indexLayoutOut,
                                  DoFDistribution& dofDistr, TLayout& elemLayout,
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

bool CreateLevelIndexLayout(	IndexLayout& layoutOut,
                                   	DoFDistribution& dofDistr,
                                   	GridLayoutMap& layoutMap,
                                   	int keyType, int level)
{
//	clear the layout
	layoutOut.clear();

//	success flag
	bool bRetVal = true;

// 	add dofs on elements
	if(dofDistr.max_dofs(VERTEX))
		if(layoutMap.has_layout<VertexBase>(keyType))
		{
			bRetVal &= AddEntriesToLevelIndexLayout(layoutOut, dofDistr,
									layoutMap.get_layout<VertexBase>(keyType).layout_on_level(level));
		}

	if(dofDistr.max_dofs(EDGE))
		if(layoutMap.has_layout<EdgeBase>(keyType))
		{
			bRetVal &= AddEntriesToLevelIndexLayout(layoutOut, dofDistr,
									layoutMap.get_layout<EdgeBase>(keyType).layout_on_level(level));
		}

	if(dofDistr.max_dofs(FACE))
		if(layoutMap.has_layout<Face>(keyType))
		{
			bRetVal &= AddEntriesToLevelIndexLayout(layoutOut, dofDistr,
									layoutMap.get_layout<Face>(keyType).layout_on_level(level));
		}

	if(dofDistr.max_dofs(VOLUME))
		if(layoutMap.has_layout<Volume>(keyType))
		{
			bRetVal &= AddEntriesToLevelIndexLayout(layoutOut, dofDistr,
									layoutMap.get_layout<Volume>(keyType).layout_on_level(level));
		}

//	we're done
	return bRetVal;
}


//	Adds dof-indices of elements in elemLayout to the specified IndexLayout.
/*
 * Make sure that TLayout holds elements of type VertexBase*, EdgeBase*,
 * Face* or Volume*.
 *
 *  \todo: replace IndexLayout with TDoFManager::IndexLayout.
 *
 * \param[in]		mg				underlying MultiGrid
 */
template <class TLayout>
bool AddEntriesToSurfaceIndexLayout(IndexLayout& indexLayoutOut,
                                    DoFDistribution& dofDistr,
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

bool CreateSurfaceIndexLayout(	IndexLayout& layoutOut,
                            	DoFDistribution& dofDistr,
                            	GridLayoutMap& layoutMap,
                            	int keyType,
                            	MultiGrid& mg, DistributedGridManager& dGrMgr)
{
//	success flag
	bool bRetVal = true;

// 	add dofs on elements
	if(dofDistr.max_dofs(VERTEX)){
		if(layoutMap.has_layout<VertexBase>(keyType)){
			for(size_t level = 0; level < layoutMap.get_layout<VertexBase>(keyType).num_levels(); ++level)
			{
				bRetVal &= AddEntriesToSurfaceIndexLayout(layoutOut, dofDistr,
										layoutMap.get_layout<VertexBase>(keyType).layout_on_level(level), mg, dGrMgr);
			}
		}
	}

	if(dofDistr.max_dofs(EDGE)){
		if(layoutMap.has_layout<EdgeBase>(keyType)){
			for(size_t level = 0; level < layoutMap.get_layout<EdgeBase>(keyType).num_levels(); ++level)
			{
				bRetVal &= AddEntriesToSurfaceIndexLayout(layoutOut, dofDistr,
										layoutMap.get_layout<EdgeBase>(keyType).layout_on_level(level), mg, dGrMgr);
			}
		}
	}

	if(dofDistr.max_dofs(FACE)){
		if(layoutMap.has_layout<Face>(keyType)){
			for(size_t level = 0; level < layoutMap.get_layout<Face>(keyType).num_levels(); ++level)
			{
				bRetVal &= AddEntriesToSurfaceIndexLayout(layoutOut, dofDistr,
										layoutMap.get_layout<Face>(keyType).layout_on_level(level), mg, dGrMgr);
			}
		}
	}

	if(dofDistr.max_dofs(VOLUME)){
		if(layoutMap.has_layout<Volume>(keyType)){
			for(size_t level = 0; level < layoutMap.get_layout<Volume>(keyType).num_levels(); ++level)
			{
				bRetVal &= AddEntriesToSurfaceIndexLayout(layoutOut, dofDistr,
										layoutMap.get_layout<Volume>(keyType).layout_on_level(level), mg, dGrMgr);
			}
		}
	}

//	we're done
	return bRetVal;
}


/// returns in a vector all appearencies of an index in a layout
void FindPositionInInterfaces(std::vector<std::pair<int, size_t> >& vIndexInterface,
                                     IndexLayout& layout, size_t index)
{
	for(IndexLayout::iterator interface_iter = layout.begin();
				interface_iter != layout.end(); ++interface_iter)
	{
	//	get interface
		IndexLayout::Interface& interface = layout.interface(interface_iter);

		int targetProc   = layout.proc_id(interface_iter);

	//	loop over indices
		int i = 0;
		for( IndexLayout::Interface::iterator iter = interface.begin();
				iter != interface.end(); ++iter, ++i)
		{
			size_t currIndex = interface.get_element(iter);

			if(currIndex == index)
				vIndexInterface.push_back(std::pair<int, size_t>(targetProc, i));
		}
	}
}

bool AddExtraProcessEntriesToSubdomainLayout(
								size_t numIDs,
								IndexLayout& processMasterLayoutIn,
								IndexLayout& processSlaveLayoutIn,
								IndexLayout& subdomainMasterLayoutInOut,
								IndexLayout& subdomainSlaveLayoutInOut)
{
	std::vector<int> vMultiplicity;
//	generate an id for each entry.
	vMultiplicity.clear();
	vMultiplicity.resize(numIDs, 0);

	int localProc = pcl::ProcRank();

	for(IndexLayout::iterator interface_iter = processMasterLayoutIn.begin();
			interface_iter != processMasterLayoutIn.end(); ++interface_iter)
	{
	//	get interface
		IndexLayout::Interface& interface = processMasterLayoutIn.interface(interface_iter);
		int targetProc = processMasterLayoutIn.proc_id(interface_iter);

	//	loop over indices
		for( IndexLayout::Interface::iterator iter = interface.begin();
				iter != interface.end(); ++iter)
		{
		//  get index
			const size_t index = interface.get_element(iter);

			UG_LOG("Checking index " << index << std::endl);

			std::vector<std::pair<int, size_t> > vIndexAppear;
			FindPositionInInterfaces(vIndexAppear, subdomainMasterLayoutInOut, index);

			if(!vIndexAppear.empty())
			{
			//	 flag
				vMultiplicity[index] = 1;

				UG_LOG("Flagging index=" << index << " on Proc " << localProc << " to target proc "<< targetProc<<std::endl);

			//	add to subdomain interface
			//	get interface
				IndexLayout::Interface& subdomInterface = subdomainMasterLayoutInOut.interface(targetProc);

				subdomInterface.push_back(index);
			}
		}
	}

//	Communicate flagged vector to slaves
	ComPol_VecCopy<std::vector<int> >	copyPol(&vMultiplicity);

	pcl::InterfaceCommunicator<IndexLayout> communicator;
	communicator.send_data(processMasterLayoutIn, copyPol);
	communicator.receive_data(processSlaveLayoutIn, copyPol);
	communicator.communicate();

//	add slaves
	for(IndexLayout::iterator interface_iter = processSlaveLayoutIn.begin();
			interface_iter != processSlaveLayoutIn.end(); ++interface_iter)
	{
	//	get interface
		IndexLayout::Interface& interface = processSlaveLayoutIn.interface(interface_iter);
		int targetProc = processSlaveLayoutIn.proc_id(interface_iter);

	//	loop over indices
		for( IndexLayout::Interface::iterator iter = interface.begin();
				iter != interface.end(); ++iter)
		{
		//  get index
			const size_t index = interface.get_element(iter);

		//	is flagged?
			if(vMultiplicity[index] > 0)
			{
			//	add to subdomain interface
			//	get interface
				IndexLayout::Interface& subdomInterface = subdomainSlaveLayoutInOut.interface(targetProc);

				UG_LOG("Adding index " << index << " on Proc "<< localProc<<" to target Proc interface "<<targetProc << "\n");

				subdomInterface.push_back(index);
			}
		}
	}

	return true;
}


// permutes an IndexLayout for the permutation of indices
void PermuteIndicesInIndexLayout(	IndexLayout& layout,
									const std::vector<size_t>& vIndNew)
{
	typedef IndexLayout::Interface 	Interface;
	typedef IndexLayout::iterator 	InterfaceIter;
	typedef Interface::iterator		IndexIter;
	typedef Interface::Element		Index;

//	iterate over all interfaces
	for(InterfaceIter iiter = layout.begin();
		iiter != layout.end(); ++iiter)
	{
	//	get interface
		Interface& interface = layout.interface(iiter);

	//	iterate over all elements
		for(IndexIter indIter = interface.begin();
			indIter != interface.end(); ++indIter)
		{
		//	get old index
			Index oldIndex = interface.get_element(indIter);

		//	check index
			UG_ASSERT(oldIndex < vIndNew.size(), "Invalid index.");

		//	replace by new index
			interface.get_element(indIter) = vIndNew[oldIndex];
		}
	}
}



/**
 *
 */
template <class TLayout>
bool AddEntriesToIndexLayout_DomainDecomposition(
							IndexLayout& processLayoutOut,
							IndexLayout& subdomainLayoutOut,
							DoFDistribution& dofDistr,
							TLayout& elemLayout,
							pcl::IDomainDecompositionInfo* ddInfoIn)
{
	typedef typename TLayout::iterator InterfaceIterator;
	typedef typename TLayout::Interface ElemInterface;
	typedef typename ElemInterface::iterator ElemIterator;

	typedef IndexLayout::Interface IndexInterface;

	int localProc = pcl::ProcRank();
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


bool CreateIndexLayouts_DomainDecomposition(
						IndexLayout& processLayoutOut,
						IndexLayout& subdomainLayoutOut,
						DoFDistribution& dofDistr,
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

} // end namespace ug
