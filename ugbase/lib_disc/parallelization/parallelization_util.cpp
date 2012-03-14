/*
 * parallelization.cpp
 *
 *  Created on: 24.03.2011
 *      Author: andreasvogel, S.Reiter
 */

#include "parallelization_util.h"

namespace ug{

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

	int localProc = pcl::GetProcRank();

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

			if(vIndexAppear.size() > 0)
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

} // end namespace ug
