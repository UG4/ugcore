/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
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

#ifndef __H__PCL__pcl_layout_util__
#define __H__PCL__pcl_layout_util__

#include <vector>

#include "common/util/hash.h"

#include "pcl_communication_structs.h"


namespace pcl{

/// \addtogroup pcl
/// \{

///	removes all empty interfaces from the given layout.
template <typename TLayout>
void RemoveEmptyInterfaces(TLayout& layout)
{
	using TInterfaceIter = typename TLayout::iterator;
	using TInterface = typename TLayout::Interface;

	for(size_t lvl = 0; lvl < layout.num_levels(); ++lvl){
		for(TInterfaceIter iter = layout.begin(lvl); iter != layout.end(lvl);)
		{
			TInterface& intfc = layout.interface(iter);
			if(intfc.empty())
				iter = layout.erase(iter, lvl);
			else
				++iter;
		}
	}
}

////////////////////////////////////////////////////////////////////////
///	collects the ids of all processes to which interfaces exist.
/**
 * Fills a vector with the process-ids, to which interfaces exist in
 * the given layout.
 *
 * TLayout has to be compatible with pcl::Layout or pcl::MultiLevelLayout.
 *
 * \returns the number of associated processes.
 */
template <typename TLayout>
size_t CollectAssociatedProcesses(std::vector<int>& procIDsOut,
								  TLayout& layout)
{
	procIDsOut.clear();

//	iterate through the levels of the layout
	for(size_t i = 0; i < layout.num_levels(); ++i){
	//	iterate through the interfaces on that level
		for(typename TLayout::iterator iIter = layout.begin(i); iIter != layout.end(i); ++iIter)
		{
			int procID = layout.proc_id(iIter);
		//	check whether the process is already contained in procIDsOut
			if(i > 0){
				if(find(procIDsOut.begin(), procIDsOut.end(), procID)
				   == procIDsOut.end())
				{
				//	the entry has not yet been added
					procIDsOut.push_back(procID);
				}
			}
			else{
			//	on level 0 each process exists only once
				procIDsOut.push_back(procID);
			}
		}
	}

	return procIDsOut.size();
}

///	writes all elements in the interfaces into the vector.
/**
 * This function extracts all elements from a layout and stores them into
 * a std::vector. Doubles may occur and are not removed. The container is
 * clear as default, before extracting.
 */
template <typename TLayout>
void CollectElements(std::vector<typename TLayout::Element>& elemsOut,
					TLayout& layout,
					bool clearContainer = true)
{
	using Interface = typename TLayout::Interface;

//	clear the return value
	if(clearContainer) elemsOut.clear();

//	iterate over all interfaces
	for(size_t lvl = 0; lvl < layout.num_levels(); ++lvl){
		for(typename TLayout::const_iterator interfaceIter = layout.begin(lvl); interfaceIter != layout.end(lvl); ++interfaceIter)
		{
		//	iterate over the entries of the interface
			const Interface& interface = layout.interface(interfaceIter);
			for(typename Interface::const_iterator iter = interface.begin(); iter != interface.end(); ++iter)
			{
			//	add elem to vector
				elemsOut.push_back(interface.get_element(iter));
			}
		}
	}
}

///	writes all elements in the interfaces into the resulting vector. avoids doubles.
template <typename TLayout>
void CollectUniqueElements(std::vector<typename TLayout::Element>& elemsOut,
						   const TLayout& layout)
{
	using Interface = typename TLayout::Interface;
	using TElem = typename TLayout::Element;

//	clear the return value
	elemsOut.clear();

//	we'll use a hash to make sure that each element only exists once
	ug::Hash<TElem, int> hash(layout.num_interface_elements());
	hash.reserve(layout.num_interface_elements());

//	iterate over all interfaces
	for(size_t lvl = 0; lvl < layout.num_levels(); ++lvl){
		for(typename TLayout::const_iterator interfaceIter = layout.begin(lvl); interfaceIter != layout.end(lvl); ++interfaceIter)
		{
		//	iterate over the entries of the interface
			const Interface& interface = layout.interface(interfaceIter);
			for(typename Interface::const_iterator iter = interface.begin(); iter != interface.end(); ++iter)
			{
			//	check whether the entry already exists in the hash
				if(!hash.has_entry(interface.get_element(iter))){
				//	we don't care about the value
					hash.insert(interface.get_element(iter), 0);
					elemsOut.push_back(interface.get_element(iter));
				}
			}
		}
	}
}

// end group pcl
/// \}

}// end of namespace

#endif
