/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__PCL_DEBUG__
#define __H__PCL_DEBUG__

#include <iostream>

#include "pcl_base.h"
#include "pcl_methods.h"
#include "pcl_communication_structs.h"
#include "pcl_interface_communicator.h"
#include "pcl_process_communicator.h"

////////////////////////////////////////////////////////////////////////
///	this allows us to print messages to the users terminal
/** Wrapper for UG_LOG.*/
#define PCLLOG(msg) UG_LOG(msg)


namespace pcl {

/// \addtogroup pcl
/// \{

////////////////////////////////////////////////////////////////////////
//	LogLayoutStructure
///	Logs the internals of a layout.
/**
 * Supported layouts are pcl::SingleLevelLayout and pcl::MultiLevelLayout.
 */
template <typename TLayout>
void LogLayoutStructure(TLayout& layout, const char* prefix = "")
{
	using namespace std;

	using Interface = typename TLayout::Interface;
	using InterfaceIter = typename TLayout::iterator;

	PCLLOG(prefix << "-- PCL_DEBUG: layout-structure --\n");
	PCLLOG(prefix << "---- num_levels: " << layout.num_levels() << endl);
	
	for(size_t lvl = 0; lvl < layout.num_levels(); ++lvl)
	{
		PCLLOG(prefix << "---- interfaces on level " << lvl << " (proc id, size): ");
				
		for(InterfaceIter iiter = layout.begin(lvl); iiter != layout.end(lvl); ++iiter)
		{
			Interface& interface = layout.interface(iiter);
			PCLLOG("(" << layout.proc_id(iiter) << ", " << interface.size() << "), ");
		}
		PCLLOG(endl);
	}
}

////////////////////////////////////////////////////////////////////////
//	LogLayoutMapStructure
///	Logs the internals of a layout-map for a given type.
/**
 * Supported layouts are pcl::SingleLevelLayout and pcl::MultiLevelLayout.
 */
template <typename TType, typename TLayoutMap>
void LogLayoutMapStructure(TLayoutMap& lm)
{
	using namespace std;

	using iterator = typename TLayoutMap::template Types<TType>::Map::iterator;
	using Layout = typename TLayoutMap::template Types<TType>::Layout;

	PCLLOG("-- PCL_DEBUG: layout-map-structure --\n");
	
	for(iterator iter = lm.template layouts_begin<TType>(); iter != lm.template layouts_end<TType>(); ++iter)
	{
	//	get the key
		PCLLOG("---- has key: " << iter->first << endl);
		
	//	log the layout structure
		Layout& layout = iter->second;
		LogLayoutStructure(layout, "----");
	}
}

// end group pcl
/// \}

}//	end of namespace

#endif
