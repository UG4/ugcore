/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__LIB_GRID__parallelization_util_impl__
#define __H__UG__LIB_GRID__parallelization_util_impl__

#include "util/compol_copy_attachment.h"
#include "pcl/pcl_interface_communicator.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
template <typename TGeomObj>
void CreateAndDistributeGlobalIDs(Grid& g, GridLayoutMap& glm,
								  AGeomObjID& aID)
{
//	make sure that aID is attached to the objects of the grid
	if(!g.has_attachment<TGeomObj>(aID)){
		g.attach_to<TGeomObj>(aID);
	}

	Grid::AttachmentAccessor<TGeomObj, AGeomObjID> aaID(g, aID);

//	set up local ids
	int rank = pcl::ProcRank();

	using GeomObjIter = typename geometry_traits<TGeomObj>::iterator;

	size_t count = 0;
	for(GeomObjIter iter = g.begin<TGeomObj>();
		iter != g.end<TGeomObj>(); ++iter, ++count)
	{
		aaID[*iter] = MakeGeomObjID(rank, count);
	}

//	distribute the ids master->slave
	using Layout = typename GridLayoutMap::Types<TGeomObj>::Layout;
	pcl::InterfaceCommunicator<Layout> com;
	ComPol_CopyAttachment<Layout, AGeomObjID> compolCopy(g, aID);

	com.exchange_data(glm, InterfaceNodeTypes::INT_H_MASTER, InterfaceNodeTypes::INT_H_SLAVE, compolCopy);
	com.communicate();

//	we copy data from vertical slaves to vertical masters and not vice versa for
//	the following reason:
//	Multiple copies of low dimensional elements may reside in vertical master
//	interfaces, and they do not necessarily have the same global id yet.
//	Elements in horizontal interfaces, however, already have the same ids yet.
//	Since all vertical slaves are either unique or in a horizontal interface,
//	we can simply copy their global ids to vertical master elements.
	com.exchange_data(glm, InterfaceNodeTypes::INT_V_SLAVE, InterfaceNodeTypes::INT_V_MASTER, compolCopy);
	com.communicate();
}

}//	end of namespace

#endif
