/*
 * Copyright (c) 2016:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG_broadcast
#define __H__UG_broadcast

#include "lib_grid/tools/selector_grid.h"
#include "lib_grid/algorithms/serialization.h"
#include "pcl/pcl_process_communicator.h"

namespace ug{

///	Broadcasts the specified Selection from 'root' to all processes in procCom
/** Selected elements and associated attachments are broadcasted from
 * 'root' to all processes in procCom (including 'root').
 *
 * @param serializer	Has to operate on 'sel.grid()'. It is responsible to serialize
 *						e.g. attachments, subset-handlers or selectors.
 *
 * @param deserializer	Has to operate on 'gridOut'. It has to contain exactly the same
 *						components as 'serializer', except that all components have
 *						to operate on 'gridOut', too.

 * \note	In order to make sure that all required sides and vertices of selected
 *			elements are broadcasted, the given selector may be adjusted.*/
void BroadcastGrid(	Grid& gridOut,
					Selector& sel,
					GridDataSerializationHandler& serializer,
					GridDataSerializationHandler& deserializer,
					int root,
					const pcl::ProcessCommunicator& procCom =
												pcl::ProcessCommunicator());
 
}//	end of namespace

#endif	//__H__UG_broadcast
