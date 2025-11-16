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

#ifndef __H__UG__distribution__
#define __H__UG__distribution__

#include <vector>
#include "lib_grid/lg_base.h"
#include "lib_grid/algorithms/serialization.h"
#include "pcl/pcl_process_communicator.h"

namespace ug{


enum InterfaceStates{
	IS_UNASSIGNED = 0,
	IS_NORMAL = 1,
	IS_VMASTER = 1<<1,
	IS_VSLAVE = 1<<2,
	IS_DUMMY = 1<<3,
	HAS_PARENT = 1<<4 // only used for dummies currently
};


///	distributes/redistributes parts of possibly distributed grids.
/**	This method is still in development... Use with care!
 *
 * Note that the method assumes that all data serialized and deserialized
 * is stored in a consistent manner - i.e. that all slaves hold the same
 * values as their associated masters.
 *
 * The method only considers the partition defined on the elements of highest
 * dimension of the given grid. Elements of lower dimension are simply sent
 * alongside those highest dimensional elements.
 *
 * The method posts the following messages to the message hub of the specified grid:
 * 	- GridMessage_Distribution(GMDT_DISTRIBUTION_STARTS) at the start of the method
 * 	- GridMessage_Creation(GMCT_CREATION_STARTS) before the local grid is cleared
 * 	- GridMessage_Creation(GMCT_CREATION_STOPS) after the local grid has been rebuilt completely
 * 	- GridMessage_Distribution(GMDT_DISTRIBUTION_STOPS) when the method is done.
 *
 * \note	??? The partition map shPartition has to fulfill the following requirement:
 * 			All siblings (all children of a parent) have to be in a common subset
 * 			(not necessarily the subset of their parent). ???
 *
 * \param	procMap is by default nullptr and thus ignored. If you specify it
 * 			make sure to specify a pointer to an array with size
 * 			shPartition.num_subsets(). All values in the array have to be
 * 			in the range [0, pcl:NumProcs()[.
 * 			The procMap associates a process rank with each subset index.
 */
bool DistributeGrid(MultiGrid& mg,
					SubsetHandler& shPartition,
					GridDataSerializationHandler& serializer,
					bool createVerticalInterfaces,
					const std::vector<int>* processMap = nullptr,
					const pcl::ProcessCommunicator& procComm =
												pcl::ProcessCommunicator());

}// end of namespace

#endif
