/*
 * Copyright (c) 2015:  G-CSC, Goethe University Frankfurt
 * Authors: Sebastian Reiter, Dmitriy Logashenko
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

#ifndef __H__UG_PCL_ATTACHMENT_OPERATIONS__
#define __H__UG_PCL_ATTACHMENT_OPERATIONS__

#include "compol_attachment_reduce.h"
#include "compol_copy_attachment.h"
#include "pcl/pcl_interface_communicator.h"

namespace ug{

/**
 * Syncronizes a given attachment using a given reduction operation.
 * Note that the syncronization takes place only for the horizontal
 * masters and slaves.
 * Use e.g. like this:
 *
 * \code
 * AttachmentAllReduce<Vertex>(grid, aValue, PCL_RO_ADD)
 * \endcode
 *
 * see pcl_methods.h for more PCL_RO_... reduce operations.
 */
template <typename TElem, typename AType>
void AttachmentAllReduce
(
	Grid& grid,
	AType aValue,
	pcl::ReduceOperation op
)
{
	using layout_t = typename GridLayoutMap::Types<TElem>::Layout;
	
	GridLayoutMap& glm = grid.distributed_grid_manager()->grid_layout_map();
	pcl::InterfaceCommunicator<layout_t> icom;

	ComPol_AttachmentReduce<layout_t, AType> cpAValue (grid, aValue, op);
	
	icom.exchange_data (glm, INT_H_SLAVE, INT_H_MASTER, cpAValue);

	icom.communicate();

	ComPol_CopyAttachment<layout_t, AType> cpCopyAValue (grid, aValue);
	
	icom.exchange_data (glm, INT_H_MASTER, INT_H_SLAVE, cpCopyAValue);
	
	icom.communicate();
}

} // end of namespace ug

#endif // __H__UG_PCL_ATTACHMENT_OPERATIONS__

/* End of File */
