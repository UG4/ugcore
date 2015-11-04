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
 */
template <typename AType>
void AttachmentAllReduce
(
	Grid& grid,
	AType aValue,
	pcl::ReduceOperation op
)
{
	typedef typename GridLayoutMap::Types<Vertex>::Layout layout_t;
	
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
