// created by Sebastian Reiter
// s.b.reiter@gmail.com
// 29.11.2012 (d,m,y)

#ifndef __H__UG__distribution__
#define __H__UG__distribution__

#include <vector>
#include "lib_grid/lg_base.h"
#include "lib_grid/algorithms/serialization.h"
#include "pcl/pcl_process_communicator.h"

namespace ug{

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
 * \param	procMap is by default NULL and thus ignored. If you specify it
 * 			make sure to specify a pointer to an array with size
 * 			shPartition.num_subsets(). All values in the array have to be
 * 			in the range [0, pcl:NumProcs()[.
 * 			The procMap associates a process rank with each subset index.
 */
bool DistributeGrid(MultiGrid& mg,
					SubsetHandler& shPartition,
					GridDataSerializationHandler& serializer,
					bool createVerticalInterfaces,
					const std::vector<int>* processMap = NULL,
					const pcl::ProcessCommunicator& procComm =
												pcl::ProcessCommunicator());

}// end of namespace

#endif
