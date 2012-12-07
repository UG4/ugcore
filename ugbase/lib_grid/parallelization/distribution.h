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
 * \note	??? The partition map shPartition has to fulfill the following requirement:
 * 			All siblings (all children of a parent) have to be in a common subset
 * 			(not necessarily the subset of their parent). ???
 *
 * \param	procMap is by default NULL and thus ignored. If you specify it
 * 			make sure to specify a pointer to an array with size
 * 			shPartition.num_subsets(). All values in the array have to be
 * 			in the range [0, pcl:GetNumProcesses()[.
 * 			The procMap associates a process rank with each subset index.
 *
 * \todo	only use one GridDataSerializationHandler.
 */
bool DistributeGrid(MultiGrid& mg,
					SubsetHandler& shPartition,
					GridDataSerializationHandler& serializer,
					GridDataSerializationHandler& deserializer,
					bool createVerticalInterfaces,
					std::vector<int>* processMap = NULL,
					const pcl::ProcessCommunicator& procComm =
												pcl::ProcessCommunicator());

}// end of namespace

#endif
