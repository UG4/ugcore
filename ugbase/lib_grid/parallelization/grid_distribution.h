//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m12 d08

#ifndef __H__LIB_GRID__GRID_DISTRIBUTION__
#define __H__LIB_GRID__GRID_DISTRIBUTION__

//	will be removed!
#include <vector>
#include "parallel_grid_layout.h"
#include "distributed_grid.h"
#include "lib_grid/algorithms/serialization.h"
#include "pcl/pcl_base.h"
#include "pcl/pcl_process_communicator.h"

namespace ug
{

/// \addtogroup lib_grid_parallelization_distribution
///	@{

////////////////////////////////////////////////////////////////////////
//	DistributeGrid
///	distributes a grid to several processes.
/**
 * Partitions the grid into parts specified by the SubsetHandler sh.
 * For partitioning the grid the algorithm CreateGridLayouts
 * is used internally. Have a look at its documentation to see how you
 * have to specify your parts in the SubsetHandler.
 * The part for the local (the calling) process (specified by localProcID)
 * won't be send through the network, instead it will be directly written
 * to pLocalGridOut, pLocalSubsetHandlerOut and pLocalGridCommSetOut - if
 * those are specified.
 * All three should be empty before calling this method. Passing NULL only makes
 * sense if you specified a processMap that has no entry for the calling
 * process.
 * You may optionally specify a process-map. This map is represented by
 * a std::vector<int> which should have as many entries as there are
 * subsets in the SubsetHandler. Each entry specifies the target process
 * for the associated subset. You may pass NULL as parameter. The grids
 * are distributed to processes 0 to sh.num_subsets()-1 in this case.
 *
 * Grids distributed through this method may be received by ReceiveGrid.
 */
bool DistributeGrid(MultiGrid& mg, ISubsetHandler& sh,
					SubsetHandler& shPartition, int localProcID,
					MultiGrid* pLocalGridOut = NULL,
					ISubsetHandler* pLocalSHOut = NULL,
					GridLayoutMap* pLocalGridLayoutMapOut = NULL,
					std::vector<int>* pProcessMap = NULL);

////////////////////////////////////////////////////////////////////////
///	copies parts of the grid to other processes.
/**
 * Copies the parts of the grid to the processes given in the subset-handler
 * (optionally to the processes in the process-map).
 * In contrary to DistributeGrid(...) the given multigrid is regarded
 * as the local grid-instance. The interfaces in the layoutMap will thus
 * be created for the given grid instance.
 * Vertical interfaces will be created that allow to communicate data
 * between the local VerticalMasters and their copies (vertical slaves).
 * All elements that are copied to other processes will be regarded as
 * vertical slaves there.
 * To make sure that those vertical interfaces are created on the
 * receiving processes, you have to call ReceiveGrid(...) with
 * createVerticalLayouts == true.
 */
bool DistributeGrid_KeepSrcGrid(MultiGrid& mg, ISubsetHandler& sh,
								GridLayoutMap& layoutMap,
								SubsetHandler& shPartition,
								int localProcID,
								bool distributeGenealogy = false,
								std::vector<int>* pProcessMap = NULL);
								
////////////////////////////////////////////////////////////////////////
//	ReceiveGrid
///	receives a part of a grid that was distributed through DistributeGrid.
/**
 * gridOut and gridLayoutOut should be empty when passed to this method.
 * srcProcID has to specify the process that distributes the grids.
 * The keys in the gridLayoutMap correspond to the type of the nodes that
 * each layout contains.
 */
bool ReceiveGrid(MultiGrid& mgOut, ISubsetHandler& shOut,
				GridLayoutMap& gridLayoutMapOut,
				int srcProcID,
				bool createVerticalLayouts = false);
				
				
///	redistributes parts of distributed grids.
/**	This method is still in development... Use with care!
 *
 * Note that the method assumes that all data serialized and deserialized
 * is stored in a consistent manner - i.e. that all slaves hold the same
 * values as their associated masters.
 *
 * \param	procMap is by default NULL and thus ignored. If you specify it
 * 			make sure to specify a pointer to an array with size
 * 			shPartition.num_subsets(). All values in the array have to be
 * 			in the range [0, pcl:GetNumProcesses()[.
 * 			The procMap associates a process rank with each subset index.
 */
bool RedistributeGrid(DistributedGridManager& distGridMgrInOut,
					  SubsetHandler& shPartition,
					  const GridDataSerializationHandler& serializer,
					  const GridDataSerializationHandler& deserializer,
					  bool createVerticalInterfaces,
					  std::vector<int>* processMap = NULL,
					  const pcl::ProcessCommunicator& procComm =
							  	  	  pcl::ProcessCommunicator());
					  
/// @}
}//	end of namespace

#endif
