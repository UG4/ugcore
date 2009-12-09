//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m12 d08

#ifndef __H__LIB_GRID__GRID_DISTRIBUTION__
#define __H__LIB_GRID__GRID_DISTRIBUTION__

//	will be removed!
#include <vector>
#include "lib_grid/lg_base.h"
#include "parallel_grid_layout.h"

namespace ug
{

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
 * to pLocalGridOut and pLocalGridCommSetOut - if those are specified.
 * Both should be empty before calling this method. Passing NULL only makes
 * sense if you specified a processMap that has no entry for the calling
 * process.
 * You may optionally specify a process-map. This map is represented by
 * a std::vector<int> which should have as many entries as there are
 * subsets in the SubsetHandler. Each entry specifies the target process
 * for the associated subset. You may pass NULL as parameter. The grids
 * are distributed to processes 0 to sh.num_subsets()-1 in this case.
 *
 * Grids distributed through this method may be received by ReveiveGrid.
 */
void DistributeGrid(MultiGrid& mg, SubsetHandler& sh, int localProcID,
					MultiGrid* pLocalGridOut = NULL,
					ParallelGridLayout* pLocalGridLayoutOut = NULL,
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
void ReceiveGrid(MultiGrid& mgOut, ParallelGridLayout& gridLayoutOut,
					int srcProcID);
}//	end of namespace

#endif
