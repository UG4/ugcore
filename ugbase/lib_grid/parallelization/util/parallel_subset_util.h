//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m07 d14

#ifndef __H__LIB_GRID__PARALLELL_SUBSET_UTIL__
#define __H__LIB_GRID__PARALLELL_SUBSET_UTIL__

#include "../distributed_grid.h"
#include "lib_grid/tools/subset_handler_multi_grid.h"
#include "pcl/pcl_process_communicator.h"

namespace ug
{

/// \addtogroup lib_grid_parallelization
/// @{

////////////////////////////////////////////////////////////////////////
template <class TElem>
void CollectSurfaceViewElements(ISubsetHandler& surfaceViewOut,
                                DistributedGridManager& distGridMgr,
								MultiGridSubsetHandler& mgsh,
								bool clearContainer);

////////////////////////////////////////////////////////////////////////
//	CreateSurfaceView
/**
 * Elements which are on the surface of the multi-grid-hierarchy
 * (elements that don't have children) are assigned to a subset of the
 * shSurfaceViewOut. The subset-index is taken from sh.
 *
 * Special care is taken for so called ghost-elements. Ghost elements
 * are elements that lie on the surface of the multi-grid but shall
 * not be considered for the suface-view (e.g. vertical-masters).
 *
 * If the distGridMgr is not operating on a multi-grid, then the method won't add
 * elements to the surface-view.
 */
template <class TSurfaceView>
void CreateSurfaceView(TSurfaceView& surfaceViewOut,
                       DistributedGridManager& distGridMgr,
						MultiGridSubsetHandler& mgsh);


////////////////////////////////////////////////////////////////////////
///	Gathers the global dimension of each subset
/**	Globally communicates the max-dimension of each subset and stores
 * it in the specified subset-property.
 * You may optionally specify a process-communicator, which shall be used for
 * communication.
 *
 * The dimension is set to -1, if the subset does not contain any elements at all.
 *
 * Note that all subset-handlers on all processes have to have the same number
 * of subsets!
 */
void UpdateGlobalMaxDimensionOfSubset(ISubsetHandler& sh,
						const std::string propertyName,
						pcl::ProcessCommunicator com = pcl::ProcessCommunicator());

/// @}
}//	end of namespace

////////////////////////////////
//	include implementation
#include "parallel_subset_util_impl.hpp"

#endif
