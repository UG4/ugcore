//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m06 d29

#ifndef __H__LIB_GRID__PARALLELIZATION_UTIL__
#define __H__LIB_GRID__PARALLELIZATION_UTIL__

#include "distributed_grid.h"
#include <boost/function.hpp>

namespace ug
{
/// \addtogroup lib_grid_parallelization
///	@{
////////////////////////////////////////////////////////////////////////
//	utility methods
///	Returns the type of associated interfaces
/**
 * \param interfaceType: a constant enumerated in ug::InterfaceNodeTypes.
 * \return the associated interface type or INT_NONE, if no associated
 *		type is known.
 */
int GetAssociatedInterfaceType(int interfaceType);


////////////////////////////////////////////////////////////////////////
///	Creates and distributes global ids for the given element type.
/**	IDs are written to the given attachment (aGeomObjID by default).
 */
template <class TGeomObj>
void CreateAndDistributeGlobalIDs(Grid& g, GridLayoutMap& glm,
								  AGeomObjID& aID = aGeomObjID);

////////////////////////////////////////////////////////////////////////
///	Checks whether the grid-layout-map on this proc is consistent with connected ones.
bool TestGridLayoutMap(MultiGrid& mg, GridLayoutMap& glm);


////////////////////////////////////////////////////////////////////////
//	default implementations of grid-partitioning methods
///	partitions the grid by repeated bisection
/**	Volume grids are partitioned in 3 dimensions, face grids are
 *	partitioned in 2 dimensions and edge grids are partitioned in
 *	one dimension.
 */
bool PartitionGrid_Bisection(SubsetHandler& partitionOut,
							  MultiGrid& mg, ISubsetHandler& sh,
							  size_t numProcs);

///	@}
}//	end of namespace


////////////////////////////////
//	include implementation
#include "parallelization_util_impl.hpp"

#endif
