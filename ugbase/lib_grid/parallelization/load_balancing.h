// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 21.04.2011 (m,d,y)

#ifndef __H__UG__load_balancing__
#define __H__UG__load_balancing__

#include <vector>
#include "lib_grid/lg_base.h"
#include "lib_grid/algorithms/geom_obj_util/geom_obj_util.h"
#include "util/parallel_callbacks.h"
#include "distributed_grid.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
//	PartitionElementsByRepeatedIntersection
///	partitions a grid into subsets that contain the same number of elements (as far as possible).
/**
 * assigns elements to subsets.
 * Tries to construct subsets that all have the same size.
 * Runs with something like O(n*log(n)).
 * It can not be guaranteed, that subsets are connected.
 */
template <class TElem, int IDimension, class TAPosition>
bool PartitionElementsByRepeatedIntersection(SubsetHandler& shOut,
										Grid& grid,
										int numSubsets,
										TAPosition& aVrtPos,
										int startDim = 0);

////////////////////////////////////////////////////////////////////////
//	PartitionElementsByRepeatedIntersection
///	partitions a grid into subsets that contain the same number of elements (as far as possible).
/**
 * assigns elements to subsets.
 * Tries to construct subsets that all have the same size.
 * Runs with something like O(n*log(n)).
 * It can not be guaranteed, that subsets are connected.
 */
template <class TElem, int IDimension, class TAPosition>
bool PartitionElementsByRepeatedIntersection(SubsetHandler& shOut,
										MultiGrid& mg,
										int level,
										int numSubsets,
										TAPosition& aVrtPos,
										int startDim = 0);

////////////////////////////////////////////////////////////////////////////////
///	Partitions the elements in the grid by sorting them into a regular grid.
/**	Let xInd, yInd be the indices of the cell in which an element lies.
 * the associated subset index is then calculated by
 * \code
 * 		subsetIndex = yInd * numCellsX + xInd;
 * \endcode
 *
 * \param bucketSubset	All elements which shall not be considered are assigned
 * 						to this subset (default -1).
 */
template <class TElem, class TIterator, class TAAPos>
bool PartitionElements_RegularGrid(SubsetHandler& shOut,
								TIterator begin, TIterator end,
								int numCellsX, int numCellsY,
								TAAPos& aaPos,
								boost::function<bool (TElem*)> cbConsiderElem
									= (bool(*)(TElem*))ConsiderAll,
								int bucketSubset = -1);

}//	end of namespace


////////////////////////////////
//	include implementation
#include "load_balancing_impl.hpp"

#endif
