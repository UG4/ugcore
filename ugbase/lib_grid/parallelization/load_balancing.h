// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 21.04.2011 (m,d,y)

#ifndef __H__UG__load_balancing__
#define __H__UG__load_balancing__

#include <vector>
#include "lib_grid/lg_base.h"
#include "lib_grid/algorithms/geom_obj_util/geom_obj_util.h"
#include "util/parallel_callbacks.h"


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
/**	This method uses Grid::mark.
 * Let xInd, yInd be the indices of the cell in which an element lies.
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
								typename Grid::traits<TElem>::callback cbConsiderElem
									= Grid::traits<TElem>::cb_consider_all,
								int bucketSubset = -1);

////////////////////////////////////////////////////////////////////////////////
///	Partitions the elements in the grid using the METIS library
/**	This method calls METIS_PartGraphKway. Note that METIS is an external library
 * developed at Karypis Labs (http://glaros.dtc.umn.edu/gkhome/)
 *
 * Note that this method is best suited for partitions with more than 8 procs.
 * For less than 8 procs metis features other, better suited methods.
 *
 * Valid template arguments are EdgeBase, Face, Volume and derived types.
 */
template <class TGeomBaseObj>
bool PartitionGrid_MetisKway(SubsetHandler& shPartitionOut,
							 Grid& grid, int numParts);

////////////////////////////////////////////////////////////////////////////////
///	Partitions the elements in the multi-grid using the METIS library
/**	This method calls METIS_PartGraphKway. Note that METIS is an external library
 * developed at Karypis Labs (http://glaros.dtc.umn.edu/gkhome/)
 *
 * Note that this method is best suited for partitions with more than 8 procs.
 * For less than 8 procs metis features other, better suited methods.
 *
 * All elements in baseLevel and higher levels will be partitioned. elements
 * below baseLevel will stay where they are and are completely ignored during
 * load balancing.
 *
 * hWeight and vWeight determine, how important it is to keep horizontal
 * and vertical neighbors on the same process as the element itself.
 * The bigger hWeight, the more attention is spend to keep neighbors together.
 * Both parameters have to be > 0. Default is 1.
 *
 * Valid template arguments are EdgeBase, Face, Volume and derived types.
 */
template <class TGeomBaseObj>
bool PartitionMultiGrid_MetisKway(SubsetHandler& shPartitionOut,
							 	  MultiGrid& grid, int numParts,
							 	  size_t baseLevel = 0,
							 	  int hWeight = 1, int vWeight = 1);


}//	end of namespace


////////////////////////////////
//	include implementation
#include "load_balancing_impl.hpp"

#endif
