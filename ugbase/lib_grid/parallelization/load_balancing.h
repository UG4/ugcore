// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 21.04.2011 (m,d,y)

#ifndef __H__UG__load_balancing__
#define __H__UG__load_balancing__

#include <vector>
#include "lib_grid/algorithms/geom_obj_util/geom_obj_util.h"
#include "lib_grid/tools/subset_handler_grid.h"
#include "lib_grid/multi_grid.h"
#include "util/parallel_callbacks.h"
#include "lib_grid/callbacks/basic_callbacks.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////////////
///	Partitions the elements in the grid by sorting them into a regular grid.
/**	This method uses Grid::mark.
 * Let xInd, yInd be the indices of the cell in which an element lies.
 * the associated subset index is then calculated by
 * \code
 * 		subsetIndex = zInd * numCellsX * numCellsY + yInd * numCellsX + xInd;
 * \endcode
 *
 * \param bucketSubset	All elements which shall not be considered are assigned
 * 						to this subset (default -1).
 */
template <class TElem, class TIterator, class TAAPos>
bool PartitionElements_RegularGrid(SubsetHandler& shOut,
								TIterator begin, TIterator end,
								int numCellsX, int numCellsY, int numCellsZ,
								TAAPos& aaPos,
								typename Grid::traits<TElem>::callback cbConsiderElem
									= ConsiderAll(),
								int bucketSubset = -1);

////////////////////////////////////////////////////////////////////////////////
///	Partitions the elements in the grid using the METIS library
/**	This method calls METIS_PartGraphKway. Note that METIS is an external library
 * developed at Karypis Labs (http://glaros.dtc.umn.edu/gkhome/)
 *
 * Note that this method is best suited for partitions with more than 8 procs.
 * For less than 8 procs metis features other, better suited methods.
 *
 * Valid template arguments are Edge, Face, Volume and derived types.
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
 * Valid template arguments are Edge, Face, Volume and derived types.
 */
template <class TGeomBaseObj>
bool PartitionMultiGrid_MetisKway(SubsetHandler& shPartitionOut,
							 	  MultiGrid& grid, int numParts,
							 	  size_t baseLevel = 0,
							 	  int hWeight = 1, int vWeight = 1);

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
 * weightFct specifies a function that attributes special weights to edges on the
 * dual graph.
 *
 * Valid template arguments are Edge, Face, Volume and derived types.
 */
template <class TGeomBaseObj>
bool PartitionMultiGrid_MetisKway(SubsetHandler& shPartitionOut,
							 	  MultiGrid& grid, int numParts, size_t baseLevel,
							 	  boost::function<int (TGeomBaseObj*, TGeomBaseObj*)>& weightFct);

////////////////////////////////////////////////////////////////////////////////
///	Partitions the elements in the multi-grid using the METIS library
/**	This method calls METIS_PartGraphKway. Note that METIS is an external library
 * developed at Karypis Labs (http://glaros.dtc.umn.edu/gkhome/)
 *
 * Note that this method is best suited for partitions with more than 8 procs.
 * For less than 8 procs metis features other, better suited methods.
 *
 * The method performs load balancing for the elements in the given level. The
 * elements are weighted according to the number of children each has.
 * Child elements will then be recursively assigned to the partitions into which
 * their parents have been assigned, starting from level+1.
 * Elements below the specified level will be assigned to the local process id.
 */
template <class TGeomBaseObj>
bool PartitionMultiGridLevel_MetisKway(SubsetHandler& shPartitionOut,
							 	  MultiGrid& mg, int numParts, size_t level);

////////////////////////////////////////////////////////////////////////////////
///	Partitions the elements in the multi-grid using the PARMETIS library
/**	This method calls METIS_PartGraphKway. Note that PARMETIS is an external library
 * developed at Karypis Labs (http://glaros.dtc.umn.edu/gkhome/)
 *
 * The method performs parallel load balancing for the elements in the given level. The
 * elements are weighted according to the number of children each has.
 * Child elements will then be recursively assigned to the partitions into which
 * their parents have been assigned, starting from level+1.
 * Elements below the specified level will be assigned to the local process id.
 */
template <class TGeomBaseObj>
bool PartitionMultiGridLevel_ParmetisKway(SubsetHandler& shPartitionOut,
							 	  	  MultiGrid& mg, int numParts, size_t level);

}//	end of namespace


////////////////////////////////
//	include implementation
#include "load_balancing_impl.hpp"

#endif
