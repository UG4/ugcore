/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

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
