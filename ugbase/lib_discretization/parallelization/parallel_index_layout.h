/*
 * parallel_index_layout.h
 *
 *  Created on: 21.5.2010
 *      Author: A. Vogel, S.Reiter
 */

#ifndef __H__LIB_DISCRETIZATION__PARALLELIZATION__PARALLEL_INDEX_LAYOUT__
#define __H__LIB_DISCRETIZATION__PARALLELIZATION__PARALLEL_INDEX_LAYOUT__

#include <vector>
#include "lib_grid/algorithms/parallelization/parallelization.h"

namespace ug
{

///	Allows communication between distributed vectors and matrices.
/**	Note that indices are stored in an std::vector in the moment.
 *	This allows fast iteration and memory allocation, if dynamic
 *	interfaces are required this may however be slower than a
 *	std::list container.
 */
/*
typedef pcl::MultiLevelLayout<
			pcl::OrderedInterface<uint, std::vector> > IndexLayout;
*/

typedef pcl::Layout<pcl::OrderedInterface<uint, std::vector> > IndexLayout;

/*
///	Holds master and slave layouts for each process
typedef pcl::LayoutMap<pcl::MultiLevelLayout,
						pcl::OrderedInterface,
						int,
						std::vector>	IndexLayoutMap;
*/

}//	end of namespace

#endif
