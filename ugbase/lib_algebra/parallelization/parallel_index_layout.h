/*
 * parallel_index_layout.h
 *
 *  Created on: 21.5.2010
 *      Author: A. Vogel, S.Reiter
 */

#ifndef __H__LIB_ALGEBRA__PARALLELIZATION__PARALLEL_INDEX_LAYOUT__
#define __H__LIB_ALGEBRA__PARALLELIZATION__PARALLEL_INDEX_LAYOUT__

#include <vector>
#include "pcl/pcl.h"

namespace ug
{

///	Allows communication between distributed vectors and matrices.
/**	Note that indices are stored in an std::vector in the moment.
 *	This allows fast iteration and memory allocation, if dynamic
 *	interfaces are required this may however be slower than a
 *	std::list container.
 */
typedef pcl::SingleLevelLayout<pcl::OrderedInterface<size_t, std::vector> >
		IndexLayout;


} // end namespace ug

#endif /* __H__LIB_ALGEBRA__PARALLELIZATION__PARALLEL_INDEX_LAYOUT__ */
