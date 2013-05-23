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

///\ingroup lib_algebra_parallelization

///	Allows communication between distributed vectors and matrices.
/**	Note that indices are stored in an std::vector in the moment.
 *	This allows fast iteration and memory allocation, if dynamic
 *	interfaces are required this may however be slower than a
 *	std::list container.
 */
typedef pcl::SingleLevelLayout<pcl::OrderedInterface<size_t, std::vector> >
		IndexLayout;

///	Logs the internals of an index layout.
/**
 * Writes information about an index interface. If depth >= 1 is passed, then
 * also the current indices in the interfaces are printed.
 */
void LogIndexLayout(IndexLayout& layout, int depth = 0);


/// logs index infos for all procs successively
void LogIndexLayoutOnAllProcs(IndexLayout& layout, int depth = 0);

/// replaces the indices in the layout based on a passed mapping
/**
 * This function replaces the indices in the layout by new indices. The mapping
 * between old and new indices is passed as newIndex = vMap[oldIndex]. The size
 * of the new index set must be smaller or equal to the old index set. If the
 * entry in the map is negative, the old index is removed from the layout.
 */
void ReplaceIndicesInLayout(IndexLayout& layout, const std::vector<int>& vMap);



inline IndexLayout::Interface::iterator find(IndexLayout::Interface &interface, size_t i)
{
	for(IndexLayout::Interface::iterator iter = interface.begin(); iter != interface.end(); ++iter)
	{
		if(interface.get_element(iter) == i)
			return iter;
	}
	return interface.end();
}

inline bool IsInInterface(IndexLayout::Interface &interface, size_t i)
{
	return find(interface, i) != interface.end();
}

inline void AddIfUnique(IndexLayout::Interface &interface, size_t i)
{
	if(IsInInterface(interface, i) == false)
		interface.push_back(i);
}



} // end namespace ug

#endif /* __H__LIB_ALGEBRA__PARALLELIZATION__PARALLEL_INDEX_LAYOUT__ */
