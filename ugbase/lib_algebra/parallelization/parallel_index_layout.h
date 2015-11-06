/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Authors: Andreas Vogel, Sebastian Reiter
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

#ifndef __H__LIB_ALGEBRA__PARALLELIZATION__PARALLEL_INDEX_LAYOUT__
#define __H__LIB_ALGEBRA__PARALLELIZATION__PARALLEL_INDEX_LAYOUT__

#include <set>
#include <vector>
#include <iostream>
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

/**
 * Marks all indices in the IndexLayout::Interface with true
 * @param mark a std::vector which has to be big enough for maximum index in interface
 * @param layout
 */
void MarkAllFromInterface(std::vector<bool> &mark, const IndexLayout::Interface &interface);

template <class T>
void MarkAllFromInterface(std::vector<T> &mark, const IndexLayout::Interface &interface, const T &default_val)
{
	for(IndexLayout::Interface::const_iterator iter = interface.begin(); iter != interface.end(); ++iter)
		mark[ interface.get_element(iter) ] = default_val;
}

/**
 * Marks all indices in the IndexLayout with true
 * @param mark a std::vector which has to be big enough for maximum index in layout
 * @param layout
 */
void MarkAllFromLayout(std::vector<bool> &mark, const IndexLayout &layout);

template <class T>
void MarkAllFromLayout(std::vector<T> &mark, const IndexLayout &layout, const T &default_val)
{
	for(IndexLayout::const_iterator iter = layout.begin(); iter != layout.end(); ++iter)
		MarkAllFromInterface<T> (mark, layout.interface(iter), default_val);
}

/**
  * @param mark a set to which all indices from the interface are added
 * @param layout
 */
void AddAllFromInterface(std::set<size_t> &s, const IndexLayout::Interface &interface);


/**
 * @param mark a set to which all indices from the interface are added
 * @param layout
 */
void AddAllFromLayout(std::set<size_t> &s, const IndexLayout &layout);


std::ostream &operator << (std::ostream &out, const IndexLayout &layout);

} // end namespace ug

#endif /* __H__LIB_ALGEBRA__PARALLELIZATION__PARALLEL_INDEX_LAYOUT__ */
