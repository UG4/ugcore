/*
 * ass_adapter.h
 *
 *  Created on: 04.02.2013
 *      Author: raphaelprohl
 */

#ifndef ASS_ADAPTER_H_
#define ASS_ADAPTER_H_

#include "lib_grid/tools/bool_marker.h"
#include "lib_grid/tools/selector_grid.h"

namespace ug{

struct AssIndex{

	///	should assemble be index-wise
	bool index_set;

	///	current index
	size_t index;
};

// AssAdapter combines tools to adapt the assemble routine
struct AssAdapter{

	///	marker used to skip elements
	BoolMarker* pBoolMarker;

	///	selector used to set a list of elements for the assembling
	Selector* 	pSelector;

	///	object for index-wise assemble routine
	AssIndex assIndex;
};

} // end namespace ug

#endif /* ASS_ADAPTER_H_ */
