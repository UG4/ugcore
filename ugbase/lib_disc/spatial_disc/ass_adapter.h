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
	bool index_set;
	size_t index; //TODO: ist der index Ÿberhaupt notwendig???
};

struct AssAdapter{

	///	marker used to skip elements
	BoolMarker* pBoolMarker;

	///	selector used to set a list of elements for the assembling
	Selector* 	pSelector;

	///	one specific index for which the assemble quantities
	/// (defect, jacobian, ...) are build up
	AssIndex assIndex;
};

} // end namespace ug

#endif /* ASS_ADAPTER_H_ */
