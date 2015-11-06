/*
 * Copyright (c) 2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG_field_util
#define __H__UG_field_util

#include "lib_grid/lg_base.h"
#include "common/util/field.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
void UG_API
CreateGridFromField(Grid& grid,
					const Field<number>& field,
					const vector2& cellSize,
					const vector2& offset,
					number noDataValue,
					Grid::VertexAttachmentAccessor<APosition> aaPos);


////////////////////////////////////////////////////////////////////////////////
void UG_API
CreateGridFromFieldBoundary(Grid& grid,
					const Field<number>& field,
					const vector2& cellSize,
					const vector2& offset,
					number noDataValue,
					Grid::VertexAttachmentAccessor<APosition> aaPos);


////////////////////////////////////////////////////////////////////////////////
///	Smoothens the field by adjusting the value of each pixel towards the average of its neighbours
/** The value type T has to support operators += and *= and = 0*/
template <class T>
void BlurField(Field<T>& field, number alpha, size_t numIterations, const T& noDataValue);


////////////////////////////////////////////////////////////////////////////////
///	eliminates invalid cells by repeatedly filling those cells with averages of neighboring cells
/** The field has to contain at least one valid cell. If it doesn't, false is returned.
 * The value type T has to support operators += and *= and = 0*/
template <class T>
bool EliminateInvalidCells(Field<T>& field, const T& noDataValue);


////////////////////////////////////////////////////////////////////////////////
///	invalidates cells that belong to a small lense
/**	Given a valid cell, the assoicated lense is the set of valid cells that can
 * be reached from that cell by only traversing valid neighbors.
 * Whether a lense is small or not is determined by the number of cells that belong
 * to a lense.*/
template <class T>
void InvalidateSmallLenses(Field<T>& field, size_t thresholdCellCount,
						   const T& noDataValue);

}//	end of namespace

////////////////////////////////////////
//	include implementation
#include "field_util_impl.h"

#endif	//__H__UG_field_util
