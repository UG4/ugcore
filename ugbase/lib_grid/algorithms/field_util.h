// created by Sebastian Reiter
// s.b.reiter@gmail.com

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
