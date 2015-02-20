// created by Sebastian Reiter
// s.b.reiter@gmail.com

#ifndef __H__UG_field_util
#define __H__UG_field_util

#include "lib_grid/lg_base.h"
#include "common/util/field.h"

namespace ug{

void UG_API
CreateGridFromField(Grid& grid,
					const Field<number>& field,
					const vector2& cellSize,
					const vector2& offset,
					number noDataValue,
					Grid::VertexAttachmentAccessor<APosition> aaPos);

}//	end of namespace

#endif	//__H__UG_field_util
