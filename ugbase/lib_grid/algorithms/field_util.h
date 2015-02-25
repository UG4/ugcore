// created by Sebastian Reiter
// s.b.reiter@gmail.com

#ifndef __H__UG_field_util
#define __H__UG_field_util

#include "lib_grid/lg_base.h"
#include "common/util/field.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
///	A heightfield represents a grid of number-values together with descriptors for
///	the cell dimensions and the total offset.
/** You can manually create heightfields or load them from a file. Given a heightfield
 * you may create a quad/tri grid or extract the heightfields boundary into a grid.
 * \sa LoadHeightfieldFromASC, CreateGridFromField, CreateGridFromFieldBoundary
 */
UG_API class Heightfield{
	public:
		Heightfield();
		
	///	returns the value of the closest entry in the m_data field.
	/** \todo: add a method 'set_interpolation_method' to also support
	 * bilinear- or spline-interpolation
	 * \{ */
		number	interpolate(number x, number y) const;
		number	interpolate(const vector2& c) const		{return interpolate(c.x(), c.y());}
	/** \} */

	///	returns the index-tuple of the closest field-entry
		std::pair<int, int> coordinate_to_index(number x, number y) const;

	///	returns the coordinate of the given cell (specified through an index-tuple)
		vector2 index_to_coordinate(int ix, int iy) const;

		Field<number>	field;
		vector2			cellSize;
		vector2			offset;
		number			noDataValue;
};



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
CreateGridFromField(Grid& grid,
					const Heightfield& hfield,
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
void UG_API
CreateGridFromFieldBoundary(Grid& grid,
					const Heightfield& hfield,
					Grid::VertexAttachmentAccessor<APosition> aaPos);

}//	end of namespace

#endif	//__H__UG_field_util
