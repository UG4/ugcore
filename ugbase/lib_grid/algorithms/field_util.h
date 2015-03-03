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

	///	returns the x- and y-extent of the heightfield
		vector2 extent() const;

		Field<number>& field()					{return m_field;}
		const Field<number>& field() const		{return m_field;}

		const vector2& cell_size() const		{return m_cellSize;}
		void set_cell_size(const vector2& s)	{m_cellSize = s;}

		const vector2& offset() const			{return m_offset;}
		void set_offset(const vector2& o)		{m_offset = o;}

		number no_data_value() const			{return m_noDataValue;}
		void set_no_data_value(number val)		{m_noDataValue = val;}

	//	moves the heightfield by altering it's offset
		void move(const vector2& v)				{m_offset += v;}

	private:
		Field<number>	m_field;
		vector2			m_cellSize;
		vector2			m_offset;
		number			m_noDataValue;
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

////////////////////////////////////////////////////////////////////////////////
///	Smoothens the field by adjusting the value of each pixel towards the average of its neighbours
/** The value type T has to support operators += and *= and = 0*/
template <class T>
void BlurField(Field<T>& field, number alpha, size_t numIterations, const T& noDataValue);

}//	end of namespace

////////////////////////////////////////
//	include implementation
#include "field_util_impl.h"

#endif	//__H__UG_field_util
