#ifndef __H__UG_heightfield_util
#define __H__UG_heightfield_util

#include "lib_grid/lg_base.h"
#include "common/util/field.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
///	A heightfield represents a grid of number-values together with descriptors for
///	the cell dimensions and the total offset.
/** You can manually create heightfields or load them from a file. Given a heightfield
 * you may create a quad/tri grid or extract the heightfields boundary into a grid.
 * Have a look at field_util.h, too, which declares some useful functions (e.g. blurring).
 * \sa LoadHeightfieldFromASC, CreateGridFromField, CreateGridFromFieldBoundary
 */
UG_API class Heightfield{
	public:
		Heightfield();
		
	///	returns the interpolated value at the given location.
	/** returns the interpolated value at the given location. Through 'interpOrder'
	 * one may specify the order of interpolation:
	 *	- 0: piecewise constant (nearest entry)
	 *	- 1: linear
	 *
	 * \todo: add a method 'set_interpolation_method' to also support
	 * bilinear- or spline-interpolation
	 * \{ */

		number	interpolate (number x, number y, int interpOrder) const;
		
		number	interpolate (const vector2& c, int interpOrder) const
		{
			return interpolate(c.x(), c.y(), interpOrder);
		}

		number	interpolate (number x, number y) const
		{
			return interpolate(x, y, 0);
		}

		number	interpolate (const vector2& c) const
		{
			return interpolate(c.x(), c.y(), 0);
		}
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

	///	Smoothens the field by adjusting the value of each pixel towards the average of its neighbours
		void blur(number alpha, size_t numIterations);

	///	eliminates invalid cells by repeatedly filling those cells with averages of neighboring cells
	/** The field has to contain at least one valid cell. If it doesn't, false is returned.*/
		bool eliminate_invalid_cells();
	
	private:
		Field<number>	m_field;
		vector2			m_cellSize;
		vector2			m_offset;
		number			m_noDataValue;
};


////////////////////////////////////////////////////////////////////////////////
void UG_API
LoadHeightfieldFromASC(Heightfield& field, const char* filename);

////////////////////////////////////////////////////////////////////////////////
void UG_API
CreateGridFromField(Grid& grid,
					const Heightfield& hfield,
					Grid::VertexAttachmentAccessor<APosition> aaPos);

////////////////////////////////////////////////////////////////////////////////
void UG_API
CreateGridFromFieldBoundary(Grid& grid,
							const Heightfield& hfield,
							Grid::VertexAttachmentAccessor<APosition> aaPos);

}//	end of namespace

#endif	//__H__UG_heightfield_util
