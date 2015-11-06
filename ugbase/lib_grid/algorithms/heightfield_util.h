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
