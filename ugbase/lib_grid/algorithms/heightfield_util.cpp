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

#include <cmath>
#include "heightfield_util.h"
#include "field_util.h"
#include "lib_grid/file_io/file_io_asc.h"
#include "lib_grid/file_io/file_io.h"
using namespace std;

namespace ug{
	
////////////////////////////////////////////////////////////////////////////////
Heightfield::Heightfield() :
	m_cellSize(1, 1),
	m_offset(0, 0),
	m_noDataValue(1e30)
{
}

number Heightfield::
interpolate(number x, number y, int interpOrder) const
{
	switch(interpOrder){
		case 0: {
			pair<int, int> ind = coordinate_to_index(x, y);
			if( ind.first >= 0 && ind.first < (int)m_field.width() &&
				ind.second >= 0 && ind.second < (int)m_field.height())
			{
				return m_field.at(ind.first, ind.second);
			}
		}break;

		case 1:{
			int cx = int((x - m_offset.x()) / m_cellSize.x());
			int cy = int((y - m_offset.y()) / m_cellSize.y());
			if(		(cx >= 0)
				&& (cy >= 0)
				&& (cx + 1 < (int)m_field.width())
				&& (cy + 1 < (int)m_field.height()) )
			{
				float vx = (x - ((number)cx * m_cellSize.x() + m_offset.x())) / m_cellSize.x();
				float vy = (y - ((number)cy * m_cellSize.y() + m_offset.y())) / m_cellSize.y();
				float ux = 1-vx;
				float uy = 1-vy;

				return	ux*uy*m_field.at(cx, cy) + vx*uy*m_field.at(cx+1,cy)
					  + ux*vy*m_field.at(cx, cy+1) + vx*vy*m_field.at(cx+1, cy+1);
			}
		}break;

		default:
			UG_THROW("Currently only interpolation orders 0 and 1 supported in Heightfield::interpolate")
	}

	return m_noDataValue;
}

std::pair<int, int> Heightfield::
coordinate_to_index(number x, number y) const
{
//	roundOffset 0: value is constant across each cell
//	roundOffset 0.5: value is constant around each node
	const number roundOffset = 0.5;

	pair<int, int> c;
	if(m_cellSize.x() != 0)
		c.first = (int)(roundOffset + (x - m_offset.x()) / m_cellSize.x());
	else
		c.first = 0;
	
	if(m_cellSize.y() != 0)
		c.second = (int)(roundOffset + (y - m_offset.y()) / m_cellSize.y());
	else
		c.second = 0;
	return c;
}

vector2 Heightfield::index_to_coordinate(int ix, int iy) const
{
	return vector2(	m_offset.x() + (number)ix * m_cellSize.x(),
					m_offset.y() + (number)iy * m_cellSize.y());
}

vector2 Heightfield::extent() const
{
	return vector2(m_cellSize.x() * (number)m_field.width(),
				   m_cellSize.y() * (number)m_field.height());
}

void Heightfield::blur(number alpha, size_t numIterations)
{
	BlurField(field(), alpha, numIterations, no_data_value());
}


bool Heightfield::eliminate_invalid_cells()
{
	return EliminateInvalidCells(field(), no_data_value());
}




////////////////////////////////////////////////////////////////////////////////
void CreateGridFromField(Grid& grid,
						 const Heightfield& hfield,
						 Grid::VertexAttachmentAccessor<APosition> aaPos)
{
	CreateGridFromField(grid, hfield.field(), hfield.cell_size(), hfield.offset(),
						hfield.no_data_value(), aaPos);
}

////////////////////////////////////////////////////////////////////////////////
void CreateGridFromFieldBoundary(Grid& grid,
								 const Heightfield& hfield,
								 Grid::VertexAttachmentAccessor<APosition> aaPos)
{
	CreateGridFromFieldBoundary(grid, hfield.field(), hfield.cell_size(), hfield.offset(),
								hfield.no_data_value(), aaPos);
}

////////////////////////////////////////////////////////////////////////////////
void LoadHeightfieldFromASC(Heightfield& hfield, const char* filename)
{
	std::string name;
	bool fileFound = FindFileInStandardGridPaths(name, filename);
	UG_COND_THROW(!fileFound, "Couldn't locate file " << filename);

	FileReaderASC reader;
	reader.set_field(&hfield.field());
	reader.load_file(name.c_str());
	hfield.set_cell_size(vector2(reader.cell_size(), reader.cell_size()));
	hfield.set_offset(vector2(reader.lower_left_corner_x(),
							  reader.lower_left_corner_y()));
	hfield.set_no_data_value(reader.no_data_value());
}

}//	end of namespace
