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
#include "common/util/file_util.h"
#include "lib_grid/file_io/file_io_asc.h"

using namespace std;

namespace ug {
	
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
			if( ind.first >= 0 && ind.first < static_cast<int>(m_field.width()) &&
				ind.second >= 0 && ind.second < static_cast<int>(m_field.height()))
			{
				return m_field.at(ind.first, ind.second);
			}
		}break;

		case 1:{
			int cx = int((x - m_offset.x()) / m_cellSize.x());
			int cy = int((y - m_offset.y()) / m_cellSize.y());
			
			const number vx = (x - (static_cast<number>(cx) * m_cellSize.x() + m_offset.x()))
							  / m_cellSize.x();
			const number vy = (y - (static_cast<number>(cy) * m_cellSize.y() + m_offset.y()))
							  / m_cellSize.y();
			const number ux = 1-vx;
			const number uy = 1-vy;

			if(		(cx >= 0)
				&&	(cy >= 0)
				&&	(cx + 1 < static_cast<int>(m_field.width()))
				&&	(cy + 1 < static_cast<int>(m_field.height())) )
			{

				return	ux*uy*m_field.at(cx, cy) + vx*uy*m_field.at(cx+1,cy)
					  + ux*vy*m_field.at(cx, cy+1) + vx*vy*m_field.at(cx+1, cy+1);
			}
			else{
				const int x0 = max(cx, 0);
				const int x1 = min(cx + 1, static_cast<int>(m_field.width()) - 1);
				const int y0 = max(cy, 0);
				const int y1 = min(cy + 1, static_cast<int>(m_field.height()) - 1);

				if(x0 >= static_cast<int>(m_field.width()) || x1 < 0 || y0 >= static_cast<int>(m_field.height()) || y1 < 0)
					return m_noDataValue;

				return	ux*uy*m_field.at(x0, y0) + vx*uy*m_field.at(x1,y0)
					  + ux*vy*m_field.at(x0, y1) + vx*vy*m_field.at(x1, y1);

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
	constexpr number roundOffset = 0.5;

	pair<int, int> c;
	if(m_cellSize.x() != 0)
		c.first = static_cast<int>(roundOffset + (x - m_offset.x()) / m_cellSize.x());
	else
		c.first = 0;
	
	if(m_cellSize.y() != 0)
		c.second = static_cast<int>(roundOffset + (y - m_offset.y()) / m_cellSize.y());
	else
		c.second = 0;
	return c;
}

vector2 Heightfield::index_to_coordinate(int ix, int iy) const
{
	return vector2(	m_offset.x() + static_cast<number>(ix) * m_cellSize.x(),
					m_offset.y() + static_cast<number>(iy) * m_cellSize.y());
}

vector2 Heightfield::extent() const
{
	return vector2(m_cellSize.x() * static_cast<number>(m_field.width()),
				   m_cellSize.y() * static_cast<number>(m_field.height()));
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
	std::string name = FindFileInStandardPaths(filename);
	UG_COND_THROW(name.empty(), "Couldn't locate file " << filename);

	FileReaderASC reader;
	reader.set_field(&hfield.field());
	reader.load_file(name.c_str());
	hfield.set_cell_size(vector2(reader.cell_size(), reader.cell_size()));
	hfield.set_offset(vector2(reader.lower_left_corner_x(),
							  reader.lower_left_corner_y()));
	hfield.set_no_data_value(reader.no_data_value());
}

}//	end of namespace
