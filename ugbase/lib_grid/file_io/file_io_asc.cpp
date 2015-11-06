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

#include <string>
#include <fstream>
#include "file_io_asc.h"
#include "common/error.h"
#include "common/util/string_util.h"

using namespace std;

namespace ug{

FileReaderASC::FileReaderASC() :
	m_cellSize(0),
	m_llcorner(0, 0),
	m_noDataValue(0)
{
	m_privateField = make_sp(new Field<number>);
	m_field = m_privateField.get();
}

FileReaderASC::~FileReaderASC()
{
}

void FileReaderASC::set_field(Field<number>* field)
{
	if(field)
		m_field = field;
	else
		m_field = m_privateField.get();
}


void FileReaderASC::load_file(const char* filename)
{
	ifstream in(filename);
	UG_COND_THROW(!in, "Couldn't access file " << filename);

	int numCols = 0, numRows = 0;
	
//	indicate whether the lower left corner was specified through cell-center
	bool llcenterx = false;
	bool llcentery = false;

//	parse header
	for(int i = 0; i < 6; ++i){
		string name;
		double value;
		in >> name >> value;
		UG_COND_THROW(!in, "Couldn't parse expected name-value pair in row " << i);

		name = ToLower(name);

		if(name.compare("ncols") == 0)
			numCols = (int)value;
		else if(name.compare("nrows") == 0)
			numRows = (int)value;
		else if(name.compare("xllcenter") == 0){
			m_llcorner.x() = value;
			llcenterx = true;
		}
		else if(name.compare("xllcorner") == 0)
			m_llcorner.x() = value;
		else if(name.compare("yllcenter") == 0){
			m_llcorner.y() = value;
			llcentery = true;
		}
		else if(name.compare("yllcorner") == 0)
			m_llcorner.y() = value;
		else if(name.compare("cellsize") == 0)
			m_cellSize = value;
		else if(name.compare("nodata_value") == 0)
			m_noDataValue = value;
		else{
			UG_THROW("unknown identifier in header: " << name);
		}
	}
	
	if(llcenterx)
		m_llcorner.x() -= 0.5 * m_cellSize;
	if(llcentery)
		m_llcorner.y() -= 0.5 * m_cellSize;

//	parse values
	Field<number>& field = *m_field;
	field.resize_no_copy(numCols, numRows);
	for(int irow = numRows - 1; irow >= 0; --irow){
		for(int icol = 0; icol < numCols; ++icol){
			in >> field.at(icol, irow);
			UG_COND_THROW(!in, "Couldn't read value at col: " << icol << ", row: " << numRows - irow);
		}
	}
}

bool LoadGridFromASC(Grid& grid, const char* filename, AVector3& aPos)
{
	if(!grid.has_vertex_attachment(aPos))
		grid.attach_to_vertices(aPos);

	
	FileReaderASC reader;
	reader.load_file(filename);

	vector2 offset(reader.lower_left_corner_x(), reader.lower_left_corner_y());

	Grid::VertexAttachmentAccessor<AVector3> aaPos(grid, aPos);
	CreateGridFromField(grid, reader.field(),
						vector2(reader.cell_size(), reader.cell_size()),
						offset,
						reader.no_data_value(),
						aaPos);
	return true;
}

}//	end of namespace
