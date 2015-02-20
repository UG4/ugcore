// created by Sebastian Reiter
// s.b.reiter@gmail.com

#include <string>
#include <fstream>
#include "file_io_asc.h"
#include "common/error.h"
#include "common/util/string_util.h"
#include "lib_grid/algorithms/field_util.h"

using namespace std;

namespace ug{

FileReaderASC::FileReaderASC() :
	m_cellSize(0),
	m_center(0, 0),
	m_noDataValue(0),
	m_numRows(0),
	m_numCols(0)
{
}

FileReaderASC::~FileReaderASC()
{
}


void FileReaderASC::load_file(const char* filename)
{
	ifstream in(filename);
	UG_COND_THROW(!in, "Couldn't read from file " << filename);

//	parse header
	for(int i = 0; i < 6; ++i){
		string name;
		double value;
		in >> name >> value;
		UG_COND_THROW(!in, "Couldn't parse expected name-value pair in row " << i);

		UG_LOG("name-value-pair: " << name << ", " << value << endl);
		name = ToLower(name);

		if(name.compare("ncols") == 0)
			m_numCols = (size_t)value;
		else if(name.compare("nrows") == 0)
			m_numRows = (size_t)value;
		else if(name.compare("xllcenter") == 0)
			m_center.x() = value;
		else if(name.compare("yllcenter") == 0)
			m_center.y() = value;
		else if(name.compare("cellsize") == 0)
			m_cellSize = value;
		else if(name.compare("nodata_value") == 0)
			m_noDataValue = value;
		else{
			UG_THROW("unknown identifier in header: " << name);
		}
	}
	
//	parse values
	m_field.resize_no_copy(m_numCols, m_numRows);
	for(size_t irow = 0; irow < m_numRows; ++irow){
		for(size_t icol = 0; icol < m_numCols; ++icol){
			in >> m_field.at(icol, irow);
			UG_COND_THROW(!in, "Couln't read value at col: " << icol << ", row: " << irow);
		}
	}
	UG_LOG("Field parsed.\n");
}

bool LoadGridFromASC(Grid& grid, const char* filename, AVector3& aPos)
{
	if(!grid.has_vertex_attachment(aPos))
		grid.attach_to_vertices(aPos);

	
	FileReaderASC reader;
	reader.load_file(filename);

	Grid::VertexAttachmentAccessor<AVector3> aaPos(grid, aPos);
	CreateGridFromField(grid, reader.field(),
						vector2(reader.cell_size(), reader.cell_size()),
						//vector2(reader.center_x(), reader.center_y()),
						vector2(0, 0),
						reader.no_data_value(),
						aaPos);
	return true;
}

}//	end of namespace
