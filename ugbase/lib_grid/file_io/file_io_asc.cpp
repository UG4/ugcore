// created by Sebastian Reiter
// s.b.reiter@gmail.com

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
	UG_COND_THROW(!in, "Couldn't read from file " << filename);

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
			UG_COND_THROW(!in, "Couln't read value at col: " << icol << ", row: " << numRows - irow);
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

void LoadHeightfieldFromASC(Heightfield& hfield, const char* filename)
{
	FileReaderASC reader;
	reader.set_field(&hfield.field());
	reader.load_file(filename);
	hfield.set_cell_size(vector2(reader.cell_size(), reader.cell_size()));
	hfield.set_offset(vector2(reader.lower_left_corner_x(),
							  reader.lower_left_corner_y()));
	hfield.set_no_data_value(reader.no_data_value());
}

}//	end of namespace
