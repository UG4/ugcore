// created by Sebastian Reiter
// s.b.reiter@gmail.com

#ifndef __H__UG_file_io_asc
#define __H__UG_file_io_asc

#include "common/util/field.h"
#include "common/math/ugmath_types.h"
#include "lib_grid/lg_base.h"

namespace ug{

class FileReaderASC{
	public:
		FileReaderASC();
		~FileReaderASC();

		void load_file(const char* filename);

		number cell_size() const		{return m_cellSize;}
		number center_x() const			{return m_center.x();}
		number center_y() const			{return m_center.y();}
		number no_data_value() const	{return m_noDataValue;}
		size_t num_rows() const			{return m_numRows;}
		size_t num_columns() const		{return m_numCols;}

		number at(size_t r, size_t c) const	{return m_field.at(r, c);}
		const Field<number>& field() const	{return m_field;}

	private:
		Field<number>	m_field;
		number			m_cellSize;
		vector2			m_center;
		number			m_noDataValue;
		size_t			m_numRows;
		size_t			m_numCols;
};

bool LoadGridFromASC(Grid& grid, const char* filename, AVector3& aPos = aPosition);

}//	end of namespace

#endif	//__H__UG_file_io_asc
