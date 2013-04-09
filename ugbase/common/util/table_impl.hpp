// created by Sebastian Reiter
// s.b.reiter@googlemail.com

#ifndef __H__TABLE_IMPL__
#define __H__TABLE_IMPL__

#include <iomanip>
#include <cassert>

namespace ug{

template <class T>
Table<T>::Table() :
	m_numRows(0),
	m_numCols(0)
{
}

template <class T>
Table<T>::Table(size_t numRows, size_t numCols) :
	m_numRows(0),
	m_numCols(0)
{
	add_rows(numRows);
	add_cols(numCols);
}

template <class T>
Table<T>::~Table()
{
	clear();
}

template <class T>
void Table<T>::clear()
{
	for(size_t i_row = 0; i_row < m_numRows; ++i_row){
		for(size_t i_col = 0; i_col < m_numCols; ++i_col){
			delete m_data[i_row][i_col];
		}
	}
	m_numRows = m_numCols = 0;
	m_data.clear();
}

template <class T>
void Table<T>::add_rows(size_t num)
{
	if(num > 0){
		size_t newSize = m_numRows + num;
		m_data.resize(newSize);

		if(m_numCols > 0){
			for(size_t i_row = m_numRows; i_row < newSize; ++i_row){
				m_data[i_row].resize(m_numCols);
				for(size_t i_col = 0; i_col < m_numCols; ++i_col)
					m_data[i_row][i_col] = new T;
			}
		}
		
		m_numRows = newSize;
	}
}

template <class T>
void Table<T>::add_cols(size_t num)
{
	if(num > 0){
		size_t newSize = m_numCols + num;
		for(size_t i_row = 0; i_row < m_numRows; ++i_row){
			m_data[i_row].resize(newSize);
			for(size_t i_col = m_numCols; i_col < newSize; ++i_col)
				m_data[i_row][i_col] = new T;
		}
		
		m_numCols = newSize;
	}
}

template <class T>
T& Table<T>::operator() (size_t rowInd, size_t colInd)
{
	if(rowInd >= num_rows())
		add_rows(rowInd - num_rows() + 1);

	if(colInd >= num_cols())
		add_cols(colInd - num_cols() + 1);

	return *m_data[rowInd][colInd];
}

template <class T>
const T& Table<T>::operator() (size_t rowInd, size_t colInd) const
{
	assert(rowInd < num_rows());
	assert(colInd < num_cols());
	return *m_data[rowInd][colInd];
}

template <class T>
size_t Table<T>::num_rows() const
{
	return m_numRows;
}

template <class T>
size_t Table<T>::num_cols() const
{
	return m_numCols;
}



template <class T>
std::ostream& operator << (std::ostream& os, const Table<T>& table)
{
//	we'll fill a fresh table with strings of the given table
//	At the same time we'll find the max-width of each column
	Table<std::string>	strTable(table.num_rows(), table.num_cols());
	std::vector<size_t> colSizes(table.num_cols(), 0);
	
	for(size_t i_row = 0; i_row < table.num_rows(); ++i_row){
		for(size_t i_col = 0; i_col < table.num_cols(); ++i_col){
			strTable(i_row, i_col) = EntryToString(table, i_row, i_col);
			colSizes[i_col] = std::max(colSizes[i_col], strTable(i_row, i_col).size());
		}
	}
	
//	now print each row
	for(size_t i_row = 0; i_row < table.num_rows(); ++i_row){
		for(size_t i_col = 0; i_col < table.num_cols(); ++i_col){
			os << std::setw(colSizes[i_col] + 2) << std::right << strTable(i_row, i_col);
		}
		os << std::endl;
	}

	return os;
}



template <class T>
std::string EntryToString(const Table<T>& table, size_t rowInd, size_t colInd)
{
	std::stringstream ss;
	ss << table(rowInd, colInd);
	return ss.str();
}

}//	end of namespace

#endif
