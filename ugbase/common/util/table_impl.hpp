/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__TABLE_IMPL__
#define __H__TABLE_IMPL__

#include <iomanip>
#include <cassert>
#include "common/error.h"
#include "common/log.h"
#include "string_util.h"

namespace ug{

template <class T>
Table<T>::Table() :
	m_numRows(0),
	m_numCols(0)
{
	m_defaultColSeperator = ' ';
	m_defaultRowSeperator = ' ';
	m_defaultColAlignment = 'l';
}

template <class T>
Table<T>::Table(size_t numRows, size_t numCols) :
	m_numRows(0),
	m_numCols(0)
{
	m_defaultColSeperator = ' ';
	m_defaultRowSeperator = ' ';
	m_defaultColAlignment = 'l';
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
		add_rows((rowInd + 1) - num_rows());

	if(colInd >= num_cols())
		add_cols((colInd + 1) - num_cols());

	return *m_data[rowInd][colInd];
}

template <class T>
const T& Table<T>::operator() (size_t rowInd, size_t colInd) const
{
	UG_COND_THROW(rowInd >= num_rows(), "Bad row index: " << rowInd << "! Only " << num_rows() << " rows exist.");
	UG_COND_THROW(colInd >= num_cols(), "Bad col index: " << colInd << "! Only " << num_cols() << " cols exist.");
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
std::string Table<T>::to_string() const
{
	std::stringstream ssOut;
	ssOut << *this;
	return ssOut.str();
}



template<typename T>
std::string Table<T>::to_latex() const
{
	Table<std::string>	strTable(num_rows(), num_cols());

	std::stringstream os;
	os << "\\begin{table}\n\\centering\n";
	os << "\\begin{tabular}{";
	for(size_t i=0; i < num_cols(); i++)
		os << get_col_sep(i) << get_col_alignment(i);
	os << get_col_sep(num_cols());
	os << "}\n";


	for(size_t i_row = 0; i_row < num_rows(); ++i_row)
	{
		if(get_row_sep(i_row) != ' ')
			os << "\\hline \n";

		for(size_t i_col = 0; i_col < num_cols(); ++i_col)
		{
			if(i_col != 0) os << " & ";
			os << EntryToString(*this, i_row, i_col);
		}

		os << " \\\\\n";
	}
	if(get_row_sep(num_rows())  != ' ')
		os << "\\hline \n";
	os << "\\end{tabular}\n";
	os << "\\end{table}\n";
	return os.str();
}


template<typename T>
std::string Table<T>::to_csv(const char *seperator) const
{
	std::stringstream os;
	for(size_t i_row = 0; i_row < num_rows(); ++i_row)
	{
		for(size_t i_col = 0; i_col < num_cols(); ++i_col)
		{
			if(i_col != 0) os << " " << seperator << " ";
			os << EntryToString(*this, i_row, i_col);
		}

		os << "\n";
	}
	return os.str();
}

template <class T>
std::ostream& operator << (std::ostream& os, const Table<T>& table)
{
	return table.stream(os);
}

template <class T>
std::ostream& Table<T>::stream(std::ostream& os) const
{
	const Table<T> &table = *this;
//	we'll fill a fresh table with strings of the given table
//	At the same time we'll find the max-width of each column
	Table<std::string>	strTable(num_rows(), num_cols());
	std::vector<size_t> colSizes(num_cols(), 0);

	for(size_t i_row = 0; i_row < num_rows(); ++i_row){
		for(size_t i_col = 0; i_col < num_cols(); ++i_col){
			strTable(i_row, i_col) = EntryToString(table, i_row, i_col);
			colSizes[i_col] = std::max(colSizes[i_col], strTable(i_row, i_col).size());
		}
	}

	size_t totalRowLength=0;
	for(size_t i_col = 0; i_col < num_cols(); ++i_col)
	{
		totalRowLength++;
		totalRowLength += colSizes[i_col] + 2;
	}
	totalRowLength++;

//	now print each row
	for(size_t i_row = 0; i_row < num_rows(); ++i_row)
	{
		if(get_row_sep(i_row) != 0x00 && get_row_sep(i_row) != ' ')
			os << repeat(get_row_sep(i_row), totalRowLength) << "\n";
		for(size_t i_col = 0; i_col < num_cols(); ++i_col)
		{
			os << get_col_sep(i_col);
			size_t l = colSizes[i_col] + 1;
			size_t s = strTable(i_row, i_col).size();
			os << " ";
			if(get_col_alignment(i_col) == 'r')
				os << std::setw(l) << std::right << strTable(i_row, i_col);
			else if(get_col_alignment(i_col) == 'l')
				os << std::setw(l) << std::left << strTable(i_row, i_col);
			else
				os << repeat(' ', (l-s)/2) << std::setw(l-(l-s)/2) << std::left << strTable(i_row, i_col);
		}
		os << get_col_sep(num_cols()) << std::endl;
	}
	if(get_row_sep(num_cols()) != 0x00 && get_row_sep(num_cols()) != ' ')
		os << repeat(get_row_sep(num_cols()), totalRowLength) << "\n";
	return os;
}




template <class T>
inline
std::string EntryToString(const Table<T>& table, size_t rowInd, size_t colInd)
{
	std::stringstream ss;
	ss << table(rowInd, colInd);
	return ss.str();
}


inline
std::string EntryToString(const Table<std::string>& table,
						  size_t rowInd, size_t colInd)
{
	return table(rowInd, colInd);
}

inline
std::string EntryToString(const Table<std::stringstream>& table,
						  size_t rowInd, size_t colInd)
{
	return table(rowInd, colInd).str();
}

}//	end of namespace

#endif
