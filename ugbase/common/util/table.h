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

#ifndef __H__TABLE
#define __H__TABLE

#include <sstream>
#include <string>
#include <vector>
#include <iostream>

namespace ug{

/// \addtogroup ugbase_common_types
/// \{

///	Useful for printing a table to the terminal or into a file.
/**	This class is for output purposes only. It automatically adds spacings, so
 * that all entries in a column are indeed displayed in the same column.
 * Note that the implementation is not optimized in regards to speed and efficiency.
 * The class thus shouldn't be used in performance critical sections of your code.
 * If it is however used in output related code, the introduced overhead should
 * be fine and probably even negligible.
 *
 * The class is most commonly used as a string table (<tt>Table\<std::string\></tt>) or
 * as a stringstream table (<tt>Table\<std::stringstream\></tt>). If you want to use it
 * for your own types, you may specialize the template method
 * \code
 * template <class T> std::string EntryToString(const Table<T>& table, size_t rowInd, size_t colInd);
 * \endcode
 * The default implementation of EntryToString may already be suited for most cases.
 *
 * Here's an example on how to use the class:
 * \code
 * #include "common/util/table.h"
 * ...
 * ug::Table<std::stringstream> table;
 * table(0, 0) << "num rows:";		table(0, 1) << table.num_rows();
 * table(1, 0) << "num columns:";	table(1, 1) << table.num_cols();
 * std::cout << table;
 * \endcode
 *
 * And this is what the output looks like:
 * \verbatim
     num rows:  2
     num columns:  2
 * \endverbatim
 *
 * Note that several typedefs exist: StringTable, StringStreamTable
 * \todo	different alignments for different columns / rows / fields
 */

template <class T>
class Table
{
	public:
		Table();
		Table(size_t numRows, size_t numCols);

		~Table();

		void clear();

		void add_rows(size_t num);
		void add_cols(size_t num);

	///	Returns a reference to the given entry.
	/**	If an entry lies outside of the tables bounds, the table will
	 * automatically be resized accordingly.*/
		T& operator() (size_t rowInd, size_t colInd);

		/// uses operator() to set an entry to a value
		void set(size_t rowInd, size_t colInd, T value)
		{
			operator()(rowInd, colInd) = value;
		}

	///	Returns a reference to the given entry.
		const T& operator() (size_t rowInd, size_t colInd) const;

	/// uses operator() to get a value
		const T &get(size_t rowInd, size_t colInd) const
		{
			return operator()(rowInd, colInd);
		}

		size_t num_rows() const;
		size_t num_cols() const;
		



		void set_default_row_seperator(const char *c)
		{
			m_defaultRowSeperator = *c;
		}
		void set_row_seperator(size_t i_row, const char *c)
		{
			if(m_rowSep.size() <= i_row) m_rowSep.resize(i_row+1, 0x00);
			m_rowSep[i_row] = *c;
		}
		void set_row_seperators(std::string s)
		{
			if(m_rowSep.size() < s.size()) m_rowSep.resize(s.size(), 0x00);
			for(size_t i=0; i<s.size(); i++)
				m_rowSep[i] = s[i];
		}
		void set_default_col_seperator(const char *c)
		{

			m_defaultColSeperator = *c;
		}
		void set_col_seperator(size_t i_col, const char *c)
		{
			if(m_colSep.size() <= i_col) m_colSep.resize(i_col+1, 0x00);
			m_colSep[i_col] = *c;
		}
		void set_col_seperators(std::string s)
		{
			if(m_colSep.size() < s.size()) m_colSep.resize(s.size(), 0x00);
			for(size_t i=0; i<s.size(); i++)
				m_colSep[i] = s[i];
		}

		void set_default_col_alignment(const char *c)
		{
			m_defaultColAlignment = *c;
		}
		void set_col_alignment(size_t i_col, const char *c)
		{
			if(m_colAlign.size() <= i_col) m_colAlign.resize(i_col+1, 0x00);
			m_colAlign[i_col] = *c;
		}
		void set_col_alignments(std::string s)
		{
			if(m_colAlign.size() < s.size()) m_colAlign.resize(s.size(), 0x00);
			for(size_t i=0; i<s.size(); i++)
				m_colAlign[i] = s[i];
		}

		std::ostream& stream(std::ostream& os) const;
		std::string to_latex() const;
		std::string to_string() const;
		std::string to_csv(const char *seperator) const;

		char get_row_sep(size_t row) const
		{
			if(row >= m_rowSep.size()) return m_defaultRowSeperator;
			return m_rowSep[row] != 0x00 ? m_rowSep[row] : m_defaultRowSeperator;
		}

		char get_col_sep(size_t col) const
		{
			if(col >= m_colSep.size()) return m_defaultColSeperator;
			return m_colSep[col] != 0x00 ? m_colSep[col] : m_defaultColSeperator;
		}

		char get_col_alignment(size_t col) const
		{
			if(col >= m_colAlign.size()) return m_defaultColAlignment;
			return m_colAlign[col] != 0x00 ? m_colAlign[col] : m_defaultColAlignment;
		}

		void transpose()
		{
			DataVec newData;
			newData.resize(m_numCols);
			for(size_t irow = 0; irow < m_numRows; ++irow){
				for(size_t icol = 0; icol < m_numCols; ++icol){
					newData[icol].push_back(m_data[irow][icol]);
				}
			}
			std::swap(m_numRows, m_numCols);
			m_data.swap(newData);
		}
	private:
		size_t m_numRows;
		size_t m_numCols;

		typedef std::vector<std::vector<T*> > DataVec;
		DataVec	m_data;

		std::vector<char> m_colAlign;
		std::vector<char> m_colSep;
		std::vector<char> m_rowSep;
		char m_defaultColSeperator, m_defaultRowSeperator, m_defaultColAlignment;

};

///	prints a table to the specified ostream.
/** Assumes that output begins at a new line.
 * The method uses the template method
 * size_t EstimateEntryWidth<T>(const Table<T>& table, size_t rowInd, size_t colInd).
 * You may overload this method to adjust it to your datatypes.
 * Note that specializations exist for std::string and std::stringstream.
 *
 * \todo	Support for line-breaks ('\n' etc) and multi-line rows.
 */
template <class T>
std::ostream& operator << (std::ostream& os, const Table<T>& table);


///	Returns a string-representation of the current entry.
/**	The default implementation prints the entry to a std::stringstream and returns
 * the associated string. Specializations exist for std::string and std::stringstream.
 */
template <class T>
std::string EntryToString(const Table<T>& table, size_t rowInd, size_t colInd);

inline
std::string EntryToString(const Table<std::string>& table,
						  size_t rowInd, size_t colInd);

inline
std::string EntryToString(const Table<std::stringstream>& table,
						  size_t rowInd, size_t colInd);


typedef Table<std::string>	StringTable;
typedef Table<std::stringstream> StringStreamTable;

// end group ugbase_common_types
/// \}

}//	end of namespace

////////////////////////////////////////
//	include implementation
#include "table_impl.hpp"

#endif
