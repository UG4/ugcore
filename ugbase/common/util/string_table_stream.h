/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
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

#ifndef STRING_TABLE_STREAM_H_
#define STRING_TABLE_STREAM_H_

#include <vector>

#include "table.h"
#include "string_util.h"
#include "stringify.h"

/// \addtogroup ugbase_common_types
/// \{
///	Useful for printing a table like a stream to the terminal or into a file.
/**
 * Here's an example on how to use the class:
 * \code
 * #include "common/util/string_table_stream.h"
 * ...
 * ug::StringTableStream sts;
 * sts << "|" << 31.3 << ", " <<  0 << "|\n";
 * sts << "|" << 0 << ", " << 768 << "|\n";
 * std::cout << table;
 * \endcode
 *
 * And this is what the output looks like:
 * \verbatim
  | 31.3 , 0   |
  | 0    , 768 |
 * \endverbatim
 *
 * This class uses internally StringTable.
 * @sa StringTable, Table
 */

namespace ug {

class StringTableStream
{
private:
	StringTable s;
	size_t m_curRow;
	size_t m_curCol;

public:
	StringTableStream()
	{
		m_curRow=m_curCol = 0;
	}

	StringTable &table() { return s; }

	/**
	 * operator <<. sets curCol++
	 * @param t anything which is ostream << -able
	 * @return *this
	 */
	template<typename T>
	StringTableStream &operator << (T t)
	{
		s(m_curRow, m_curCol++) = ToString(t);
		return *this;
	}

	/**
	 * operator <<
	 * @param c a const char *. If contains a \n, sets curRow++, curCol=0
	 * @return *this
	 */
	StringTableStream &operator << (const char *c)
	{
		size_t l = strlen(c);
//		size_t begin=0;
		std::vector<char> vc;
		for(size_t i=0; i<l; i++)
		{
			if(c[i] == '\n')
			{
				vc.push_back(0);
				s(m_curRow, m_curCol) = &vc[0];
				vc.clear();
				m_curRow++;
				m_curCol=0;
				if(i == l-1) return *this;
			}
			else
				vc.push_back(c[i]);
		}
		vc.push_back(0);
		s(m_curRow, m_curCol++) = &vc[0];
		return *this;
	}



	struct RepeatedCol
	{
		RepeatedCol(const std::string &_content, size_t _number) :
			content(_content), number(_number) {}
		std::string content;
		size_t number;
	};

	StringTableStream &operator << (RepeatedCol c)
	{
		for(size_t i=0; i<c.number; i++)
			*this << c.content;
		return *this;
	}


	/**
	 * \code
	 * ug::StringTableStream sts;
	 * sts << 1 << "\n";
	 * sts << sts.emtpy_col(1) << 1 << "\n";
	 * sts << sts.empty_col(2) << 1 << "\n";
	 * \endcode
	 * prints
	 * \verbatim
	 * 1
	 *   1
	 *     1
	 * \endverbatim
	 * @param nrRepeat			number of times to repeat
	 * @param contentToRepeat	what to repeat
	 * @return RepeatedCol-struct for <<-ing into sts
	 */
	RepeatedCol empty_col(size_t i)
	{
		return RepeatedCol("", i);
	}

	/**
	 * \code
	 * ug::StringTableStream sts;
	 * sts << 1 << sts.repeat_col(2, 0) << "\n";
	 * sts << 0 << 1 << 0 << "\n";
	 * sts << sts.repeat(2, 0) << 1 << "\n";
	 * \endcode
	 * prints
	 * \verbatim
	 * 1 0 0
	 * 0 1 0
	 * 0 0 1
	 * \endverbatim
	 * @param nrRepeat			number of times to repeat
	 * @param contentToRepeat	what to repeat
	 * @return RepeatedCol-struct for <<-ing into sts
	 */
	template<typename T>
	RepeatedCol repeat_col(size_t nrRepeat, T contentToRepeat)
	{
		return RepeatedCol(ToString(contentToRepeat), nrRepeat);
	}

	template<typename T>
	void set(size_t r, size_t c, T t)
	{
		s(r, c) = ToString(t);
	}

	StringTableStream &operator << (std::string str)
	{
		*this << str.c_str();
		return *this;
	}
	std::string to_string() const
	{
		return s.to_string();
	}

	void clear()
	{
		m_curRow=m_curCol = 0;
		s.clear();
	}
};

inline std::ostream& operator << (std::ostream& os, const StringTableStream &sts)
{
	os << sts.to_string();
	return os;
}

}
#endif
