// created by Sebastian Reiter
// s.b.reiter@googlemail.com

#ifndef __H__TABLE
#define __H__TABLE

#include <sstream>
#include <string>
#include <vector>
#include <iostream>

namespace ug{
///	Useful for printing a table into the console or into a file.
/**	This class is for output purposes only. It automatically adds spacings, so
 * that all entries in a column are indeed displayed in the same column.
 * Note that the implementation is not optimized in regards to speed and efficiency.
 * The class thus shouldn't be used in performance critical sections of your code.
 * If it is however used in ouput related code, the introduced overhead should
 * be fine and probably even neglectible.
 *
 * The class is most commonly used as a string table (Table<std::string>) or
 * as a stringstream table (Table<std::stringstream>). If you want to use it
 * for your own types, you may specialize the template method
 * \code
 * template <class T> std::string EntryToString(const Table<T>& table,
 *												size_t rowInd, size_t colInd);
 * \endcode
 * The default implementation of EntryToString may already be suited for most cases.
 *
 * Here's an example on how to use the class:
 * \code
 * ug::Table<std::stringstream> table(2, 2);
 * table(0, 0) << "num rows:";		table(0, 1) << table.num_rows();
 * table(1, 0) << "num columns:";	table(1, 1) << table.num_cols();
 * \endcode
 *
 * \todo	different alignements for different columns / rows / fields
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
		
		T& operator() (size_t rowInd, size_t colInd);
		const T& operator() (size_t rowInd, size_t colInd) const;
		
		size_t num_rows() const;
		size_t num_cols() const;
		
	private:
		size_t m_numRows;
		size_t m_numCols;
		
		typedef std::vector<std::vector<T*> > DataVec;
		DataVec	m_data;
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

std::string EntryToString(const Table<std::string>& table,
						  size_t rowInd, size_t colInd);
std::string EntryToString(const Table<std::stringstream>& table,
						  size_t rowInd, size_t colInd);


}//	end of namespace

////////////////////////////////////////
//	include implementation
#include "table_impl.hpp"

#endif
