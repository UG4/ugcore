// created by Sebastian Reiter
// s.b.reiter@gmail.com
// 13.11.2012 (d,m,y)

#include "table.h"

namespace ug{

std::string EntryToString(const Table<std::string>& table,
						  size_t rowInd, size_t colInd)
{
	return table(rowInd, colInd);
}

std::string EntryToString(const Table<std::stringstream>& table,
						  size_t rowInd, size_t colInd)
{
	return table(rowInd, colInd).str();
}

}// end of namespace
