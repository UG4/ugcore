/**
 * \file amg_debug.h
 *
 * \author Martin Rupp
 *
 * \date 16.06.10
 *
 * Goethe-Center for Scientific Computing 2009-2010.
 */



#ifndef __H__UG__LIB_DISC__AMG__AMG_DEBUG_H__
#define __H__UG__LIB_DISC__AMG__AMG_DEBUG_H__

#include <fstream>
#include "amg_debug_helper.h"
#include "lib_algebra/common/graph/graph.h"
#include "lib_algebra/lib_algebra.h"
#include "amg_profiling.h"
#include "postscript.h"

#include "lib_algebra/common/connection_viewer_output.h"

namespace ug {
// WriteToFile
//--------------------------------------------------
//! writes to a file in somewhat SparseMatrix-market format (for connection viewer)
template<typename T>
void AMGWriteToFile(const SparseMatrix<T> &A, int fromlevel, int tolevel, const char *filename, const cAMG_helper &h)
{
	AMG_PROFILE_FUNC();
	if(h.has_positions() == false)
	{
		UG_LOG("AWriteToFile not possible: no positions available.")
		return;
	}

	std::fstream file(filename, std::ios::out);
	file << 1 << std::endl; // connection viewer version

	int minlevel = std::min(fromlevel, tolevel);
	WritePositionsToStream(file, h.positions[minlevel], h.dimension);

	file << 1 << std::endl;
	for(size_t i=0; i < A.num_rows(); i++)
	{
		for(typename SparseMatrix<T>::const_row_iterator conn = A.begin_row(i); conn != A.end_row(i); ++conn)
			if(conn.value() != 0.0)
				file << h.GetOriginalIndex(tolevel, minlevel, i) << " " << h.GetOriginalIndex(fromlevel, minlevel, conn.index()) << " " << conn.value() << std::endl;
	}
}

// writeToFile
//--------------------------------------------------
//! writes to a file in somewhat SparseMatrix-market format (for connection viewer)
template<typename T>
void AMGWriteToFilePS(const SparseMatrix<T> &A, int fromlevel, int tolevel, const char *filename, const cAMG_helper &h)
{
	AMG_PROFILE_FUNC();
	if(h.has_positions() == false)
	{
		UG_LOG("AWriteToFilePS not possible: no positions available.")
		return;
	}

	postscript ps;
	ps.create(filename);
	int minlevel = std::min(fromlevel, tolevel);
	for(size_t i=0; i < A.num_rows(); i++)
	{
		int from = h.GetOriginalIndex(tolevel, minlevel, i);
		ps.move_to(h.positions[minlevel][from].x, h.positions[minlevel][from].y);
		ps.print_text( std::string("0") + ToString(i));

		for(typename SparseMatrix<T>::const_row_iterator conn = A.begin_row(i); conn != A.end_row(i); ++conn)
		{
			if(conn.value() != 0.0)
			{
				if(conn.index() != i)
				{
					int to = h.GetOriginalIndex(fromlevel, minlevel, conn.index());
					ps.move_to(h.positions[minlevel][from].x, h.positions[minlevel][from].y);
					ps.line_to(h.positions[minlevel][to].x, h.positions[minlevel][to].y);

				}
			}
		}
	}

	std::cout << std::endl;
}

// could be in cpp
inline void WriteAMGGraphToFile(cgraph &G, const char *filename, const cAMG_helper &h, int level)
{
	AMG_PROFILE_FUNC();
	if(h.has_positions() == false)
	{
		UG_LOG("WriteAMGGraphToFile not possible: no positions available.")
		return;
	}

	std::fstream file(filename, std::ios::out);
	file << /*CONNECTION_VIEWER_VERSION*/ 1 << std::endl;

	WritePositionsToStream(file, h.positions[level], h.dimension);
	file << 1 << std::endl;
	for(size_t i=0; i < G.size(); i++)
	{
		for(cgraph::const_row_iterator it = G.begin_row(i); it != G.end_row(i); ++it)
			file << i << " " << (*it) << "  " << std::endl;
	}
}

}



#endif // __H__UG__LIB_DISC__AMG_SOLVER__AMG_DEBUG_H__
