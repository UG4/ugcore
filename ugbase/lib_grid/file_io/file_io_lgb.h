//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m11 d11

#ifndef __H__LIBGRID__FILE_IO_BIN__
#define __H__LIBGRID__FILE_IO_BIN__

#include "lib_grid/grid/grid.h"
#include "lib_grid/tools/subset_handler_interface.h"
#include "lib_grid/common_attachments.h"

namespace ug
{
/**
 * Saves a grid to LibGridBinary-format.
 * Awaits a list of subset-handler-pointers and the number
 * of subset-handlers that shall be written.
 */
bool SaveGridToLGB(Grid& grid, const char* filename,
				   ISubsetHandler** ppSH, int numSHs,
				   APosition aPos = aPosition);


/**
 * Loads a grid from LibGridBinary-format.
 * Awaits a list of subset-handler-pointers and the number
 * of subset-handlers that shall be read.
 * Make sure that all passed subset-handlers are already registered
 * at the grid.
 */
bool LoadGridFromLGB(Grid& grid, const char* filename,
				   ISubsetHandler** ppSH, int numSHs,
				   APosition aPos = aPosition);

}//	end of namespace

#endif
