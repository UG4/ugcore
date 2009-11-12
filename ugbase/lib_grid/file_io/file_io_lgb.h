//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m11 d11

#ifndef __H__LIBGRID__FILE_IO_BIN__
#define __H__LIBGRID__FILE_IO_BIN__

#include "lib_grid/lg_base.h"

namespace ug
{
/**
 * Saves a grid to LibGridBinary-format.
 */
bool SaveGridToLGB(Grid& grid, const char* filename,
				   SubsetHandler* pSH = NULL, APosition aPos = aPosition);


/**
 * Loads a grid from LibGridBinary-format.
 */
bool LoadGridFromLGB(Grid& grid, const char* filename,
				   SubsetHandler* pSH = NULL, APosition aPos = aPosition);

}//	end of namespace

#endif
