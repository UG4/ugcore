//	created by Sebastian Reiter
//	y09 m08 d03
//	s.b.reiter@googlemail.com

#ifndef __LIBGRID__FILE_IO_STL__
#define __LIBGRID__FILE_IO_STL__

#include "lib_grid/grid/grid.h"
#include "lib_grid/tools/subset_handler_interface.h"
#include "lib_grid/common_attachments.h"

namespace ug
{

///	loads stl-ascii-files.
/**
 * if aPosition is not attached to the vertices of the grid,
 * it will be automatically attached.
 * If however aNormal is not attached to the faces of the grid,
 * it will be ignored.
 */
bool LoadGridFromSTL(Grid& grid, const char* filename,
					ISubsetHandler* pSH = NULL,
					AVector3& aPos = aPosition,
					AVector3& aNormFACE = aNormal);

bool STLFileHasASCIIFormat(const char* filename);

bool LoadGridFromSTL_ASCII(Grid& grid, const char* filename,
						   ISubsetHandler* pSH = NULL,
						   AVector3& aPos = aPosition,
						   AVector3& aNormFACE = aNormal);

bool LoadGridFromSTL_BINARY(Grid& grid, const char* filename,
							ISubsetHandler* pSH = NULL,
							AVector3& aPos = aPosition,
							AVector3& aNormFACE = aNormal);

bool SaveGridToSTL(Grid& grid, const char* filename,
					ISubsetHandler* pSH = NULL,
					AVector3& aPos = aPosition);

}

#endif
