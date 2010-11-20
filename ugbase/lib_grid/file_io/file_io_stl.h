//	created by Sebastian Reiter
//	y09 m08 d03
//	s.b.reiter@googlemail.com

#ifndef __LIBGRID__FILE_IO_STL__
#define __LIBGRID__FILE_IO_STL__

#include "lib_grid/lg_base.h"

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

}

#endif
