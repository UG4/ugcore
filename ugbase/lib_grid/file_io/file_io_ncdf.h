//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m07 d09

#ifndef __H__LIBGRID__FILE_IO_NCDF__
#define __H__LIBGRID__FILE_IO_NCDF__

#include "lib_grid/lg_base.h"

namespace ug
{
/**
 * Saves a grid to NCDF-ASCII-format (EXODUS).
 */
bool SaveGridToNCDF(Grid& grid, const char* filename,
					SubsetHandler* pSH,
					APosition aPos = aPosition);


}//	end of namespace

#endif
