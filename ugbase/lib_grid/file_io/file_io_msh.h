//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m09 d10

#ifndef __H__LIB_GRID__FILE_IO_MSH__
#define __H__LIB_GRID__FILE_IO_MSH__

#include "lib_grid/grid/grid.h"
#include "lib_grid/tools/subset_handler_interface.h"
#include "lib_grid/common_attachments.h"

namespace ug
{
////////////////////////////////////////////////////////////////////////
///	loads a grid from the GMSH ascii .msh format
/**	Please check the GMSH manual for syntax information. */
bool LoadGridFromMSH(Grid&grid, const char* filename,
					 ISubsetHandler* psh = NULL,
					 AVector3& aPos = aPosition);

};

#endif
