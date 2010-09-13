//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m09 d10

#include "lib_grid/lg_base.h"

#ifndef __H__LIB_GRID__FILE_IO_MSH__
#define __H__LIB_GRID__FILE_IO_MSH__

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
