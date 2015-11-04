#ifndef __H__LIB_GRID__FILE_IO_ART__
#define __H__LIB_GRID__FILE_IO_ART__

#include "../grid/grid.h"
#include "lib_grid/subset_handler.h"
#include "../common_attachments.h"

namespace ug
{
////////////////////////////////////////////////////////////////////////
///	loads a grid from art
bool LoadGridFromART(Grid&grid, const char* filename,
					 ISubsetHandler* pSH = NULL,
					 AVector3& aPos = aPosition);

////////////////////////////////////////////////////////////////////////
///	saves a grid to art
bool SaveGridToART(Grid& grid, const char* filename,
				   ISubsetHandler* pSH = NULL,
				   AVector3& aPos = aPosition);

};

#endif
