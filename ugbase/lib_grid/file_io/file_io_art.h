//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m12 d01

#ifndef __H__LIB_GRID__FILE_IO_ART__
#define __H__LIB_GRID__FILE_IO_ART__

#include "../grid/grid.h"
#include "../common_attachments.h"

namespace ug
{
////////////////////////////////////////////////////////////////////////
//	predeclarations
class SubsetHandler;

////////////////////////////////////////////////////////////////////////
///	loads a grid from art
bool LoadGridFromART(Grid&grid, const char* filename,
					 SubsetHandler* pSH = NULL,
					 AVector3& aPos = aPosition);
/*
////////////////////////////////////////////////////////////////////////
///	saves a grid to art
bool SaveGridToART(Grid& grid, const char* filename,
				   SubsetHandler* pSH = NULL,
				   AVector3& aPos = aPosition);
*/
};

#endif
