//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m11 d13

#include "../grid/grid.h"
#include "../common_attachments.h"

#ifndef __H__LIB_GRID__FILE_IO_TXT__
#define __H__LIB_GRID__FILE_IO_TXT__

namespace ug
{
////////////////////////////////////////////////////////////////////////
///	loads a grid from txt
bool LoadGridFromTXT(Grid&grid, const char* filename, AVector3& aPos = aPosition);

////////////////////////////////////////////////////////////////////////
///	saves a grid to txt
bool SaveGridToTXT(Grid& grid, const char* filename, AVector3& aPos = aPosition);





};

#endif
