//	created by Sebastian Reiter
//	y09 m08 d03
//	s.b.reiter@googlemail.com

#ifndef __LIBGRID__FILE_IO_DUMP__
#define __LIBGRID__FILE_IO_DUMP__

#include "lib_grid/grid/grid.h"
#include "lib_grid/tools/subset_handler_interface.h"
#include "lib_grid/common_attachments.h"

namespace ug
{

///	loads dump-files.
/**
 * A dump-file is a very simple format that describes a triangular grid.
 * The format is an ASCII format.
 * a # at the beginning of a line marks a comment.
 * other lines should either be empty or should describe a triangle.
 * a triangle is described be three coordinate-triples. Space is used as
 * seperator. a file could look like this:
 *
 * # start of .dump file
 * # a single triangle
 *
 * 0 0 0 1 0 0 0 1 0
 *
 * # end of .dump file
 *
 */
bool LoadGridFromDUMP(Grid& grid, const char* filename,
					ISubsetHandler* pSH = NULL, AVector3& aPos = aPosition);

}

#endif
