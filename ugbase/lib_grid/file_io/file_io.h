//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m11 d13

#include "file_io_txt.h"
#include "file_io_tetgen.h"
#include "file_io_obj.h"
#include "file_io_lgm.h"
#include "file_io_lgb.h"
#include "file_io_ng.h"
#include "file_io_ug.h"
#include "file_io_dump.h"
#include "lib_grid/lg_base.h"

#ifndef __H__LIB_GRID__FILE_IO__
#define __H__LIB_GRID__FILE_IO__

namespace ug
{
//	methods that load a file by automatically choosing the right method
//	from the filenames suffix should be added here.

////////////////////////////////////////////////////////////////////////
//	LoadGridFromFile
///	Loads a grid from a file. The importer is chosen by the filenames suffix.
bool LoadGridFromFile(Grid& grid, const char* filename, AVector3& aPos = aPosition);

////////////////////////////////////////////////////////////////////////
//	SaveGridToFile
///	Saves a grid from a file. The exporter is chosen by the filenames suffix.
bool SaveGridToFile(Grid& grid, const char* filename, AVector3& aPos = aPosition);

////////////////////////////////////////////////////////////////////////
//	LoadGridFromFile
///	Loads a grid from a file. The importer is chosen by the filenames suffix.
/**
 * Supports a subset-handler.
 * If the file does not support subsets, all elements of highest dimension
 * will automatically be added to subset 0.
 */
bool LoadGridFromFile(Grid& grid, const char* filename,
						SubsetHandler& sh,
						AVector3& aPos = aPosition);

////////////////////////////////////////////////////////////////////////
//	SaveGridToFile
///	Saves a grid from a file. The exporter is chosen by the filenames suffix.
/**
 * Supports a subset-handler.
 */
bool SaveGridToFile(Grid& grid, const char* filename,
					SubsetHandler& sh,
					AVector3& aPos = aPosition);

};

#endif
