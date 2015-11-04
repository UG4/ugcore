#ifndef __H__UG_file_io_tikz
#define __H__UG_file_io_tikz

#include "lib_grid/grid/grid.h"
#include "lib_grid/tools/subset_handler_grid.h"
#include "lib_grid/common_attachments.h"

namespace ug{

struct TikzExportDesc{
	TikzExportDesc();

	number	vrtRadius;
	number	vrtRimWidth;
	vector3	vrtColor;
	number	edgeWidth;
	vector3	edgeColor;
	number	faceRimWidth;
	vector3	faceColor;
	number	smallestVal;
};

bool ExportGridToTIKZ(Grid& grid, const char* filename, const ISubsetHandler* psh,
					  APosition aPos, TikzExportDesc desc = TikzExportDesc());

}//	end of namespace

#endif	//__H__file_io_tikz
