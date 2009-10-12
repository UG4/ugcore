//	created by Sebastian Reiter, Martin Stepniewski
//	s.b.reiter@googlemail.com
//	y09 m05 d27

#ifndef __H__LIBGRID__FILE_IO_UG__
#define __H__LIBGRID__FILE_IO_UG__

#include "lib_grid/lg_base.h"

namespace ug
{

/**
 * sh serves as surface handler for both the faces and volumes of the grid.
 * Please make sure that the faces are assigned to the subsets from
 * 0 to numFaceSubsets and the volumes are assigned from 0 to numVolumeSubsets.
 *
 * lgmName, problemName and convex correlate to the parameters that appear
 * at the beginning of each lgm-file.
 */
bool ExportGridToUG(Grid& grid, SubsetHandler& sh, const char* fileNamePrefix,
					const char* lgmName, const char* problemName, int convex);

}//	end of namespace

#endif
