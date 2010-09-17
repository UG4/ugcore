//	created by Sebastian Reiter, Martin Stepniewski
//	s.b.reiter@googlemail.com
//	y09 m05 d27

#ifndef __H__LIBGRID__FILE_IO_UG__
#define __H__LIBGRID__FILE_IO_UG__

#include "lib_grid/lg_base.h"

namespace ug
{
////////////////////////////////////////////////////////////////////////
/**
 * shFace has to contain the interface-faces that separate the different
 * volume-subsets (specified in shVolume).
 * 
 * lgmName, problemName and convex correlate to the parameters that appear
 * at the beginning of each lgm-file.
 */
bool ExportGridToUG(const Grid& g, const SubsetHandler& shFace,
					const SubsetHandler& shVolume,
					const char* fileNamePrefix, const char* lgmName,
					const char* problemName, int convex);

////////////////////////////////////////////////////////////////////////
///	exports a grid to a 2d lgm / ng file combination.
/**	Please note that edge-subsets are ignored.
 *	Only x and y coordinates are written.
 */
bool ExportGridToUG_2D(Grid& grid, const char* fileName, const char* lgmName,
					   const char* problemName, int convex,
					   SubsetHandler* psh = NULL);
}//	end of namespace

#endif
