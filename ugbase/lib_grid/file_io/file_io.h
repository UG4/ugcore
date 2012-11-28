//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m11 d13

#include "file_io_art.h"
#include "file_io_txt.h"
#include "file_io_tetgen.h"
#include "file_io_obj.h"
#include "file_io_lgm.h"
#include "file_io_lgb.h"
#include "file_io_ng.h"
#include "file_io_ug.h"
#include "file_io_dump.h"
#include "file_io_ncdf.h"
#include "file_io_ugx.h"
#include "file_io_msh.h"
#include "file_io_stl.h"

#ifndef __H__LIB_GRID__FILE_IO__
#define __H__LIB_GRID__FILE_IO__

namespace ug
{

//	predeclaration of surface view (declared in lib_grid/tools/surface_view.h)
class SurfaceView;

////////////////////////////////////////////////////////////////////////////////
///	Loads a grid from a file. Position data is written to the specified attachment.
/**
 * Make sure that the given position attachment is either of type AVector1,
 * AVector2 or AVector3.
 *
 * If the given file can't be found, LoadGridFromFile will looks for it reative
 * to the following additional places:
 * 	- PathProvider::get_current_path()
 * 	- PathProvider::get_path(GRID_PATH)
 * \{
 */
template <class TAPos>
UG_API
bool LoadGridFromFile(Grid& grid, ISubsetHandler& sh,
					  const char* filename, TAPos& aPos);

template <class TAPos>
UG_API
bool LoadGridFromFile(Grid& grid, const char* filename, TAPos& aPos);
/**	\} */

////////////////////////////////////////////////////////////////////////////////
///	Loads a grid from a file. Position data is written to aPosition.
/**
 * If the given file can't be found, LoadGridFromFile will looks for it reative
 * to the following additional places:
 * 	- PathProvider::get_current_path()
 * 	- PathProvider::get_path(GRID_PATH)
 *
 * \{ */
UG_API
bool LoadGridFromFile(Grid& grid, ISubsetHandler& sh, const char* filename);

UG_API
bool LoadGridFromFile(Grid& grid, const char* filename);
/** \} */


////////////////////////////////////////////////////////////////////////////////
///	Saves a grid to a file. Position data is read from the specified attachment.
/**
 * Make sure that the given position attachment is either of type AVector1,
 * AVector2 or AVector3.
 *
 * \{
 */
template <class TAPos>
UG_API
bool SaveGridToFile(Grid& grid, ISubsetHandler& sh,
					const char* filename, TAPos& aPos);

template <class TAPos>
UG_API
bool SaveGridToFile(Grid& grid, const char* filename, TAPos& aPos);
/**	\} */

////////////////////////////////////////////////////////////////////////////////
///	Saves a grid to a file. Position data is read from aPosition.
/** \{ */
UG_API
bool SaveGridToFile(Grid& grid, ISubsetHandler& sh, const char* filename);

UG_API
bool SaveGridToFile(Grid& grid, const char* filename);
/** \} */


///	Saves a grid hierarchy by offsetting levels along the z-axis.
/**	\todo:	Support any type of position attachment*/
bool SaveGridHierarchyTransformed(MultiGrid& mg, ISubsetHandler& sh,
								  const char* filename, number offset);

///	Saves a grid hierarchy by offsetting levels along the z-axis.
/**	\todo:	Support any type of position attachment*/
bool SaveGridHierarchyTransformed(MultiGrid& mg, const char* filename,
									  number offset);

///	Saves the grid-layout of parallel multi-grids
/**	\todo:	Support any type of position attachment*/
bool SaveParallelGridLayout(MultiGrid& mg, const char* filename,
							number offset = 0.1);


///	Saves a grid and assigns elements to subsets based on their surface-view-state.
/**	\todo:	Support any type of position attachment*/
bool SaveSurfaceViewTransformed(MultiGrid& mg, const SurfaceView& sv,
								const char* filename, number offset = 0.1);

};

#endif
