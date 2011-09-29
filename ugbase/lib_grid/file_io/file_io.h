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
#include "lib_grid/lg_base.h"

#ifndef __H__LIB_GRID__FILE_IO__
#define __H__LIB_GRID__FILE_IO__

namespace ug
{

////////////////////////////////////////////////////////////////////////////////
///	Loads a grid from a file. Position data is written to the specified attachment.
/**
 * Make sure that the given position attachment is either of type AVector1,
 * AVector2 or AVector3.
 *
 * \{
 */
template <class TAPos>
bool LoadGridFromFile(Grid& grid, ISubsetHandler& sh,
					  const char* filename, TAPos& aPos);

template <class TAPos>
bool LoadGridFromFile(Grid& grid, const char* filename, TAPos& aPos);
/**	\} */

////////////////////////////////////////////////////////////////////////////////
///	Loads a grid from a file. Position data is written to aPosition.
/** \{ */
bool LoadGridFromFile(Grid& grid, ISubsetHandler& sh, const char* filename);

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
bool SaveGridToFile(Grid& grid, SubsetHandler& sh,
					const char* filename, TAPos& aPos);

template <class TAPos>
bool SaveGridToFile(Grid& grid, const char* filename, TAPos& aPos);
/**	\} */

////////////////////////////////////////////////////////////////////////////////
///	Saves a grid to a file. Position data is read from aPosition.
/** \{ */
bool SaveGridToFile(Grid& grid, SubsetHandler& sh, const char* filename);

bool SaveGridToFile(Grid& grid, const char* filename);
/** \} */

};

#endif
