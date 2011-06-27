// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 13.01.2011 (m,d,y)

#ifndef __H__UG__LIB_GIRD__DEBUG_UTIL__
#define __H__UG__LIB_GIRD__DEBUG_UTIL__

#include "lib_grid/lg_base.h"

namespace ug
{

/**
 * Methods that log informations on grids, subset-handlers, ...
 * \defgroup lib_grid_algorithms_log_util log util
 * \ingroup lib_grid_algorithms
 * @{
 */

///	prints how many elements of each type exist in the goc.
void PrintElementNumbers(const GeometricObjectCollection& goc);

///	prints how many elements of each type exist in the grid.
void PrintGridElementNumbers(Grid& grid);

///	prints how many elements of each type exist in the multi grid.
void PrintGridElementNumbers(MultiGrid& mg);

///	prints how many elements of each type exist in the subset handler.
void PrintGridElementNumbers(GridSubsetHandler& sh);


///	prints information on all attachments of the specified grid
void PrintAttachmentInfo(Grid& grid);

/**@}*/ // end of doxygen defgroup command


}//	end of namespace

#endif
