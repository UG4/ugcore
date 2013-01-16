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



///	Returns the center of the given element (SLOW - for debugging only!)
/**	Caution: This method is pretty slow and should only be used for debugging
 * purposes. It only works if one of the standard position attachments
 * aPosition, aPosition2 or aPosition1 is attached to the vertices of g.
 *
 * The method returns the center of the specified element in a vector3 structure,
 * regardless of the actual dimension of the position attachment. Unused
 * coordinates are set to 0.
 *
 * TElem has to be derived from GeometricObject
 *
 * \{
 */
template <class TElem>
vector3 GetGeometricObjectCenter(Grid& g, TElem* elem);

/** \} */


/// Writes level and center of each object together with a custom value to a file
template <class TElem, class TAValue>
void WriteDebugValuesToFile(const char* filename, Grid& grid,
							TAValue& aVal, bool levelWise);


///	returns the index of the given element in the given grid.
/**	Runtime O(n). Returns -1 if the element could not be found.
 * \note: The all elements of the same base-type are considered during counting.*/
template <class TElem>
int GetGeometricObjectIndex(Grid& g, TElem* elem);


///	checks whether all constraining and constrained objects are correctly connected.
bool CheckHangingNodeConsistency(Grid& g);

///	checks whether a multigrid is a valid haging vertex grid.
/**	This algorithm e.g. checks, whether adaptivity is represented by
 * constrained / constraining objects.
 *
 * Calls CheckHangingVertexConsistency(Grid&)
 */
bool CheckHangingNodeConsistency(MultiGrid& mg);

/**@}*/ // end of doxygen defgroup command

}//	end of namespace

////////////////////////////////////////
//	include implementation
#include "debug_util_impl.hpp"

#endif
