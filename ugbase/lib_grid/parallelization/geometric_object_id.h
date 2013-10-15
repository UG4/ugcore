// created by Sebastian Reiter
// s.b.reiter@gmail.com
// May 2013

#ifndef __H__LIB_GRID__GEOMETRIC_OBJECT_ID__
#define __H__LIB_GRID__GEOMETRIC_OBJECT_ID__

#include <algorithm>
#include <ostream>
#include "lib_grid/grid/geometric_base_objects.h"
#include "common/util/hash_function.h"

namespace ug{

/// \addtogroup lib_grid_parallelization
/// @{

////////////////////////////////////////////////////////////////////////
///	The global id can be used to uniquely identify distributed objects.
/**	Note that normally no global IDs are associated with geometric objects.
 * However, methods exist, which can assign global ids based on the
 * current layouts.
 */
typedef std::pair<int, size_t>	GeomObjID;

UG_API std::ostream& operator<<(std::ostream& out, const GeomObjID& goId);

UG_API bool operator<(const GeomObjID& gid1, const GeomObjID& gid2);

////////////////////////////////////////////////////////////////////////
///	Can be used to construct a GeomObjID from a proc-rank and a local id.
inline GeomObjID MakeGeomObjID(int procRank, size_t localGeomObjID)
{
	return std::make_pair(procRank, localGeomObjID);
}

////////////////////////////////////////////////////////////////////////
///	An attachment which can store GeomObjIDs
typedef Attachment<GeomObjID>	AGeomObjID;

///	This attachment instance should be used to store global ids
extern AGeomObjID aGeomObjID;

///	generates a hash key for a GeomObjID.
/**	\todo Check distribution quality.*/
template <>
size_t hash_key<GeomObjID>(const GeomObjID& key);

}// end of namespace

#endif
