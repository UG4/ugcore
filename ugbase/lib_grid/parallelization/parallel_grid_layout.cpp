// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 08.03.2011 (m,d,y)
 
#include "parallel_grid_layout.h"

namespace ug{
AGeomObjID 	aGeomObjID("globalID", false);

template <>
unsigned long hash_key<GeomObjID>(const GeomObjID& key)
{
//	of course this hash does not completly avoid collisions.
//	One should check whether the chosen key is fine.
	return (unsigned long)(99971 * key.first + key.second * key.second);
}

}// end of namespace
