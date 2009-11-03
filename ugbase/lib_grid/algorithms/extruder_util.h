//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m07 d21

#ifndef __H__LIB_GRID__EXTRUDER_UTIL__
#define __H__LIB_GRID__EXTRUDER_UTIL__

#include "lib_grid/lg_base.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
template <class TIterator>
bool repeated_vertex_extrusion(Grid& grid,
							TIterator iterBegin, TIterator iterEnd,
							int numSteps, const vector3& stepDir);

////////////////////////////////////////////////////////////////////////
template <class TIterator>
bool repeated_edge_extrusion(Grid& grid,
							TIterator iterBegin, TIterator iterEnd,
							int numSteps, const vector3& stepDir);

////////////////////////////////////////////////////////////////////////
template <class TIterator>
bool repeated_face_extrusion(Grid& grid,
							TIterator iterBegin, TIterator iterEnd,
							int numSteps, const vector3& stepDir);

}//	end of namespace

////////////////////////////////////////////////
//	include implementation of template methods.
#include "extruder_util_impl.hpp"

#endif
