//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m02 d02

#ifndef __H__LIB_GRID__VOLUME_UTIL__
#define __H__LIB_GRID__VOLUME_UTIL__

#include <vector>
#include "lib_grid/lg_base.h"

namespace ug
{
////////////////////////////////////////////////////////////////////////
//	CalculateMinTetrahedronHeight - mstepnie
/// calculates the minimal height of a tetrahedral element
//double CalculateMinTetrahedronHeight(Grid& grid, Volume& v);
number CalculateMinTetrahedronHeight(const vector3& a, const vector3& b, 
									 const vector3& c, const vector3& d);


////////////////////////////////////////////////////////////////////////
//	CalculateMinTetrahedronEdge - mstepnie
/// calculates the shortest edge of a tetrahedral element
number CalculateMaxTetrahedronEdgelength(Grid& grid, Volume& v);


////////////////////////////////////////////////////////////////////////
//	CalculateTetrahedronAspectRatio - mstepnie
/// calculates the aspect ratio of a tetrahedral element
number CalculateTetrahedronAspectRatio(Grid& grid, Volume& v);


////////////////////////////////////////////////////////////////////////
//	CalculateTetrahedronVolume - mstepnie
/// calculates the volume of a tetrahedral element
number CalculateTetrahedronVolume(const vector3& a, const vector3& b,
								  const vector3& c, const vector3& d);
	
}//	end of namespace

#endif
