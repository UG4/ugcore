//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m07 d21

#ifndef __H__LIB_GRID__EXTRUDER_UTIL_IMPL__
#define __H__LIB_GRID__EXTRUDER_UTIL_IMPL__

#include <vector>
#include "extruder_util.h"
#include "extrusion/extrusion.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
template <class TIterator>
bool repeated_vertex_extrusion(Grid& grid,
							TIterator iterBegin, TIterator iterEnd,
							int numSteps, const vector3& stepDir)
{
	std::vector<VertexBase*> vElems;
	vElems.assign(iterBegin, iterEnd);

	for(int i = 0; i < numSteps; ++i)
		Extrude(grid, &vElems, NULL, NULL, stepDir);

	return true;
}

////////////////////////////////////////////////////////////////////////
template <class TIterator>
bool repeated_edge_extrusion(Grid& grid,
							TIterator iterBegin, TIterator iterEnd,
							int numSteps, const vector3& stepDir)
{
	std::vector<EdgeBase*> vElems;
	vElems.assign(iterBegin, iterEnd);

	for(int i = 0; i < numSteps; ++i)
		Extrude(grid, NULL, &vElems, NULL, stepDir);

	return true;
}

////////////////////////////////////////////////////////////////////////
template <class TIterator>
bool repeated_face_extrusion(Grid& grid,
							TIterator iterBegin, TIterator iterEnd,
							int numSteps, const vector3& stepDir)
{
	std::vector<Face*> vElems;
	vElems.assign(iterBegin, iterEnd);

	for(int i = 0; i < numSteps; ++i)
		Extrude(grid, NULL, NULL, &vElems, stepDir);

	return true;
}

}//	end of namespace

#endif
