//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m11 d19

#include <list>
#include <vector>
#include <cassert>
#include <map>
#include <stack>
#include "lib_grid/lg_base.h"
#include "common/ug_config.h"

#ifndef __H__UG4__LIB_GRID__TRIANGLE_FILL_SWEEP_LINE__
#define __H__UG4__LIB_GRID__TRIANGLE_FILL_SWEEP_LINE__

namespace ug
{

///	Performs triangulation of a polygon and resoves additional inner edges.
/**	This algorithm triangulates a set of edges.
 *	It is important, that a closed outer boundary exists.
 *	All edges are resolved in the final geometry.*/
UG_API
bool TriangleFill_SweepLine(std::vector<int>& facesOut,
							const std::vector<vector2>& srcVrts,
							/*const */std::vector<int>& srcEdges);

///	Performs triangulation of a 3d polygon and resolves inner edges.
/**	The polygon should lie in a 2d hyperplane.*/
UG_API
bool TriangleFill_SweepLine(std::vector<int>& facesOut,
							const std::vector<vector3>& srcVrts,
							/*const */std::vector<int>& srcEdges);

///	Performs triangulation of a polygon and resoves additional inner edges.
/**
 *	This algortighm uses Grid::mark.
 *
 *	The polygon should lie in a 2d hyperplane.
 */
template <class TIterator>
bool TriangleFill_SweepLine(Grid& grid, TIterator edgesBegin,
							TIterator edgesEnd, APosition& aPosVRT,
							AInt& aIntVRT,
							SubsetHandler* pSH = NULL,
							int newSubsetIndex = -1);

}//	namespace ug

////////////////////////////////
//	include implementation
#include "triangle_fill_sweep_line_impl.hpp"

#endif
