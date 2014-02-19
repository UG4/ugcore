// created by Sebastian Reiter
// s.b.reiter@gmail.com
// Jun 21, 2013 (d,m,y)

#ifndef __H__UG__resolve_intersections__
#define __H__UG__resolve_intersections__

#include "lib_grid/grid/grid.h"

namespace ug{

template <class TAAPosVRT>
Vertex* ResolveVertexEdgeIntersection(Grid& grid, Vertex* v,
										   EdgeBase* e, TAAPosVRT& aaPos,
										   number snapThreshold);

/**
 * No support for volumes in the current version.
 * \todo Instead of manually refining the face, an external function SplitFace
 *		 should be used, which can take care of volumes, too.
 */
template <class TAAPosVRT>
bool ResolveVertexFaceIntersection(Grid& grid, Vertex* v,
								   Face* f, TAAPosVRT& aaPos,
								   number snapThreshold);

/**
 * This method does not resolve intersections between close, parallel edges or
 * between degenerate edges. You can treat such cases with
 * ReolveVertexEdgeIntersection.
 */
template <class TAAPosVRT>
Vertex* ResolveEdgeEdgeIntersection(Grid& grid, EdgeBase* e1, EdgeBase* e2,
										TAAPosVRT& aaPos, number snapThreshold);

/**
 * No support for volumes in the current version.
 * \todo Instead of manually refining the face, an external function SplitFace
 *		 should be used, which can take care of volume, too.
 */
template <class TAAPosVRT>
bool ResolveEdgeFaceIntersection(Grid& grid, EdgeBase* e, Face* f,
								 TAAPosVRT& aaPos, number snapThreshold);

/**
 *	Projects vertices in elems onto close edges in elems.
 *	Though this method can be used to remove degenerated triangles,
 *	it is not guaranteed, that no degenerated triangles will remain
 *	(indeed, new degenerated triangles may be introduced).
 */
template <class TAAPosVRT>
bool ProjectVerticesToCloseEdges(Grid& grid,
								 GridObjectCollection elems,
								 TAAPosVRT& aaPos,
								 number snapThreshold);

/**
 *	Projects vertices in elems onto close faces in elems.
 */
template <class TAAPosVRT>
bool ProjectVerticesToCloseFaces(Grid& grid,
								 GridObjectCollection elems,
								 TAAPosVRT& aaPos,
								 number snapThreshold);

/**THIS METHOD USES Grid::mark.
 * Intersects all edges in elems which are closer to each other
 * than snapThreshold.*/
template <class TAAPosVRT>
bool IntersectCloseEdges(Grid& grid,
						 GridObjectCollection elems,
						 TAAPosVRT& aaPos,
						 number snapThreshold);


///	returns the index of the first vertex closer to p than snapThreshold.
/**	returns -1 if nothing was found.*/
template <class TAAPosVRT>
int FindCloseVertexInArray(std::vector<Vertex*>& array,
							const typename TAAPosVRT::ValueType& p,
							TAAPosVRT& aaPos, number snapThreshold);

////////////////////////////////////////////////////////////////////////
/**	This method uses Grid::mark
 */
template <class TAAPosVRT>
bool ResolveGridIntersections(Grid& grid, TriangleIterator trisBegin,
							  TriangleIterator trisEnd, number snapThreshold,
							  TAAPosVRT& aaPos);

}// end of namespace

////////////////////////////////////////
//	include implementation
#include "resolve_intersections_impl.hpp"
#endif
