/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 ��7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., N��gel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef __H__UG__resolve_intersections__
#define __H__UG__resolve_intersections__

#include "lib_grid/grid/grid.h"

namespace ug{

template <class TAAPosVRT>
Vertex* ResolveVertexEdgeIntersection(Grid& grid, Vertex* v,
										   Edge* e, TAAPosVRT& aaPos,
										   number snapThreshold);

/**
 * No support for volumes in the current version.
 * \todo Instead of manually refining the face, an external function SplitFace
 *		 should be used, which can take care of volumes, too.
 */
template <class TAAPosVRT>
bool ResolveVertexFaceIntersection(Grid& grid, Vertex* v,
								   Face* f, TAAPosVRT& aaPos,
								   number snapThreshold,
								   std::vector<Face*>* pNewFacesOut);

/**
 * This method does not resolve intersections between close, parallel edges or
 * between degenerate edges. You can treat such cases with
 * ReolveVertexEdgeIntersection.
 */
template <class TAAPosVRT>
Vertex* ResolveEdgeEdgeIntersection(Grid& grid, Edge* e1, Edge* e2,
										TAAPosVRT& aaPos, number snapThreshold);

/**
 * No support for volumes in the current version.
 * \todo Instead of manually refining the face, an external function SplitFace
 *		 should be used, which can take care of volume, too.
 */
template <class TAAPosVRT>
bool ResolveEdgeFaceIntersection(Grid& grid, Edge* e, Face* f,
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
template <class TObjectCollection, class TAPos>
bool ProjectVerticesToCloseFaces(Grid& grid,
								 TObjectCollection& elems,
								 TAPos& aPos,
								 number snapThreshold);

/**THIS METHOD USES Grid::mark.
 * Intersects all edges in elems which are closer to each other
 * than snapThreshold.*/
template <class TObjectCollection, class TAAPosVRT>
bool IntersectCloseEdges(Grid& grid,
						 TObjectCollection& elems,
						 TAAPosVRT& aaPos,
						 number snapThreshold);


///	returns the index of the first vertex closer to p than snapThreshold.
/**	returns -1 if nothing was found.*/
template <class TAAPosVRT>
int FindCloseVertexInArray(std::vector<Vertex*>& array,
							const typename TAAPosVRT::ValueType& p,
							TAAPosVRT& aaPos, number snapThreshold);

/// returns the index of the closest vertex to p if closer than snapThreshold.
template <class TAAPosVRT>
int FindClosestVertexInArray(std::vector<Vertex*>& array,
							const typename TAAPosVRT::ValueType& p,
							TAAPosVRT& aaPos, number snapThreshold);

////////////////////////////////////////////////////////////////////////
/**	This method uses Grid::mark
 */
template <class TAPos>
bool ResolveTriangleIntersections(Grid& grid, TriangleIterator trisBegin,
							  TriangleIterator trisEnd, number snapThreshold,
							  TAPos& aPos);

}// end of namespace

////////////////////////////////////////
//	include implementation
#include "resolve_intersections_impl.hpp"
#endif
