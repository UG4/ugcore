/*
 * Copyright (c) 2014-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
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
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef __H__UG_mark_util
#define __H__UG_mark_util

namespace ug{

////////////////////////////////////////////////////////////////////////
/**
 * paFaceNormal is ignored in the current implementation.
 * In the moment normals are calculated on the fly and not stored.
 * That means that the normal of each single face is calculated up to
 * four times. This can be improved!
 */
template <class TEdgeIterator>
UG_API 
void MarkCreaseEdges(Grid& grid, ISubsetHandler& sh,
					TEdgeIterator edgesBegin, TEdgeIterator edgesEnd,
					int subsetIndex, number angle,
					APosition& aPos = aPosition,
					ANormal* paFaceNormal = nullptr);

////////////////////////////////////////////////////////////////////////
/**	Assigns vertices between vrtsBegin and vrtsEnd to the specified subsetIndex
 * if they are adjacent to more than 2 path edges or to exactly 1 path edge or.
 * If a vertex is adjacent to exactly 2 path edges, it will be assigned if the
 * angle between those edges is smaller than the given threshold-angle.*/
template <class TVertexIterator, class TAPosition>
UG_API 
void MarkCorners(Grid& grid, ISubsetHandler& sh,
					TVertexIterator vrtsBegin, TVertexIterator vrtsEnd,
					Grid::edge_traits::callback cbPathEdge,
					int subsetIndex, number angle,
					TAPosition& aPos);


}//	end of namespace

////////////////////////////////////////
//	include implementation
#include "mark_util_impl.h"

#endif	//__H__mark_util
