// created by Sebastian Reiter
// s.b.reiter@gmail.com

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
					ANormal* paFaceNormal = NULL);

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
