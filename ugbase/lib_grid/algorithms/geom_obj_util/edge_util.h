//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m02 d02

#ifndef __H__LIB_GRID__EDGE_UTIL__
#define __H__LIB_GRID__EDGE_UTIL__

#include "lib_grid/lg_base.h"

namespace ug
{
/** \defgroup edgeUtil Edge Util
 * @{
 */

////////////////////////////////////////////////////////////////////////
//	GetEdgeIndex
///	returns the index at which edge e is found in the given object
/**
 * returns -1 if the edge was not found.
 */
int GetEdgeIndex(Face* f, EdgeBase* e);

////////////////////////////////////////////////////////////////////////
//	GetEdgeIndex
///	returns the index at which edge e is found in the given object
/**
 * returns -1 if the edge was not found.
 */
int GetEdgeIndex(Volume* vol, EdgeBase* e);

////////////////////////////////////////////////////////////////////////
///	returns whether an edge lies on the boundary of a 2D grid.
/**	An edge is regarded as a boundary edge if it is adjacent
 *	to exactly one face.
 *	if EDGEOPT_STORE_ASSOCIATED_FACES is enabled, the algorithm will be faster.
 */
bool IsBoundaryEdge2D(Grid& grid, EdgeBase* e);

////////////////////////////////////////////////////////////////////////
//	CollapseEdge
///	Collapses the specified edge performs local grid restructuring.
/**
 * The edge e will be replaced by newVrt.
 * Before calling this method you should check if an edge-collapse
 * won't destroy the topology of your grid. You can do this by calling
 * EdgeCollapseIsValid.
 * During an edge-collapse all adjacent faces will be deleted or
 * replaced by new ones. Same for volumes. Several edges will be
 * deleted as well.
 * if eraseUnusedVertices is set to true, vertices of the collapsed edge
 * that are not used by other objects after the collapse will be removed.
 */
bool CollapseEdge(Grid& grid, EdgeBase* e, VertexBase* newVrt);

////////////////////////////////////////////////////////////////////////
//	EdgeCollapseIsValid
///	Checks if an edge-collapse would invalidate the current topology.
/**
 * returns true if the topology would not be affected by the collapse.
 * returns false if the topology would be affected.
 */
bool EdgeCollapseIsValid(Grid& grid, EdgeBase* e);


////////////////////////////////////////////////////////////////////////
//	SplitEdge
///	inserts new triangles and one new vertex by splitting the specified edge.
/**
 * returns the newly created vertex if everything went right, NULL if not.
 * The vertex that will be created will be of type TVertex.
 * If bConservative == true then SplitEdge will replace e and its adjacent
 * geometry by the newly generated geometry.
 */
template<class TVertex>
TVertex* SplitEdge(Grid& grid, EdgeBase* e, bool bConservative = false);

////////////////////////////////////////////////////////////////////////
//	SplitEdge
///	inserts new triangles and one new vertex by splitting the specified edge.
/**
 * returns the newly created vertex.
 * The vertex that will be created will be of type TVertex.
 * The new vertex and triangles are copied to destGrid.
 * e has to be a member of srcGrid.
 * If bConservative == true then SplitEdge will replace e and its adjacent
 * geometry by the newly generated geometry.
 * paAssociatedVertices has to be specified if destGrid and srcGrid do not match.
 * If destGrid and srcGrid do match, paAssociatedVertices may be specified optionally.
 * paAssociatedVertices has to be a vertex-attachment of srcGrid, that stores for each
 * vertex in srcGrid the associated vertex of destGrid. NULL indicates that
 * no associated vertex exists in destGrid. New ones will be automatically
 * constructed in this case.
 */
template<class TVertex>
TVertex* SplitEdge(Grid& destGrid, Grid& srcGrid, EdgeBase* e,
						AVertexBase* paAssociatedVertices = NULL,
						bool bConservative = false);

////////////////////////////////////////////////////////////////////////
//	CreateEdgeSplitGeometry
///	given an edge and a vertex (the split-vertex) this method constructs the split-geometry.
/**
 * The new triangles are copied to destGrid.
 * e has to be a member of srcGrid.
 * The old edge (e) will not be deleted.
 * paAssociatedVertices has to be specified if destGrid and srcGrid do not match.
 * If destGrid and srcGrid do match, paAssociatedVertices may be specified optionally.
 * paAssociatedVertices has to be a vertex-attachment of srcGrid, that stores for each
 * vertex in srcGrid the associated vertex of destGrid. NULL indicates that
 * no associated vertex exists in destGrid. New ones will be automatically
 * constructed in this case by cloning the associated ones in srcGrid.
 */
bool CreateEdgeSplitGeometry(Grid& destGrid, Grid& srcGrid, EdgeBase* e, VertexBase* newVertex, AVertexBase* paAssociatedVertices = NULL);


/**@}*/ // end of doxygen defgroup command

}//	end of namespace

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	include template-methods implementations
#include "edge_util_impl.hpp"

#endif
