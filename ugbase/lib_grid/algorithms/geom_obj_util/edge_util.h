//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m02 d02

#ifndef __H__LIB_GRID__EDGE_UTIL__
#define __H__LIB_GRID__EDGE_UTIL__

#include "lib_grid/lg_base.h"
#include "face_util.h"

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
//	GetAssociatedFaces
///	writes associated faces of e to facesOut.
/**
 * This method uses ug::Grid::mark.
 *
 * Associated faces of e are written to facesOut.
 * facesOut has to be an array of size maxNumFaces.
 * If there are more then maxNumFaces associated faces, they are not
 * written to facesOut.
 *
 * The method returns the number of total number of associated faces.
 */
int GetAssociatedFaces(Face** facesOut, Grid& grid,
						EdgeBase* e, int maxNumFaces);

////////////////////////////////////////////////////////////////////////
//	CalculateNormal
///	Calculates the normal of the given edge
/**
 * This method indirectly uses ug::Grid::mark.
 *
 * The normal is calculated as the normized sum of associated face normals.
 * If there are no associated faces, the normal is set to (0, 0, 0) and
 * 0 is returned.
 *
 * \param paaNormFACE: An optional parameter that allows to specify an
 *						accessor for precalculated face normals.
 *
 * \returns the number of faces that are associated with the edge.
 */
int CalculateNormal(vector3& vNormOut, Grid& grid, EdgeBase* e,
					Grid::VertexAttachmentAccessor<APosition>& aaPos,
					Grid::FaceAttachmentAccessor<ANormal>* paaNormFACE = NULL);
					
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
 * If bConservative == false then SplitEdge will replace e and its adjacent
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
//	SwapEdge
///	swaps e and thus reconnects its two adjacent triangles.
/**
 *	A swap is only allowed if e has exactly 2 adjacent triangles.
 *	The method erases the old edge and triangles and constructs new ones.
 *	Old elements are passed as parents to the grids creation method.
 *
 *	The swapped edge is returned.
 */
EdgeBase* SwapEdge(Grid& grid, EdgeBase* e);

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


////////////////////////////////////////////////////////////////////////
/**
 * paFaceNormal is ignored in the current implementation.
 * In the moment normals are calculated on the fly and not stored.
 * That means that the normal of each single face is calculated up to
 * four times. This can be improved!
 */
template <class TEdgeIterator>
void MarkCreaseEdges(Grid& grid, ISubsetHandler& sh,
					TEdgeIterator edgesBegin, TEdgeIterator edgesEnd,
					int subsetIndex, number angle,
					APosition& aPos = aPosition,
					ANormal* paFaceNormal = NULL);


template<class TVertexPositionAttachmentAccessor>
typename TVertexPositionAttachmentAccessor::ValueType
CalculateCenter(EdgeBase* e, TVertexPositionAttachmentAccessor& aaPosVRT);

								
/**@}*/ // end of doxygen defgroup command

}//	end of namespace

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	include template-methods implementations
#include "edge_util_impl.hpp"

#endif
