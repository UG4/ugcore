//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m02 d02

#ifndef __H__LIB_GRID__EDGE_UTIL__
#define __H__LIB_GRID__EDGE_UTIL__

#include "face_util.h"
#include "lib_grid/grid/grid.h"
#include "lib_grid/grid/grid_util.h"
#include "lib_grid/grid_objects/grid_objects.h"
#include "lib_grid/common_attachments.h"
#include "lib_grid/tools/selector_grid.h"
#include "lib_grid/tools/subset_handler_interface.h"

namespace ug
{
/**
 * \brief contains methods to manipulate edges
 * 
 * \defgroup lib_grid_algorithms_edge_util edge util
 * \ingroup lib_grid_algorithms
 * @{
 */

////////////////////////////////////////////////////////////////////////
//	GetEdgeIndex
///	returns the index at which edge e is found in the given object
/**
 * returns -1 if the edge was not found.
 */
UG_API int GetEdgeIndex(Face* f, Edge* e);

////////////////////////////////////////////////////////////////////////
//	GetEdgeIndex
///	returns the index at which edge e is found in the given object
/**
 * returns -1 if the edge was not found.
 */
UG_API int GetEdgeIndex(Volume* vol, Edge* e);


///	returns true if the edge is connected to exactly one surface face.
/**	If the given callback returns true for a given face, then the face is
 * considered to be part of a surface. If exactly one of the faces adjacent to
 * the given edge is part of a surface, then the edge is considered to be
 * a boundary edge and true is returned.
 * Take a look at existing standard callbacks, if you want to use this method.
 */
UG_API 
bool IsBoundaryEdge(Grid& grid, Edge* e,
					Grid::face_traits::callback funcIsSurfFace);


////////////////////////////////////////////////////////////////////////
///	returns whether an edge lies on the boundary of a 2D grid.
/**	An edge is regarded as a boundary edge if it is adjacent
 *	to exactly one face.
 *	if EDGEOPT_STORE_ASSOCIATED_FACES is enabled, the algorithm will be faster.
 */
UG_API bool IsBoundaryEdge2D(Grid& grid, Edge* e);

////////////////////////////////////////////////////////////////////////
///	returns whether an edge lies on the boundary of a 3D grid.
/**	An edge is regarded as a boundary edge in 3d if it is adjacent
 *	to at least one boundary face.
 *	if EDGEOPT_STORE_ASSOCIATED_FACES is enabled, the algorithm will be faster.
 *
 *	Please Note: This algorithm requires the grid option VOLOPT_AUTOGENERATE_FACES.
 *	If it is not enabled, this algorithm will enable it.
 *	\todo This algorithm should work without VOLOPT_AUTOGENERATE_FACES, too.
 */
UG_API bool IsBoundaryEdge3D(Grid& grid, Edge* e);

////////////////////////////////////////////////////////////////////////
///	returns true, if the edge lies on a 2d or 3d boundary
UG_API bool LiesOnBoundary(Grid& grid, Edge* e);

////////////////////////////////////////////////////////////////////////
//	GetAssociatedFaces
///	writes associated faces of e to facesOut.
/**
 * Associated faces of e are written to facesOut.
 * facesOut has to be an array of size maxNumFaces.
 * If there are more then maxNumFaces associated faces, they are not
 * written to facesOut.
 *
 * The method returns the number of total number of associated faces.
 */
UG_API 
int GetAssociatedFaces(Face** facesOut, Grid& grid,
						Edge* e, int maxNumFaces);

////////////////////////////////////////////////////////////////////////
//	NumAssociatedFaces
///	returns the number of associated faces of the given edge
/**
 * This method uses ug::Grid::mark.
 *
 * The method returns the number of total number of associated faces.
 */
UG_API 
int NumAssociatedFaces(Grid& grid, Edge* e);
						
////////////////////////////////////////////////////////////////////////
///	Calculates the squared length of the given edge
/**	The specified accessor has to access a MathVector compatible type
 * in the vertices of the underlying grid.*/
template <class TAAPosVRT>
UG_API 
inline number EdgeLengthSq(const EdgeVertices* e, TAAPosVRT& aaPos);

////////////////////////////////////////////////////////////////////////
///	Calculates the length of the given edge
/**	The specified accessor has to access a MathVector compatible type
 * in the vertices of the underlying grid.*/
template <class TAAPosVRT>
UG_API 
inline number EdgeLength(const EdgeVertices* e, TAAPosVRT& aaPos);

////////////////////////////////////////////////////////////////////////
//	CalculateNormal
///	Calculates the normal of the given edge
/**
 * This method indirectly uses ug::Grid::mark.
 *
 * The normal is calculated as the normized sum of associated face normals.
 * If there are no associated faces, it is assumed, that the edge lies
 * in the xy plane. Its normal will be calculated to the right.
 *
 * \param vNormOut		normal
 * \param grid			Grid
 * \param e				RegularEdge
 * \param aaPos			vertec attachment accessor
 * \param paaNormFACE 	An optional parameter that allows to specify an
 *						accessor for precalculated face normals.
 *
 * \returns the number of faces that are associated with the edge.
 */
UG_API 
int CalculateNormal(vector3& vNormOut, Grid& grid, Edge* e,
					Grid::AttachmentAccessor<Vertex, APosition>& aaPos,
					Grid::AttachmentAccessor<Face, ANormal>* paaNormFACE = NULL);

////////////////////////////////////////////////////////////////////////
//	CalculateNormal
///	Calculates the normal of the given edge
/**
 * This method indirectly uses ug::Grid::mark.
 *
 * The normal is calculated as the normized sum of associated face normals.
 * If there are no associated faces, it is assumed, that the edge lies
 * in the xy plane. Its normal will be calculated to the right.
 *
 * \param vNormOut		normal
 * \param grid			Grid
 * \param e				RegularEdge
 * \param aaPos			vertec attachment accessor
 * \param paaNormFACE: An optional parameter that allows to specify an
 *						accessor for precalculated face normals.
 *
 * \returns the number of faces that are associated with the edge.
 */
UG_API 
int CalculateNormalNoNormalize(vector3& vNormOut, Grid& grid, Edge* e,
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
UG_API 
bool CollapseEdge(Grid& grid, Edge* e, Vertex* newVrt);

////////////////////////////////////////////////////////////////////////
//	EdgeCollapseIsValid
///	Checks if an edge-collapse would invalidate the current topology.
/**
 * returns true if the topology would not be affected by the collapse.
 * returns false if the topology would be affected.
 */
UG_API 
bool EdgeCollapseIsValid(Grid& grid, Edge* e);


////////////////////////////////////////////////////////////////////////
//	SplitEdge
///	inserts new triangles and one new vertex by splitting the specified edge.
/**
 * returns the newly created vertex if everything went right, NULL if not.
 * The vertex that will be created will be of type TVertex.
 * If bConservative == false then SplitEdge will replace e and its adjacent
 * geometry by the newly generated geometry.
 */
template<class TVertex>
UG_API 
TVertex* SplitEdge(Grid& grid, Edge* e, bool bConservative = false);

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
UG_API 
TVertex* SplitEdge(Grid& destGrid, Grid& srcGrid, Edge* e,
						AVertex* paAssociatedVertices = NULL,
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
UG_API 
Edge* SwapEdge(Grid& grid, Edge* e);

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
UG_API 
bool CreateEdgeSplitGeometry(Grid& destGrid, Grid& srcGrid, Edge* e,
							 Vertex* newVertex,
							 AVertex* paAssociatedVertices = NULL);


////////////////////////////////////////////////////////////////////////
///	Calculates the center of an edge
template<class TVertexPositionAttachmentAccessor>
UG_API 
typename TVertexPositionAttachmentAccessor::ValueType
CalculateCenter(const Edge* e, TVertexPositionAttachmentAccessor& aaPosVRT);


////////////////////////////////////////////////////////////////////////
///	returns the weighted center of the vertices of the given edge
/** TAAWeightVRT has to be an attachment to the vertices of the grid in which
 * e is contained, with ValueType number (or compatible).
 */
template<class TAAPosVRT, class TAAWeightVRT>
UG_API
typename TAAPosVRT::ValueType
CalculateCenter(const EdgeVertices* e, TAAPosVRT& aaPos, TAAWeightVRT& aaWeight);


////////////////////////////////////////////////////////////////////////
///	refines all edges in sel which cut the given plane.
/**	New vertices are inserted on the plane.
 *	When the method is done, sel will contain all refined elements.*/
UG_API 
bool CutEdgesWithPlane(Selector& sel, const vector3& p, const vector3& n,
						APosition& aPos = aPosition);

////////////////////////////////////////////////////////////////////////
//	FixOrientation
///	creates uniform orientation of neighboured edges.
/** This algorithm uses Grid::mark
 *
 * swaps orientation of edges so that all neighboured
 * edges share the same.
 *
 * Value type of TEdgeIterator has to be compatible with Edge*.
 *
 * Note that all edges between edgesBegin and edgesEnd have to be members
 * of the specified grid.
 *
 * The orientation can only be successfully fixed, if vertices between
 * the given edges share at most 2 edges between edgesBegin and edgesEnd.
 */
template <class TEdgeIterator>
UG_API 
void FixEdgeOrientation(Grid& grid, TEdgeIterator edgesBegin,
						TEdgeIterator edgesEnd);

////////////////////////////////////////////////////////////////////////////////
///	Orientates boundary edges in the given edge set to the orientation of associated faces.
/**	The orientation of boundary edges will match the orientation of associated
 * faces, after this algorithm terminates.
 */
template <class TEdgeIterator>
UG_API
void AdjustEdgeOrientationToFaceOrientation(Grid& grid, TEdgeIterator edgesBegin,
						   	   	   	   	    TEdgeIterator edgesEnd);

////////////////////////////////////////////////////////////////////////////////
///	Returns the shortest edge in a list of edges
/**	TEdgeIterator has to point to values of type Edge* and has to be an
 * stl-iterator compatible iterator. TAAPosVRT has to be an attachment accessor
 * which provides a position value (vector1, vector2, vector3, ...) for the
 * vertices of the edge.
 * If the specified list is empty, NULL is returned.
 * If multiple shortest edges exist, the first one is returned.
 */
template <class TEdgeIterator, class TAAPosVRT>
UG_API 
Edge* FindShortestEdge(TEdgeIterator edgesBegin, TEdgeIterator edgesEnd,
							TAAPosVRT& aaPos);

//////////////////////////////////////////////////////////////////////////////////
/////	Removes edges that connect the same two vertices as another edge.
///**	THIS ALGORITHM USES Grid::mark*/
//template <class TEdgeIterator>
//UG_API 
//void RemoveDoubleEdges(Grid& grid, TEdgeIterator edgesBegin, TEdgeIterator edgesEnd);

///	Transforms the given edge-set so that the sum of the length the edges is minimized.
/**
 * Currently only works for 3d position attachments.
 *
 * \todo	add support for 2d position attachments.
 */
template <class EdgeIterator, class TAAPos>
UG_API 
void MinimizeEdgeLength_SwapsOnly(Grid& grid, EdgeIterator edgesBegin,
								  EdgeIterator edgesEnd, TAAPos& aaPos);

////////////////////////////////////////////////////////////////////////
///	Returns true if the given point lies on the given edge.
/**	\note	The method only works properly, if the point and the edge are located
 * 			on a line parallel to the x-axis.
 */
template <class vector_t, class TAAPos>
UG_API bool
ContainsPoint(const EdgeVertices* e, const vector_t& p, TAAPos aaPos);

/// @} // end of doxygen defgroup command

}//	end of namespace

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	include template-methods implementations
#include "edge_util_impl.hpp"

#endif
