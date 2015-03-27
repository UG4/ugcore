//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m02 d02

#ifndef __H__LIB_GRID__VERTEX_UTIL__
#define __H__LIB_GRID__VERTEX_UTIL__

#include <vector>
#include "lib_grid/grid/grid.h"
#include "lib_grid/grid_objects/grid_objects.h"
#include "lib_grid/common_attachments.h"
#include "lib_grid/tools/subset_handler_grid.h"
#include "lib_grid/callbacks/basic_callbacks.h"
#include "common/math/ugmath.h"

namespace ug
{
/**
 * \brief contains methods to manipulate vertices
 *
 * \defgroup lib_grid_algorithms_vertex_util vertex util
 * \ingroup lib_grid_algorithms
 * @{
 */


////////////////////////////////////////////////////////////////////////
//	GetVertexIndex
///	returns the index at which vertex v is found in the given edge
/**
 * returns -1 if the vertex was not found.
 */
UG_API int GetVertexIndex(EdgeVertices* e, Vertex* v);

////////////////////////////////////////////////////////////////////////
//	GetVertexIndex
///	returns the index at which vertex v is found in the given face
/**
 * returns -1 if the vertex was not found.
 */
UG_API int GetVertexIndex(FaceVertices* f, Vertex* v);

////////////////////////////////////////////////////////////////////////
//	GetVertexIndex
///	returns the index at which vertex v is found in the given volume
/**
 * returns -1 if the vertex was not found.
 */
UG_API int GetVertexIndex(VolumeVertices* vol, Vertex* v);

////////////////////////////////////////////////////////////////////////
//	GetConnectedVertex
///	returns the vertex that is connected to v via e.
/**
 * returns NULL if v is not contained in e.
 */
UG_API Vertex* GetConnectedVertex(Edge* e, Vertex* v);

////////////////////////////////////////////////////////////////////////
//	GetConnectedVertex
///	returns the index of the first vertex that is contained in f and is not contained in e.
UG_API Vertex* GetConnectedVertex(EdgeVertices* e, Face* f);

////////////////////////////////////////////////////////////////////////
//	GetConnectedVertexIndex
///	returns the index of the first vertex that is contained in the specified face and is not contained in the given edge.
UG_API int GetConnectedVertexIndex(Face* f, const EdgeDescriptor& ed);

////////////////////////////////////////////////////////////////////////
///	returns the edge in the triangle tri, which does not contain vrt.
/**	Make sure that tri is a triangle!*/
UG_API Edge* GetConnectedEdge(Grid& g, Vertex* vrt, Face* tri);

////////////////////////////////////////////////////////////////////////
///	returns the number of associated edges of the given vertex
UG_API
int NumAssociatedEdges(Grid& grid, Vertex* v);

////////////////////////////////////////////////////////////////////////
///	returns the number of associated faces of the given vertex
UG_API
int NumAssociatedFaces(Grid& grid, Vertex* v);

////////////////////////////////////////////////////////////////////////
//	CollectSurfaceNeighborsSorted
///	Collects neighbor-vertices in either clockwise or counter clockwise order.
/**	Please note: This algorithm uses Grid::mark.
 *
 *	This method will only work if the triangles connected to the given
 *	vertex are homeomorphic to the unit-disc.
 *
 *	Current implementation requires FACEOPT_AUTOGENERATE_EDGES (could be avoided).
 */
UG_API 
bool CollectSurfaceNeighborsSorted(std::vector<Vertex*>& vNeighborsOut,
								   Grid& grid, Vertex* v);
									
////////////////////////////////////////////////////////////////////////
///	Returns the squared distance between two vertices
template <class TAAPos>
UG_API 
number VertexDistanceSq(Vertex* v0, Vertex* v1, TAAPos& aaPos);

////////////////////////////////////////////////////////////////////////
///	Returns the distance between two vertices
template <class TAAPos>
UG_API 
number VertexDistance(Vertex* v0, Vertex* v1, TAAPos& aaPos);

////////////////////////////////////////////////////////////////////////
//	FindVertexByCoordiante
///	returns the vertex that is the closest to the given coordinate
/**
 * returns NULL if no vertex was found (if iterBegin == iterEnd).
 */
UG_API 
Vertex* FindVertexByCoordiante(vector3& coord, VertexIterator iterBegin,
									VertexIterator iterEnd,
									Grid::VertexAttachmentAccessor<APosition>& aaPos);

////////////////////////////////////////////////////////////////////////
///	calculates the normal of a vertex using associated faces
/**
 * TAAPosVRT has to be an attachment accessor for the vector3 type that
 * works on the vertices in grid.
 */
template <class TAAPosVRT>
UG_API 
void CalculateVertexNormal(vector3& nOut, Grid& grid, Vertex* vrt,
						   TAAPosVRT& aaPos);

////////////////////////////////////////////////////////////////////////
///	calculates the normal of a boundary vertex using associated faces
/**
 * TAAPosVRT has to be an attachment accessor for the vector2 or vector3 type
 * that works on the vertices in grid.
 *
 * Note that this method makes only sense if called for a boundary vertex,
 * which is connected to some volume elements.
 *
 * The returned normal is normalized and will point outwards of the area
 * defined by the associated volumes.
 *
 * Note that this method assumes, that all faces do lie in the x-y-plane.
 */
template <class TAAPosVRT>
UG_API 
void CalculateBoundaryVertexNormal2D(typename TAAPosVRT::ValueType& nOut,
									 Grid& grid, Vertex* vrt,
						   	   	     TAAPosVRT& aaPos);

////////////////////////////////////////////////////////////////////////
///	calculates the normal of a boundary vertex using associated volumes
/**
 * TAAPosVRT has to be an attachment accessor for the vector3 type that
 * works on the vertices in grid.
 *
 * Note that this method makes only sense if called for a boundary vertex,
 * which is connected to some volume elements.
 *
 * The returned normal is normalized and will point outwards of the area
 * defined by the associated volumes.
 */
template <class TAAPosVRT>
UG_API 
void CalculateBoundaryVertexNormal3D(vector3& nOut, Grid& grid, Vertex* vrt,
						   	   	     TAAPosVRT& aaPos);

////////////////////////////////////////////////////////////////////////
//	CalculateVertexNormals
///	calculates the normals of all vertices in grid and stores them in aNorm.
/**
 * aPos has to be attached to grid.
 * If some attachments were not attached correctly, the method returns false.
 * \{
 */
UG_API 
bool CalculateVertexNormals(Grid& grid, APosition& aPos, ANormal& aNorm);

UG_API 
bool CalculateVertexNormals(Grid& grid,
							Grid::AttachmentAccessor<Vertex, APosition>& aaPos,
							Grid::AttachmentAccessor<Vertex, ANormal>& aaNorm);

/** \} */

////////////////////////////////////////////////////////////////////////
//	CalculateBoundingBox
/// calculates the BoundingBox

template <class TVrtIter, class TAPosition>
UG_API 
void
CalculateBoundingBox(typename TAPosition::ValueType& vMinOut,
					 typename TAPosition::ValueType& vMaxOut,
					 TVrtIter vrtsBegin, TVrtIter vrtsEnd,
					 Grid::AttachmentAccessor<Vertex, TAPosition>& aaPos);

////////////////////////////////////////////////////////////////////////
//	CalculateCenter
/// calculates the center of a set of vertices.
/**	The difference to CalculateBarycenter is that this method
 * returns the center of the bounding box which contains the
 * given set of vertices.*/
template <class TVrtIter, class TAPosition>
UG_API 
typename TAPosition::ValueType
CalculateCenter(TVrtIter vrtsBegin, TVrtIter vrtsEnd,
				Grid::AttachmentAccessor<Vertex, TAPosition>& aaPos);

////////////////////////////////////////////////////////////////////////
//	CalculateBarycenter
/// calculates the barycenter of a set of vertices
template <class TVrtIter, class TAPosition>
UG_API 
typename TAPosition::ValueType
CalculateBarycenter(TVrtIter vrtsBegin, TVrtIter vrtsEnd,
					Grid::VertexAttachmentAccessor<TAPosition>& aaPos);

////////////////////////////////////////////////////////////////////////
//	MergeVertices
///	merges two vertices and restructures the adjacent elements.
/**
 * Since vertex v2 has to be removed in the process, the associated elements
 * of this vertex have to be replaced by new ones. Values attached to
 * old elements are passed on to the new ones using grid::pass_on_values.
 */
UG_API 
void MergeVertices(Grid& grid, Vertex* v1, Vertex* v2);

////////////////////////////////////////////////////////////////////////
///	Merges all vertices between the given iterators into a single vertex.
/**	Note that connected elements may be removed or replaced during this process.
 * The method returns the remaining vertex in the given list (*vrtsBegin).*/
template <class TVrtIterator>
UG_API 
Vertex* MergeMultipleVertices(Grid& grid, TVrtIterator vrtsBegin,
						  	  	  TVrtIterator vrtsEnd);

////////////////////////////////////////////////////////////////////////
//	RemoveDoubles
///	merges all vertices that are closer to each other than the specified threshold.
/**	The current implementation sadly enforces some restrictions to the
 *	container from which iterBegin and iterEnd stem. Only Grid, MultiGrid,
 *	Selector, MGSelector, SubsetHandler, MGSubsetHandler and similar containers
 *	are allowed. This is due to an implementation detail in the algorithm
 *	that should be removed in future revisins.
 *
 *	\todo	remove container restrictions as described above.
 * \{ */
template <int dim, class TVrtIterator>
UG_API 
void RemoveDoubles(Grid& grid, const TVrtIterator& iterBegin,
					const TVrtIterator& iterEnd, Attachment<MathVector<dim> >& aPos,
					number threshold);

template <int dim, class TVrtIterator, class TAAPos>
UG_API 
void RemoveDoubles(Grid& grid, const TVrtIterator& iterBegin,
					const TVrtIterator& iterEnd,
					TAAPos aaPos,
					number threshold);
/** \} */

////////////////////////////////////////////////////////////////////////
///	returns whether a vertex lies on the boundary of a polygonal chain.
/** The polygonal chain may be part of a bigger grid containing faces
 *	and volume elements. To distinguish which edges should be part of
 *	the polygonal chain, you may specify a callback to identify them.
 */
UG_API 
bool IsBoundaryVertex1D(Grid& grid, Vertex* v,
						Grid::edge_traits::callback cbConsiderEdge = ConsiderAll());

////////////////////////////////////////////////////////////////////////
///	returns whether a vertex lies on the boundary of a 2D grid.
/** A vertex is regarded as a 2d boundary vertex if it lies on a
 * 2d boundary edge.
 * if EDGEOPT_STORE_ASSOCIATED_FACES and VRTOPT_STORE_ASSOCIATED_EDGES
 * are enabled, the algorithm will be faster.
 */
UG_API bool IsBoundaryVertex2D(Grid& grid, Vertex* v);

////////////////////////////////////////////////////////////////////////
///	returns true if a vertex lies on the boundary of a 3D grid.
/** A vertex is regarded as a 3d boundary vertex if it lies on a
 * 3d boundary face.
 * if FACEOPT_STORE_ASSOCIATED_VOLUMES and VRTOPT_STORE_ASSOCIATED_FACES
 * are enabled, the algorithm will be faster.
*/
UG_API bool IsBoundaryVertex3D(Grid& grid, Vertex* v);

////////////////////////////////////////////////////////////////////////////////
///	returns true, if the vertex lies on a 1d, 2d or 3d boundary
UG_API bool LiesOnBoundary(Grid& grid, Vertex* v);

////////////////////////////////////////////////////////////////////////
//	IsRegularSurfaceVertex
///	returns true if the vertex lies inside a regular surface
/**
 * This algorithm indirectly uses Grid::mark.
 *	
 * The vertex is regarded as a regular surface vertex, if all associated
 * edges are connected to exactly 2 faces.
 */
UG_API bool IsRegularSurfaceVertex(Grid& grid, Vertex* v);

////////////////////////////////////////////////////////////////////////
/**
 * Uses Grid::mark()
 *
 * Vertices that are adjacent with more than two crease-edges are
 * regarded as a fixed vertex.
 */
UG_API 
void MarkFixedCreaseVertices(Grid& grid, SubsetHandler& sh,
							int creaseSI, int fixedSI);

////////////////////////////////////////////////////////////////////////
///	returns the position of the vertex.
/**	Main purpose is to allow the use of vertices in template-methods
 *	that call CalculateCenter
 *	\{ */
template<class TVertexPositionAttachmentAccessor>
UG_API 
inline
typename TVertexPositionAttachmentAccessor::ValueType
CalculateCenter(const Vertex* v, TVertexPositionAttachmentAccessor& aaPosVRT);

template<class TAAPosVRT, class TAAWeightVRT>
UG_API
typename TAAPosVRT::ValueType
CalculateCenter(const Vertex* v, TAAPosVRT& aaPos, TAAWeightVRT& aaWeight);
/** \} */


////////////////////////////////////////////////////////////////////////
///	transforms a vertex by a given matrix
template<class TAAPos> inline
UG_API 
void TransformVertex(Vertex* vrt, matrix33& m, TAAPos& aaPos);

////////////////////////////////////////////////////////////////////////
///	transforms all given vertices by a given matrix
template<class TIterator, class TAAPos> inline
UG_API 
void TransformVertices(TIterator vrtsBegin, TIterator vrtsEnd,
					   matrix33& m, TAAPos& aaPos);

////////////////////////////////////////////////////////////////////////
///	moves vertices by the specified offset
template<class TIterator, class TAAPos> inline
void MoveVertices(TIterator vrtsBegin, TIterator vrtsEnd, TAAPos aaPos,
				  const typename TAAPos::ValueType& offset);

////////////////////////////////////////////////////////////////////////
///	Returns true if the given point lies on the given vertex.
/**	The method simply checks, whether coordinates are identical. This method
 * is of no particular use, except for completeness regarding the use
 * of ConatinsPoint in template-programming.
 */
template <class vector_t, class TAAPos>
UG_API bool
ContainsPoint(const Vertex* v, const vector_t& p, TAAPos aaPos);

/// @} // end of doxygen defgroup command

}//	end of namespace

////////////////////////////////
//	include implementation
#include "vertex_util_impl.hpp"

#endif
