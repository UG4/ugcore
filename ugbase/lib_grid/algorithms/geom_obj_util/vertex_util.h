//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m02 d02

#ifndef __H__LIB_GRID__VERTEX_UTIL__
#define __H__LIB_GRID__VERTEX_UTIL__

#include <vector>
#include "lib_grid/lg_base.h"
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
int GetVertexIndex(EdgeBase* e, VertexBase* v);

////////////////////////////////////////////////////////////////////////
//	GetVertexIndex
///	returns the index at which vertex v is found in the given face
/**
 * returns -1 if the vertex was not found.
 */
int GetVertexIndex(Face* f, VertexBase* v);

////////////////////////////////////////////////////////////////////////
//	GetVertexIndex
///	returns the index at which vertex v is found in the given volume
/**
 * returns -1 if the vertex was not found.
 */
int GetVertexIndex(Volume* vol, VertexBase* v);

////////////////////////////////////////////////////////////////////////
//	GetConnectedVertex
///	returns the vertex that is connected to v via e.
/**
 * returns NULL if v is not contained in e.
 */
VertexBase* GetConnectedVertex(EdgeBase* e, VertexBase* v);

////////////////////////////////////////////////////////////////////////
//	GetConnectedVertex
///	returns the index of the first vertex that is contained in f and is not contained in e.
VertexBase* GetConnectedVertex(EdgeVertices* e, Face* f);

////////////////////////////////////////////////////////////////////////
//	GetConnectedVertexIndex
///	returns the index of the first vertex that is contained in the specified face and is not contained in the given edge.
int GetConnectedVertexIndex(Face* f, const EdgeDescriptor& ed);

////////////////////////////////////////////////////////////////////////
///	returns the edge in the triangle tri, which does not contain vrt.
/**	Make sure that tri is a triangle!*/
EdgeBase* GetConnectedEdge(Grid& g, VertexBase* vrt, Face* tri);

////////////////////////////////////////////////////////////////////////
//	CollectNeighbours
///	fills an array with all neighbour-vertices of v.
/**
 * v will not be contained in vNeighboursOut.
 * requires grid-option GRIDOPT_STANDARD_INTERCONNECTION.
 * This method is fast if grid-options FACEOPT_AUTOGENERATE_EDGES
 * and VOLOPT_AUTOGENERATE_EDGES are enabled - if there are any
 * faces and volumes.
 * It works without these options, too. The method however will
 * require more time in this case.
 */
void CollectNeighbours(std::vector<VertexBase*>& vNeighborsOut, Grid& grid, VertexBase* v);

////////////////////////////////////////////////////////////////////////
//	FindVertexByCoordiante
///	returns the vertex that is the closest to the given coordinate
/**
 * returns NULL if no vertex was found (if iterBegin == iterEnd).
 */
VertexBase* FindVertexByCoordiante(vector3& coord, VertexBaseIterator iterBegin, VertexIterator iterEnd,
									Grid::VertexAttachmentAccessor<APosition>& aaPos);

////////////////////////////////////////////////////////////////////////
///	calculates the normal of a vertex using associated faces
/**
 * TAAPosVRT has to be an attachment accessor for the vector3 type that
 * works on the triangles in grid.
 */
template <class TAAPosVRT>
void CalculateVertexNormal(vector3& nOut, Grid& grid, VertexBase* vrt,
						   TAAPosVRT& aaPos);

////////////////////////////////////////////////////////////////////////
//	CalculateVertexNormals
///	calculates the normals of all vertices in grid and stores them in aNorm.
/**
 * aPos has to be attached to grid.
 * If some attachments were not attached correctly, the method returns false.
 */
bool CalculateVertexNormals(Grid& grid, APosition& aPos, ANormal& aNorm);

////////////////////////////////////////////////////////////////////////
//	CalculateBoundingBox
/// calculates the BoundingBox
void CalculateBoundingBox(vector3& vMinOut, vector3& vMaxOut, VertexBaseIterator vrtsBegin,
						  VertexBaseIterator vrtsEnd, Grid::VertexAttachmentAccessor<AVector3>& aaPos);

////////////////////////////////////////////////////////////////////////
//	CalculateBarycenter - mstepnie
/// calculates the barycenter of a set of vertices
vector3 CalculateBarycenter(VertexBaseIterator vrtsBegin, VertexBaseIterator vrtsEnd,
							Grid::VertexAttachmentAccessor<AVector3>& aaPos);

////////////////////////////////////////////////////////////////////////
//	MergeVertices
///	merges two vertices and restructures the adjacent elements.
/**
 * Since vertex v2 has to be removed in the process, the associated elements
 * of this vertex have to be replaced by new ones. Values attached to
 * old elements are passed on to the new ones using grid::pass_on_values.
 */
void MergeVertices(Grid& grid, VertexBase* v1, VertexBase* v2);

////////////////////////////////////////////////////////////////////////
//	RemoveDoubles
///	merges all vertices that are closer to each other than the specified threshold.
template <int dim>
void RemoveDoubles(Grid& grid, const VertexBaseIterator& iterBegin,
					const VertexBaseIterator& iterEnd, Attachment<MathVector<dim> >& aPos,
					number threshold);

////////////////////////////////////////////////////////////////////////
///	returns whether a vertex lies on the boundary of a 2D grid.
/** A vertex is regarded as a 2d boundary vertex if it lies on a
 * 2d boundary edge.
 * if EDGEOPT_STORE_ASSOCIATED_FACES and VRTOPT_STORE_ASSOCIATED_EDGES
 * are enabled, the algorithm will be faster.
 */
bool IsBoundaryVertex2D(Grid& grid, VertexBase* v);

////////////////////////////////////////////////////////////////////////
///	returns true if a vertex lies on the boundary of a 3D grid.
/** A vertex is regarded as a 3d boundary vertex if it lies on a
 * 3d boundary face.
 * if FACEOPT_STORE_ASSOCIATED_VOLUMES and VRTOPT_STORE_ASSOCIATED_FACES
 * are enabled, the algorithm will be faster.
*/
bool IsBoundaryVertex3D(Grid& grid, VertexBase* v);

////////////////////////////////////////////////////////////////////////
/**
 * Uses Grid::mark()
 *
 * Vertices that are adjacent with more than two crease-edges are
 * regarded as a fixed vertex.
 */
void MarkFixedCreaseVertices(Grid& grid, SubsetHandler& sh,
							int creaseSI, int fixedSI);

////////////////////////////////////////////////////////////////////////
template <class TIterator, class AAPosVRT>
void LaplacianSmooth(Grid& grid, TIterator vrtsBegin,
					TIterator vrtsEnd, AAPosVRT& aaPos,
					number alpha, int numIterations);

////////////////////////////////////////////////////////////////////////
///	returns the position of the vertex.
/**	Main purpose is to allow the use of vertices in template-methods
 *	that call CalculateCenter*/
template<class TVertexPositionAttachmentAccessor>
inline
typename TVertexPositionAttachmentAccessor::ValueType
CalculateCenter(VertexBase* v, TVertexPositionAttachmentAccessor& aaPosVRT);

/// @} // end of doxygen defgroup command

}//	end of namespace

////////////////////////////////
//	include implementation
#include "vertex_util_impl.hpp"

#endif
