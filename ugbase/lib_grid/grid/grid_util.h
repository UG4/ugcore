//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m10 d16

#ifndef __H__LIB_GRID__GRID_UTIL__
#define __H__LIB_GRID__GRID_UTIL__

#include <vector>
#include "geometric_base_objects.h"

namespace ug
{

class Grid;


////////////////////////////////////////////////////////////////////////
//	CompareVertices
///	Checks whether ev1 and ev2 contain the same vertices.
/**
 * \ingroup lib_grid_algorithms_edge_util
 *
 * Can be used to compare EdgeBase with EdgeBase,
 * EdgeDescriptor with EdgeDescriptor or
 * EdgeBase with EdgeDescriptor.
 */
bool CompareVertices(const EdgeVertices* ev1,
					const EdgeVertices* ev2);

///	Checks whether fv1 and fv2 contain the same vertices.
/**
 * \ingroup lib_grid_algorithms_face_util
 *
 * Can be used to compare Face with Face,
 * FaceDescriptor with FaceDescriptor or
 * Face with FaceDescriptor.
 *
 * Before calling this method one should consider to compare the
 * hashes of fv1 and fv2 (if(hash_key(fv1) == hash_key(fv2))...)
 */
bool CompareVertices(const FaceVertices* fv1,
					const FaceVertices* fv2);

///	Checks whether vv1 and vv2 contain the same vertices.
/**
 * \ingroup lib_grid_algorithms_volume_util
 *
 * Can be used to compare Volume with Volume,
 * VolumeDescriptor with VolumeDescriptor or
 * Volume with VolumeDescriptor.
 *
 * Before calling this method one should consider to compare the
 * hashes of vv1 and vv2 (if(hash_key(vv1) == hash_key(vv2))...)
 */
bool CompareVertices(const VolumeVertices* vv1,
					const VolumeVertices* vv2);


////////////////////////////////////////////////////////////////////////
//	CompareVertexContainer
///	compares vertices in a container
/**
 * \ingroup lib_grid_algorithms_vertex_util
 *
 *	If you want to compare EdgeBase, Face, Volume, EdgeDescriptor,
 *	FaceDescriptor or VolumeDescriptor with other EdgeBase, ... then please
 *	use CompareVertices instead.
 *
 *	TVrtContainer has to feature the following methods (or compatible ones):
 *	int size(): has to return the numbe of vertices in the container.
 *	VertexBase* operator[](int index): has to return the i-th vertex.
 *
 *	Good types would be:
 *	EdgeVertices, FaceVertices, VolumeVertices,
 *	EdgeBase, Face, Volume, ...,
 *	std::vector<VertexBase*>, ...
 *
 *	returns true if con1 and con2 contain the same vertices.
 */
template <class TVrtContainer1, class TVrtContainer2>
bool CompareVertexContainer(const TVrtContainer1& con1,
					const TVrtContainer2& con2);

////////////////////////////////////////////////////////////////////////
//	CollectVertices
///	Dummy-method. Puts the vertex itself into the given vector.
/**
 * \ingroup lib_grid_algorithms_vertex_util
 *
 * This function simply returns the vertex itself. It is added for completeness,
 * such that the function can be used in template code.
 * \{
 */
void CollectVertices(std::vector<VertexBase*>& vVertexOut, Grid& grid, VertexBase* v, bool clearContainer = true);

inline void CollectAssociated(std::vector<VertexBase*>& vVertexOut,
					  Grid& grid, VertexBase* v, bool clearContainer = true);
/** \} */

////////////////////////////////////////////////////////////////////////
//	CollectVertices
///	Collects all vertices that are part of the given edge
/**
 * \ingroup lib_grid_algorithms_edge_util
 *
 * This function returns a std::vector of pointers to all vertices,
 * that are part of the given edge.
 * \{
 */
void CollectVertices(std::vector<VertexBase*>& vVertexOut, Grid& grid, EdgeBase* e, bool clearContainer = true);

inline void CollectAssociated(std::vector<VertexBase*>& vVertexOut,
					   Grid& grid, EdgeBase* e, bool clearContainer = true);
/** \} */

////////////////////////////////////////////////////////////////////////
//	CollectVertices
///	Collects all vertices that are part of the given face
/**
 * \ingroup lib_grid_algorithms_face_util
 *
 * This function returns a std::vector of pointers to all vertices,
 * that are part of the given face.
 * \{
 */
void CollectVertices(std::vector<VertexBase*>& vVertexOut, Grid& grid, Face* f, bool clearContainer = true);

inline void CollectAssociated(std::vector<VertexBase*>& vVertexOut,
					   Grid& grid, Face* f, bool clearContainer = true);
/** \} */

////////////////////////////////////////////////////////////////////////
//	CollectVertices
///	Collects all vertices that are part of the given volume
/**
 * \ingroup lib_grid_algorithms_volume_util
 *
 * This function returns a std::vector of pointers to all vertices,
 * that are part of the given volume.
 * \{
 */
void CollectVertices(std::vector<VertexBase*>& vVertexOut, Grid& grid, Volume* v, bool clearContainer = true);

inline void CollectAssociated(std::vector<VertexBase*>& vVertexOut,
					   Grid& grid, Volume* v, bool clearContainer = true);
/** \} */

////////////////////////////////////////////////////////////////////////
//	GetVertex
///	returns the i'th vertex of a vertex
/**
 * \ingroup lib_grid_algorithms_vertex_util
 *
 * This function simply returns the vertex itself. It is added for completeness,
 * such that the function can be used in template code.
 */
inline VertexBase* GetVertex(VertexBase* v, size_t i);

///	returns the i'th vertex of an edge
/**
 * \ingroup lib_grid_algorithms_edge_util
 *
 * This function simply returns the i'th vertex of an edge
 */
inline VertexBase* GetVertex(EdgeBase* e, size_t i);

///	returns the i'th vertex of a face
/**
 * \ingroup lib_grid_algorithms_face_util
 *
 * This function simply returns the i'th vertex of a face
 */
inline VertexBase* GetVertex(Face* f, size_t i);

///	returns the i'th vertex of a volume
/**
 * \ingroup lib_grid_algorithms_volume_util
 *
 * This function simply returns the i'th vertex of a volume
 */
inline VertexBase* GetVertex(Volume* v, size_t i);

////////////////////////////////////////////////////////////////////////
//	CollectEdgesSorted
///	Collects all edges that exist in the given grid are part of the given edge in the order defined by the reference elements.
/**
 * \ingroup lib_grid_algorithms_edge_util
 *
 * This function simply returns the edge itself. It is added for completeness,
 * such that the function can be used in template code.
 */
void CollectEdgesSorted(std::vector<EdgeBase*>& vEdgesOut, Grid& grid, EdgeBase* e, bool clearContainer = true);

///	Collects all edges that exist in the given grid are part of the given face in the order defined by the reference elements.
/**
 * \ingroup lib_grid_algorithms_face_util
 *
 * This function returns the associated edges of a face in an std::vector. The order of the edges in the vector
 * is equal to the numbering of edges in the reference element for the face.
 */
void CollectEdgesSorted(std::vector<EdgeBase*>& vEdgesOut, Grid& grid, Face* f, bool clearContainer = true);

///	Collects all edges that exist in the given grid are part of the given volume in the order defined by the reference elements
/**
 * \ingroup lib_grid_algorithms_volume_util
 *
 * This function returns the associated edges of a volume in an std::vector. The order of the edges in the vector
 * is equal to the numbering of edges in the reference element for the volume.
 */
void CollectEdgesSorted(std::vector<EdgeBase*>& vEdgesOut, Grid& grid, Volume* v, bool clearContainer = true);

////////////////////////////////////////////////////////////////////////
//	CollectEdges
///	Collects all edges that exist in the given grid are part of the given edge.
/**
 * \ingroup lib_grid_algorithms_edge_util
 * \{
 */
void CollectEdges(std::vector<EdgeBase*>& vEdgesOut, Grid& grid, VertexBase* vrt, bool clearContainer = true);

inline void CollectAssociated(std::vector<EdgeBase*>& vEdgesOut,
					Grid& grid, VertexBase* vrt, bool clearContainer = true);
/** \} */

///	Collects all edges that exist in the given grid are part of the given edge.
/**
 * \ingroup lib_grid_algorithms_edge_util
 *
 * This function simply returns the edge itself. It is added for completeness,
 * such that the function can be used in template code.
 * \{
 */
void CollectEdges(std::vector<EdgeBase*>& vEdgesOut, Grid& grid, EdgeBase* e, bool clearContainer = true);

inline void CollectAssociated(std::vector<EdgeBase*>& vEdgesOut,
					Grid& grid, EdgeBase* e, bool clearContainer = true);
/** \} */

////////////////////////////////////////////////////////////////////////
//	CollectEdges
///	Collects all edges that exist in the given grid are part of the given face.
/**
 * \ingroup lib_grid_algorithms_face_util
 *
 * if FACEOPT_STORE_ASSOCIATED_EDGES is enabled then the algorithm uses
 * associated_edges_begin(f), associated_edges_end(f)
 * in order to find the edges.
 * If not, GetEdge is used. This possibly involves auto-enabling of
 * VRTOPT_STORE_ASSOCIATED_EDGES.
 * The second option performs worse!
 * \{
 */
void CollectEdges(std::vector<EdgeBase*>& vEdgesOut, Grid& grid, Face* f, bool clearContainer = true);

inline void CollectAssociated(std::vector<EdgeBase*>& vEdgesOut,
						Grid& grid, Face* f, bool clearContainer = true);
/** \} */

///	Collects all edges that exist in the given grid are part of the given volume.
/**
 * \ingroup lib_grid_algorithms_volume_util
 *
 * if VOLOPT_STORE_ASSOCIATED_EDGES is enabled then the algorithm uses
 * associated_edges_begin(v), associated_edges_end(v)
 * in order to find the edges.
 * If not, GetEdge is used. This possibly involves auto-enabling of
 * VRTOPT_STORE_ASSOCIATED_EDGES.
 * The second option performs worse!
 * \{
 */
void CollectEdges(std::vector<EdgeBase*>& vEdgesOut, Grid& grid, Volume* v, bool clearContainer = true);

inline void CollectAssociated(std::vector<EdgeBase*>& vEdgesOut,
						Grid& grid, Volume* v, bool clearContainer = true);
/** \} */

////////////////////////////////////////////////////////////////////////
//	EdgeContains
///	returns true if the given edge contains the given vertex
/// \ingroup lib_grid_algorithms_edge_util
inline bool EdgeContains(EdgeVertices* e, VertexBase* vrt);

///	returns true if the given edge contains the given vertices
/// \ingroup lib_grid_algorithms_edge_util
inline bool EdgeContains(EdgeVertices* e, VertexBase* vrt1, VertexBase* vrt2);

////////////////////////////////////////////////////////////////////////
//	CollectFacesSorted
///	Collects all face that exist in the given grid are part of the given edge in the order defined by the reference elements.
/**
 * \ingroup lib_grid_algorithms_edge_util
 */
void CollectFacesSorted(std::vector<Face*>& vFacesOut, Grid& grid, EdgeBase* e, bool clearContainer = true);

///	Collects all face that exist in the given grid are part of the given face in the order defined by the reference elements.
/**
 * \ingroup lib_grid_algorithms_face_util
 */
void CollectFacesSorted(std::vector<Face*>& vFacesOut, Grid& grid, Face* f, bool clearContainer = true);

///	Collects all face that exist in the given grid are part of the given volume in the order defined by the reference elements.
/**
 * \ingroup lib_grid_algorithms_volume_util
 *
 * This function returns the associated faces of a volume in an std::vector. The order of the faces in the vector
 * is equal to the numbering of faces in the reference element for the volume.
 */
void CollectFacesSorted(std::vector<Face*>& vFacesOut, Grid& grid, Volume* v, bool clearContainer = true);

////////////////////////////////////////////////////////////////////////
//	CollectFaces
/**
 * \ingroup lib_grid_algorithms_face_util
 *
 * \{
 */
void CollectFaces(std::vector<Face*>& vFacesOut, Grid& grid, VertexBase* vrt, bool clearContainer = true);

inline void CollectAssociated(std::vector<Face*>& vFacesOut,
					Grid& grid, VertexBase* vrt, bool clearContainer = true);
/** \} */

///	Collects all faces that exist in the given grid which contain the given edge.
/**
 * \ingroup lib_grid_algorithms_edge_util
 *
 * if EDGEOPT_STORE_ASSOCIATED_FACES is enabled then the algorithm uses
 * associated_faces_begin(e), associated_faces_end(e)
 * in order to find the faces.
 * if not the algorithm iterates over all faces associated with one of the edges
 * end-points and returns each face which contains the edge.
 * \{
 */
void CollectFaces(std::vector<Face*>& vFacesOut, Grid& grid, EdgeBase* e, bool clearContainer = true);

inline void CollectAssociated(std::vector<Face*>& vFacesOut,
					Grid& grid, EdgeBase* e, bool clearContainer = true);
/** \} */

///	Collects all faces that exist in the given grid are part of the given volume.
/**
 * \ingroup lib_grid_algorithms_volume_util
 *
 * if VOLOPT_STORE_ASSOCIATED_VOLUMES is enabled then the algorithm uses
 * associated_faces_begin(v), associated_faces_end(v)
 * in order to find the faces.
 * If not, GetFace is used. This possibly involves auto-enabling of
 * VRTOPT_STORE_ASSOCIATED_FACES.
 * The second option performs worse!
 * \{
 */
void CollectFaces(std::vector<Face*>& vFacesOut, Grid& grid, Volume* v, bool clearContainer = true);

inline void CollectAssociated(std::vector<Face*>& vFacesOut,
					Grid& grid, Volume* v, bool clearContainer = true);

/** \} */

////////////////////////////////////////////////////////////////////////
//	FaceMatches
///	returns true if the given face contains exactly the same points as the given descriptor.
//bool FaceMatches(Face* f, FaceDescriptor& fd);

////////////////////////////////////////////////////////////////////////
//	FaceContains
///	returns true if the given face contains the given vertex
/// \ingroup lib_grid_algorithms_face_util
bool FaceContains(Face* f, VertexBase* v);

////////////////////////////////////////////////////////////////////////
//	FaceContains
///	returns true if the given face contains the given edge
/// \ingroup lib_grid_algorithms_face_util
bool FaceContains(Face* f, EdgeVertices* ev);

////////////////////////////////////////////////////////////////////////
//	CollectVolumes
///	Collects all volumes that exist in the given grid which contain the given vertex.
/**
 * \ingroup lib_grid_algorithms_vertex_util
 */
void CollectVolumes(std::vector<Volume*>& vVolumesOut, Grid& grid, VertexBase* vrt, bool clearContainer = true);

inline void CollectAssociated(std::vector<Volume*>& vVolumesOut,
					Grid& grid, VertexBase* vrt, bool clearContainer = true);
/** \} */

////////////////////////////////////////////////////////////////////////
//	CollectVolumes
///	Collects all volumes that exist in the given grid which contain the given edge.
/**
 * \ingroup lib_grid_algorithms_edge_util
 *
 * if EDGEOPT_STORE_ASSOCIATED_VOLUMES is enabled then the algorithm uses
 * associated_volumes_begin(e), associated_volumes_end(e)
 * in order to find the volumes.
 * if not the algorithm iterates over all volumes associated with one of the edges
 * end-points and returns each volume that contains the edge.
 */
void CollectVolumes(std::vector<Volume*>& vVolumesOut, Grid& grid, EdgeBase* e, bool clearContainer = true);

inline void CollectAssociated(std::vector<Volume*>& vVolumesOut,
					Grid& grid, EdgeBase* e, bool clearContainer = true);
/** \} */

/// @{
///	Collects all volumes that exist in the given grid which contain the given face.
/**
 * \ingroup lib_grid_algorithms_face_util
 *
 * if FACEOPT_STORE_ASSOCIATED_VOLUMES is enabled then the algorithm uses
 * associated_volumes_begin(f), associated_volumes_end(f)
 * in order to find the volumes.
 * if not the algorithm iterates over all volumes associated with one of the faces
 * end-points and returns each volume that contains the face.
 */
void CollectVolumes(std::vector<Volume*>& vVolumesOut, Grid& grid, Face* f, bool clearContainer = true, bool ignoreAssociatedVolumes = false);
inline void CollectAssociated(std::vector<Volume*>& vVolumesOut,
					Grid& grid, Face* f, bool clearContainer = true,
					bool ignoreAssociatedVolumes = false);

void CollectVolumes(std::vector<Volume*>& vVolumesOut, Grid& grid, FaceDescriptor& fd, bool clearContainer = true);
inline void CollectAssociated(std::vector<Volume*>& vVolumesOut,
					Grid& grid, FaceDescriptor& fd, bool clearContainer = true);
/// @}

////////////////////////////////////////////////////////////////////////
//	VolumeContains
/// \ingroup lib_grid_algorithms_volume_util
///	returns true if the given volume contains the given vertex
bool VolumeContains(Volume* v, VertexBase* vrt);

/// \ingroup lib_grid_algorithms_volume_util
///	returns true if the given volume contains the given edge
bool VolumeContains(Volume* v, EdgeVertices* ev);

/// \ingroup lib_grid_algorithms_volume_util
///	returns true if the given volume contains the given face
bool VolumeContains(Volume* v, FaceVertices* f);

}//	end of namespace libGrid

//	include template implementation
#include "grid_util_impl.hpp"

#endif
