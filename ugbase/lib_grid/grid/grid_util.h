//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m10 d16

#ifndef __H__LIB_GRID__GRID_UTIL__
#define __H__LIB_GRID__GRID_UTIL__

#include <vector>
#include "grid.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
//	IsVolumeBoundaryFace
///	returns true if the given face is a boundary face.
bool IsVolumeBoundaryFace(Grid& grid, Face* f);

////////////////////////////////////////////////////////////////////////
///	returns whether an edge lies on the boundary of a 2D grid.
/**
 * if EDGEOPT_STORE_ASSOCIATED_FACES is enabled, the algorithm will be faster.
 */
bool IsBoundaryEdge2D(Grid& grid, EdgeBase* e);

////////////////////////////////////////////////////////////////////////
///	returns whether a vertex lies on the boundary of a 2D grid.
/**
 * if EDGEOPT_STORE_ASSOCIATED_FACES and VRTOPT_STORE_ASSOCIATED_EDGES
 * are enabled, the algorithm will be faster.
 */
bool IsBoundaryVertex2D(Grid& grid, VertexBase* v);

////////////////////////////////////////////////////////////////////////
//	GetEdge
///	returns the edge of the given grid, that connects the given points.
/**
 * This method requires the grid-option VRTOPT_STORE_ASSOCIATED_EDGES.
 * returns NULL if grid does not contain a matching edge.
 */
EdgeBase* FindEdge(Grid& grid, VertexBase* vrt1, VertexBase* vrt2);

inline EdgeBase* FindEdge(Grid& grid, EdgeDescriptor& ed)
	{
		return FindEdge(grid, ed.vertex(0), ed.vertex(1));
	}
///	returns the i-th edge of the given face
/**
 * This method requires the grid-option VRTOPT_STORE_ASSOCIATED_EDGES.
 * returns NULL if grid does not contain a matching edge.
 */
EdgeBase* FindEdge(Grid& grid, Face* f, uint ind);

///	returns the i-th edge of the given face
/**
 * This method requires the grid-option VRTOPT_STORE_ASSOCIATED_EDGES.
 * returns NULL if grid does not contain a matching edge.
 */
EdgeBase* FindEdge(Grid& grid, Volume* v, uint ind);

////////////////////////////////////////////////////////////////////////
//	CollectEdges
///	Collects all edges that exist in the given grid are part of the given face.
/**
 * if FACEOPT_STORE_ASSOCIATED_EDGES is enabled then the algorithm uses
 * associated_edges_begin(f), associated_edges_end(f)
 * in order to find the edges.
 * If not, GetEdge is used. This possibly involves auto-enabling of
 * VRTOPT_STORE_ASSOCIATED_EDGES.
 * The second option performs worse!
 */
void CollectEdges(std::vector<EdgeBase*>& vEdgesOut, Grid& grid, Face* f, bool clearContainer = true);

///	Collects all edges that exist in the given grid are part of the given volume.
/**
 * if VOLOPT_STORE_ASSOCIATED_EDGES is enabled then the algorithm uses
 * associated_edges_begin(v), associated_edges_end(v)
 * in order to find the edges.
 * If not, GetEdge is used. This possibly involves auto-enabling of
 * VRTOPT_STORE_ASSOCIATED_EDGES.
 * The second option performs worse!
 */
void CollectEdges(std::vector<EdgeBase*>& vEdgesOut, Grid& grid, Volume* v, bool clearContainer = true);

////////////////////////////////////////////////////////////////////////
//	EdgeContains
///	returns true if the given edge contains the given vertices
bool EdgeContains(EdgeBase* e, VertexBase* vrt1, VertexBase* vrt2);

////////////////////////////////////////////////////////////////////////
//	EdgeMatches
///	returns true if the given edge matches the given EdgeDescriptor
bool EdgeMatches(EdgeBase* e, EdgeDescriptor& ed);

////////////////////////////////////////////////////////////////////////
//	FindFace
///	returns the face that matches the face descriptor.
/**
 * This method requires the grid-option VRTOPT_STORE_ASSOCIATED_FACES
 * returns NULL if grid does not contain a matching face.
 * calculates hash values for faces and thus reduces calls to MatchFace.
 */
Face* FindFace(Grid& grid, FaceDescriptor& fd);

///	returns the i-th face from the given volume
/**
 * This method requires the grid-option VRTOPT_STORE_ASSOCIATED_FACES
 * returns NULL if grid does not contain a matching face.
 */
Face* FindFace(Grid& grid, Volume* v, uint ind, bool ignoreAssociatedFaces = false);

////////////////////////////////////////////////////////////////////////
//	CollectFaces
///	Collects all faces that exist in the given grid which contain the given edge.
/**
 * if EDGEOPT_STORE_ASSOCIATED_FACES is enabled then the algorithm uses
 * associated_faces_begin(e), associated_faces_end(e)
 * in order to find the faces.
 * if not the algorithm iterates over all faces associated with one of the edges
 * end-points and returns each face which contains the edge.
 */
void CollectFaces(std::vector<Face*>& vFacesOut, Grid& grid, EdgeBase* e, bool clearContainer = true);

///	Collects all faces that exist in the given grid are part of the given volume.
/**
 * if VOLOPT_STORE_ASSOCIATED_VOLUMES is enabled then the algorithm uses
 * associated_faces_begin(v), associated_faces_end(v)
 * in order to find the faces.
 * If not, GetFace is used. This possibly involves auto-enabling of
 * VRTOPT_STORE_ASSOCIATED_FACES.
 * The second option performs worse!
 */
void CollectFaces(std::vector<Face*>& vFacesOut, Grid& grid, Volume* v, bool clearContainer = true);

////////////////////////////////////////////////////////////////////////
//	FaceMatches
///	returns true if the given face contains exactly the same points as the given descriptor.
bool FaceMatches(Face* f, FaceDescriptor& fd);

////////////////////////////////////////////////////////////////////////
//	FaceContains
///	returns true if the given face contains the given vertex
bool FaceContains(Face* f, VertexBase* v);

////////////////////////////////////////////////////////////////////////
//	FaceContains
///	returns true if the given face contains the given edge
bool FaceContains(Face* f, EdgeBase* e);

////////////////////////////////////////////////////////////////////////
//	CollectVolumes
///	Collects all volumes that exist in the given grid which contain the given edge.
/**
 * if EDGEOPT_STORE_ASSOCIATED_VOLUMES is enabled then the algorithm uses
 * associated_volumes_begin(e), associated_volumes_end(e)
 * in order to find the volumes.
 * if not the algorithm iterates over all volumes associated with one of the edges
 * end-points and returns each volume that contains the edge.
 */
void CollectVolumes(std::vector<Volume*>& vVolumesOut, Grid& grid, EdgeBase* e, bool clearContainer = true);

////////////////////////////////////////////////////////////////////////
//	FindVolume
///	returns the volume that matches the volume descriptor.
/**
 * This method requires the grid-option VRTOPT_STORE_ASSOCIATED_VOLUMES
 * returns NULL if grid does not contain a matching volume.
 * calculates hash values for volumes and thus reduces calls to MatchVolume.
 */
Volume* FindVolume(Grid& grid, VolumeDescriptor& vd);

///	Collects all volumes that exist in the given grid which contain the given face.
/**
 * if FACEOPT_STORE_ASSOCIATED_VOLUMES is enabled then the algorithm uses
 * associated_volumes_begin(f), associated_volumes_end(f)
 * in order to find the volumes.
 * if not the algorithm iterates over all volumes associated with one of the faces
 * end-points and returns each volume that contains the face.
 */
void CollectVolumes(std::vector<Volume*>& vVolumesOut, Grid& grid, Face* f, bool clearContainer = true, bool ignoreAssociatedVolumes = false);

void CollectVolumes(std::vector<Volume*>& vVolumesOut, Grid& grid, FaceDescriptor& fd, bool clearContainer = true);

////////////////////////////////////////////////////////////////////////
//	VolumeContains
///	returns true if the given volume contains the given vertex
bool VolumeContains(Volume* v, VertexBase* vrt);

///	returns true if the given volume contains the given edge
bool VolumeContains(Volume* v, EdgeBase* e);

///	returns true if the given volume contains the given face
bool VolumeContains(Volume* v, Face* f);
bool VolumeContains(Volume* v, FaceDescriptor& fd);

}//	end of namespace libGrid

#endif
