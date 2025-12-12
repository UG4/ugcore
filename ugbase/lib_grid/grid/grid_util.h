/*
 * Copyright (c) 2009-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__LIB_GRID__GRID_UTIL__
#define __H__LIB_GRID__GRID_UTIL__

#include <vector>

#include "grid_base_objects.h"

namespace ug {

class Grid;


////////////////////////////////////////////////////////////////////////
//	CompareVertices
///	Checks whether ev1 and ev2 contain the same vertices.
/**
 * \ingroup lib_grid_algorithms_edge_util
 *
 * Can be used to compare Edge with Edge,
 * EdgeDescriptor with EdgeDescriptor or
 * Edge with EdgeDescriptor.
 */
UG_API 
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
UG_API 
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
UG_API 
bool CompareVertices(const VolumeVertices* vv1,
					const VolumeVertices* vv2);


/**
 * @brief Checks whether two vertex groups contain the same vertices.
 *
 * Typical use case: The first vertex group is an element, the other
 * is a descriptor. In that capacity, this function works like the
 * CompareVertex functions; but it is a tad more general and works
 * also with, e.g., CustomVertexGroups.
 *
 * @note One might think about merging VertexGroupsMatch and CompareVertices.
 */
UG_API
bool VertexGroupsMatch(const IVertexGroup* elem, const IVertexGroup& desc);

/// specialization of VertexGroupsMatch for single Vertex
UG_API
bool VertexGroupsMatch(const Vertex* elem, const IVertexGroup& desc);


////////////////////////////////////////////////////////////////////////
//	CompareVertexContainer
///	compares vertices in a container
/**
 * \ingroup lib_grid_algorithms_vertex_util
 *
 *	If you want to compare Edge, Face, Volume, EdgeDescriptor,
 *	FaceDescriptor or VolumeDescriptor with other Edge, ... then please
 *	use CompareVertices instead.
 *
 *	TVrtContainer has to feature the following methods (or compatible ones):
 *	int size(): has to return the numbe of vertices in the container.
 *	Vertex* operator [] (int index): has to return the i-th vertex.
 *
 *	Good types would be:
 *	EdgeVertices, FaceVertices, VolumeVertices,
 *	Edge, Face, Volume, ...,
 *	std::vector<Vertex*>, ...
 *
 *	returns true if con1 and con2 contain the same vertices.
 */
template <typename TVrtContainer1, typename TVrtContainer2>
UG_API 
bool CompareVertexContainer(const TVrtContainer1& con1,
					const TVrtContainer2& con2);


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///	\addtogroup lib_grid_algorithms_neighborhood_util
///	@{

////////////////////////////////////////////////////////////////////////////////
//	CollectVertices
////////////////////////////////////////////////////////////////////////////////
UG_API 
inline void CollectAssociated(std::vector<Vertex*>& vVertexOut, Grid& grid,
                              GridObject* obj, bool clearContainer = true);
UG_API 
inline void CollectVertices(std::vector<Vertex*>& vVertexOut, Grid& grid,
                            GridObject* obj, bool clearContainer = true);

///	Dummy-method. Puts the vertex itself into the given vector.
/** This function simply returns the vertex itself. It is added for completeness,
 * such that the function can be used in template code.
 * \{
 */
UG_API 
void CollectVertices(std::vector<Vertex*>& vVertexOut, Grid& grid,
                     	 	 	 	 Vertex* v, bool clearContainer = true);

UG_API 
inline void CollectAssociated(std::vector<Vertex*>& vVertexOut,
					  Grid& grid, Vertex* v, bool clearContainer = true);
/** \} */

///	Collects all vertices that are part of the given edge
/** This function returns a std::vector of pointers to all vertices,
 * that are part of the given edge.
 * \{
 */
UG_API 
void CollectVertices(std::vector<Vertex*>& vVertexOut, Grid& grid,
                     	 	 	 	 Edge* e, bool clearContainer = true);

UG_API 
inline void CollectAssociated(std::vector<Vertex*>& vVertexOut,
                          Grid& grid, Edge* e, bool clearContainer = true);
/** \} */

///	Collects all vertices that are part of the given face
/** This function returns a std::vector of pointers to all vertices,
 * that are part of the given face.
 * \{
 */
UG_API 
void CollectVertices(std::vector<Vertex*>& vVertexOut, Grid& grid,
                     	 	 	 	 	 Face* f, bool clearContainer = true);

UG_API 
inline void CollectAssociated(std::vector<Vertex*>& vVertexOut,
                              Grid& grid, Face* f, bool clearContainer = true);
/** \} */

///	Collects all vertices that are part of the given volume
/** This function returns a std::vector of pointers to all vertices,
 * that are part of the given volume.
 * \{
 */
UG_API 
void CollectVertices(std::vector<Vertex*>& vVertexOut, Grid& grid,
                     	 	 	 	 	 Volume* v, bool clearContainer = true);

inline void CollectAssociated(std::vector<Vertex*>& vVertexOut,
					        Grid& grid, Volume* v, bool clearContainer = true);
/** \} */

//////////////////////////////////////////////////////////////////////////////
//	GetVertex
//////////////////////////////////////////////////////////////////////////////

///	returns the i'th vertex of a vertex
/**
 * \ingroup lib_grid_algorithms_vertex_util
 *
 * This function simply returns the vertex itself. It is added for completeness,
 * such that the function can be used in template code.
 */
inline Vertex* GetVertex(Vertex* v, size_t i);

///	returns the i'th vertex of an edge
/**
 * \ingroup lib_grid_algorithms_edge_util
 *
 * This function simply returns the i'th vertex of an edge
 */
inline Vertex* GetVertex(Edge* e, size_t i);

///	returns the i'th vertex of a face
/**
 * \ingroup lib_grid_algorithms_face_util
 *
 * This function simply returns the i'th vertex of a face
 */
inline Vertex* GetVertex(Face* f, size_t i);

///	returns the i'th vertex of a volume
/**
 * \ingroup lib_grid_algorithms_volume_util
 *
 * This function simply returns the i'th vertex of a volume
 */
inline Vertex* GetVertex(Volume* v, size_t i);

//////////////////////////////////////////////////////////////////////////////
// NumVertices
//////////////////////////////////////////////////////////////////////////////

///	Returns the number of vertices of the given geometric object.
/** \ingroup lib_grid_algorithms_volume_util
 * Valid template arguments are Vertex, Edge, Face, Volume and
 * derived types.
 *
 * unifies access to geometric objects.
 * \{*/
inline size_t NumVertices(Vertex* elem);
inline size_t NumVertices(Edge* elem);
inline size_t NumVertices(Face* elem);
inline size_t NumVertices(Volume* elem);
/**	\} */


//////////////////////////////////////////////////////////////////////////////
//	CollectEdgesSorted
//////////////////////////////////////////////////////////////////////////////

inline void CollectEdgesSorted(std::vector<Edge*>& vEdgesOut, Grid& grid,
                               GridObject* obj, bool clearContainer = true);

UG_API 
void CollectEdgesSorted(std::vector<Edge*>& vEdgesOut, Grid& grid,
                        			Vertex* v, bool clearContainer = true);

///	Collects all edges that exist in the given grid are part of the given edge in the order defined by the reference elements.
/** This function simply returns the edge itself. It is added for completeness,
 * such that the function can be used in template code.
 */
UG_API 
void CollectEdgesSorted(std::vector<Edge*>& vEdgesOut, Grid& grid,
                        				Edge* e, bool clearContainer = true);

///	Collects all edges that exist in the given grid are part of the given face in the order defined by the reference elements.
/** This function returns the associated edges of a face in an std::vector. The order of the edges in the vector
 * is equal to the numbering of edges in the reference element for the face.
 */
UG_API 
void CollectEdgesSorted(std::vector<Edge*>& vEdgesOut, Grid& grid,
                        					Face* f, bool clearContainer = true);

///	Collects all edges that exist in the given grid are part of the given volume in the order defined by the reference elements
/** This function returns the associated edges of a volume in an std::vector. The order of the edges in the vector
 * is equal to the numbering of edges in the reference element for the volume.
 */
UG_API 
void CollectEdgesSorted(std::vector<Edge*>& vEdgesOut, Grid& grid,
                        				Volume* v, bool clearContainer = true);

///////////////////////////////////////////////////////////////////////////////
//	CollectEdges
///////////////////////////////////////////////////////////////////////////////

inline void CollectAssociated(std::vector<Edge*>& vEdgesOut, Grid& grid,
                              GridObject* obj, bool clearContainer = true);

///	Collects all edges that exist in the given grid are part of the given edge.
/** \{
 */
UG_API 
void CollectEdges(std::vector<Edge*>& vEdgesOut, Grid& grid,
                  	  	  	  	  Vertex* vrt, bool clearContainer = true);

inline void CollectAssociated(std::vector<Edge*>& vEdgesOut,
					Grid& grid, Vertex* vrt, bool clearContainer = true);
/** \} */

///	Collects all edges that exist in the given grid are part of the given edge.
/** This function simply returns the edge itself. It is added for completeness,
 * such that the function can be used in template code.
 * \{
 */
UG_API 
void CollectEdges(std::vector<Edge*>& vEdgesOut, Grid& grid,
                  	  	  	  	  	  Edge* e, bool clearContainer = true);

inline void CollectAssociated(std::vector<Edge*>& vEdgesOut,
					Grid& grid, Edge* e, bool clearContainer = true);
/** \} */

///	Collects all edges that exist in the given grid are part of the given face.
/** if FACEOPT_STORE_ASSOCIATED_EDGES is enabled then the algorithm uses
 * associated_edges_begin(f), associated_edges_end(f)
 * in order to find the edges.
 * If not, GetEdge is used. This possibly involves auto-enabling of
 * VRTOPT_STORE_ASSOCIATED_EDGES.
 * The second option performs worse!
 * \{
 */
UG_API 
void CollectEdges(std::vector<Edge*>& vEdgesOut, Grid& grid,
                  	  	  	  	  	  	  Face* f, bool clearContainer = true);

inline void CollectAssociated(std::vector<Edge*>& vEdgesOut,
						Grid& grid, Face* f, bool clearContainer = true);
/** \} */

///	Collects all edges that exist in the given grid are part of the given volume.
/** if VOLOPT_STORE_ASSOCIATED_EDGES is enabled then the algorithm uses
 * associated_edges_begin(v), associated_edges_end(v)
 * in order to find the edges.
 * If not, GetEdge is used. This possibly involves auto-enabling of
 * VRTOPT_STORE_ASSOCIATED_EDGES.
 * The second option performs worse!
 * \{
 */
UG_API 
void CollectEdges(std::vector<Edge*>& vEdgesOut, Grid& grid,
                  	  	  	  	  	  	  Volume* v, bool clearContainer = true);

inline void CollectAssociated(std::vector<Edge*>& vEdgesOut,
						Grid& grid, Volume* v, bool clearContainer = true);
/** \} */

///////////////////////////////////////////////////////////////////////////////
//	EdgeContains
///////////////////////////////////////////////////////////////////////////////

///	returns true if the given edge contains the given vertex
/// \ingroup lib_grid_algorithms_edge_util
inline bool EdgeContains(EdgeVertices* e, Vertex* vrt);

///	returns true if the given edge contains the given vertices
/// \ingroup lib_grid_algorithms_edge_util
inline bool EdgeContains(EdgeVertices* e, Vertex* vrt1, Vertex* vrt2);

///////////////////////////////////////////////////////////////////////////////
//	CollectFacesSorted
///////////////////////////////////////////////////////////////////////////////

inline void CollectFacesSorted(std::vector<Face*>& vFacesOut, Grid& grid,
                            				GridObject* obj, bool clearContainer = true);


UG_API 
void CollectFacesSorted(std::vector<Face*>& vFacesOut, Grid& grid,
                        			Vertex* v, bool clearContainer = true);

///	Collects all face that exist in the given grid are part of the given edge in the order defined by the reference elements.
UG_API 
void CollectFacesSorted(std::vector<Face*>& vFacesOut, Grid& grid,
                        				Edge* e, bool clearContainer = true);

///	Collects all face that exist in the given grid are part of the given face in the order defined by the reference elements.
UG_API 
void CollectFacesSorted(std::vector<Face*>& vFacesOut, Grid& grid,
                        					Face* f, bool clearContainer = true);

///	Collects all face that exist in the given grid are part of the given volume in the order defined by the reference elements.
/** This function returns the associated faces of a volume in an std::vector. The order of the faces in the vector
 * is equal to the numbering of faces in the reference element for the volume.
 */
UG_API 
void CollectFacesSorted(std::vector<Face*>& vFacesOut, Grid& grid,
                        				Volume* v, bool clearContainer = true);

///////////////////////////////////////////////////////////////////////////////
//	CollectFaces
///////////////////////////////////////////////////////////////////////////////

inline void CollectAssociated(std::vector<Face*>& vFacesOut, Grid& grid,
                              GridObject* obj, bool clearContainer = true);

UG_API 
void CollectFaces(std::vector<Face*>& vFacesOut, Grid& grid,
                  	  	  	  	 Vertex* vrt, bool clearContainer = true);

inline void CollectAssociated(std::vector<Face*>& vFacesOut,
					Grid& grid, Vertex* vrt, bool clearContainer = true);

///	Collects all faces that exist in the given grid which contain the given edge.
/** if EDGEOPT_STORE_ASSOCIATED_FACES is enabled then the algorithm uses
 * associated_faces_begin(e), associated_faces_end(e)
 * in order to find the faces.
 * if not the algorithm iterates over all faces associated with one of the edges
 * end-points and returns each face which contains the edge.
 * \{
 */
UG_API 
void CollectFaces(std::vector<Face*>& vFacesOut, Grid& grid,
                  	  	  	  	  	  Edge* e, bool clearContainer = true);

inline void CollectAssociated(std::vector<Face*>& vFacesOut,
					Grid& grid, Edge* e, bool clearContainer = true);
/** \} */

///	Collects all faces. (Returns the face itself)
UG_API 
void CollectFaces(std::vector<Face*>& vFacesOut, Grid& grid,
                  	  	  	  	  	  	  Face* v, bool clearContainer = true);

///	Collects all faces. (Returns the face itself)
inline void CollectAssociated(std::vector<Face*>& vFacesOut,
					Grid& grid, Face* vol, bool clearContainer = true);

///	Collects all faces that exist in the given grid are part of the given volume.
/** if VOLOPT_STORE_ASSOCIATED_VOLUMES is enabled then the algorithm uses
 * associated_faces_begin(v), associated_faces_end(v)
 * in order to find the faces.
 * If not, GetFace is used. This possibly involves auto-enabling of
 * VRTOPT_STORE_ASSOCIATED_FACES.
 * The second option performs worse!
 * \{
 */
UG_API 
void CollectFaces(std::vector<Face*>& vFacesOut, Grid& grid,
                  	  	  	  	  	  	  Volume* v, bool clearContainer = true);

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
bool FaceContains(FaceVertices* f, Vertex* v);

////////////////////////////////////////////////////////////////////////
//	FaceContains
///	returns true if the given face contains the given edge
/// \ingroup lib_grid_algorithms_face_util
UG_API 
bool FaceContains(Face* f, EdgeVertices* ev);

///////////////////////////////////////////////////////////////////////////////
//	CollectVolumes
///////////////////////////////////////////////////////////////////////////////

//	CollectVolumes
inline void CollectAssociated(std::vector<Volume*>& vVolumesOut, Grid& grid,
                              GridObject* obj, bool clearContainer = true);

///	Collects all volumes that exist in the given grid which contain the given vertex.
UG_API 
void CollectVolumes(std::vector<Volume*>& vVolumesOut, Grid& grid, Vertex* vrt, bool clearContainer = true);

inline void CollectAssociated(std::vector<Volume*>& vVolumesOut,
					Grid& grid, Vertex* vrt, bool clearContainer = true);

///	Collects all volumes that exist in the given grid which contain the given edge.
/** if EDGEOPT_STORE_ASSOCIATED_VOLUMES is enabled then the algorithm uses
 * associated_volumes_begin(e), associated_volumes_end(e)
 * in order to find the volumes.
 * if not the algorithm iterates over all volumes associated with one of the edges
 * end-points and returns each volume that contains the edge.
 * \{
 */
UG_API 
void CollectVolumes(std::vector<Volume*>& vVolumesOut, Grid& grid, Edge* e, bool clearContainer = true);

inline void CollectAssociated(std::vector<Volume*>& vVolumesOut,
					Grid& grid, Edge* e, bool clearContainer = true);
/** \} */

///	Collects all volumes that exist in the given grid which contain the given face.
/** if FACEOPT_STORE_ASSOCIATED_VOLUMES is enabled then the algorithm uses
 * associated_volumes_begin(f), associated_volumes_end(f)
 * in order to find the volumes.
 * if not the algorithm iterates over all volumes associated with one of the faces
 * end-points and returns each volume that contains the face.
 * \{
 */
void CollectVolumes(std::vector<Volume*>& vVolumesOut, Grid& grid, Face* f, bool clearContainer = true, bool ignoreAssociatedVolumes = false);
inline void CollectAssociated(std::vector<Volume*>& vVolumesOut,
					Grid& grid, Face* f, bool clearContainer = true,
					bool ignoreAssociatedVolumes = false);

UG_API 
void CollectVolumes(std::vector<Volume*>& vVolumesOut, Grid& grid, FaceDescriptor& fd, bool clearContainer = true);
inline void CollectAssociated(std::vector<Volume*>& vVolumesOut,
					Grid& grid, FaceDescriptor& fd, bool clearContainer = true);
/// \}

///	Collects all volumes. (Returns the volume itself)
UG_API 
void CollectVolumes(std::vector<Volume*>& vVolumesOut, Grid& grid, Volume* v, bool clearContainer = true);

///	Collects all volumes. (Returns the volume itself)
inline void CollectAssociated(std::vector<Volume*>& vVolumesOut,
					Grid& grid, Volume* vol, bool clearContainer = true);


////////////////////////////////////////////////////////////////////////
//	VolumeContains
/// \ingroup lib_grid_algorithms_volume_util
///	returns true if the given volume contains the given vertex
bool VolumeContains(VolumeVertices* v, Vertex* vrt);

/// \ingroup lib_grid_algorithms_volume_util
///	returns true if the given volume contains the given edge
UG_API 
bool VolumeContains(Volume* v, EdgeVertices* ev);

/// \ingroup lib_grid_algorithms_volume_util
///	returns true if the given volume contains the given face
UG_API 
bool VolumeContains(Volume* v, FaceVertices* f);

// end of group lib_grid_algorithms_neighborhood_util
///	\}
}//	end of namespace libGrid

//	include template implementation
#include "grid_util_impl.hpp"

#endif
