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

#ifndef __H__LIB_GRID__FACE_UTIL__
#define __H__LIB_GRID__FACE_UTIL__

#include "lib_grid/grid/grid.h"
#include "lib_grid/grid_objects/grid_objects.h"
#include "lib_grid/common_attachments.h"

namespace ug
{

/**
 * \brief contains methods to manipulate faces
 * 
 * \defgroup lib_grid_algorithms_face_util face util
 * \ingroup lib_grid_algorithms
 * @{
 */


////////////////////////////////////////////////////////////////////////
//	GetFaceIndex
///	returns the index at which face f is found in the given object
/**
 * returns -1 if the face was not found.
 */
UG_API int GetFaceIndex(Volume* vol, Face* f);


////////////////////////////////////////////////////////////////////////
//	CalculateNormal
///	calculates the normal of the given face
/**
 * aaPos has to be a valid positition-attachment-accessor to the vertices
 * of the grid which contains the face.
 * if the face contains less than 3 vertices (0, 0, 0) will be written
 * to vNormOut (this should never happen!).
 * For triangles the normal is calculated using the standard cross-product.
 * for quadrilaterals the normals of the two sub-triangles (0, 1, 2) and
 * (2, 3, 0) is calculated and averaged.
 * If the face contains more than 4 vertices the normal of the first
 * sub-triangle is returned.
 *
 * Performs normalization on the calcluated normals.
 */
UG_API 
void CalculateNormal(vector3& vNormOut, const FaceVertices* face,
					Grid::AttachmentAccessor<Vertex, APosition>& aaPos);

////////////////////////////////////////////////////////////////////////
//	CalculateNormal
///	calculates the normal of the given face
/**
 * aaPos has to be a valid positition-attachment-accessor to the vertices
 * of the grid which contains the face.
 * if the face contains less than 3 vertices (0, 0, 0) will be written
 * to vNormOut (this should never happen!).
 * For triangles the normal is calculated using the standard cross-product.
 * for quadrilaterals the normals of the two sub-triangles (0, 1, 2) and
 * (2, 3, 0) is calculated and averaged.
 * If the face contains more than 4 vertices the normal of the first
 * sub-triangle is returned.
 *
 * performs no normalization on the calculated normals
 */
UG_API 
void CalculateNormalNoNormalize(vector3& vNormOut, FaceVertices* face,
						Grid::AttachmentAccessor<Vertex, APosition>& aaPos);
					
////////////////////////////////////////////////////////////////////////
//	CalculateFaceNormals
///	calculates the normal of each face. Presumes that all faces are flat.
/**
 * aPos has to be attached to the vertices of the grid.
 * aPos should contain the position data.
 * Normals will be written to aNorm (face attachment).
 */
UG_API 
void CalculateFaceNormals(Grid& grid, const FaceIterator& facesBegin,
						const FaceIterator& facesEnd,
						AVector3& aPos, AVector3& aNorm);

////////////////////////////////////////////////////////////////////////
/// returns the number of associated volumes of the specified face
int NumAssociatedVolumes(Grid& grid, Face* f);

////////////////////////////////////////////////////////////////////////
//	IsVolumeBoundaryFace
///	returns true if the given face is a boundary face.
/**	A face is regarded as a boundary face if it is adjacent
 *	to exactly one volume.
 *
 *	Please note that overloads of this function for Constrained- and
 *	ConstrainingFaces exist.*/
UG_API 
bool IsVolumeBoundaryFace(Grid& grid, Face* f);

////////////////////////////////////////////////////////////////////////
//	IsVolumeBoundaryFace
///	returns true if the given face is a boundary face.
/**	Overload for ConstrainedFace.
 *	A face is regarded as a boundary face if it is adjacent
 *	to exactly one volume.*/
UG_API 
bool IsVolumeBoundaryFace(Grid& grid, ConstrainedFace* f);

////////////////////////////////////////////////////////////////////////
//	IsVolumeBoundaryFace
///	returns true if the given face is a boundary face.
/**	Overload for ConstrainedFace.
 *	A face is regarded as a boundary face if it is adjacent
 *	to exactly one volume.*/
UG_API 
bool IsVolumeBoundaryFace(Grid& grid, ConstrainingFace* f);

////////////////////////////////////////////////////////////////////////
///	A wrapper for IsVolumeBoundaryFace
template <class TFace>
UG_API 
bool IsBoundaryFace3D(Grid& grid, TFace* f)
{return IsVolumeBoundaryFace(grid, f);}

////////////////////////////////////////////////////////////////////////
///	A wrapper for IsVolumeBoundaryFace
template <class TFace>
UG_API 
bool LiesOnBoundary(Grid& grid, TFace* f)
{return IsBoundaryFace3D(grid, f);}

////////////////////////////////////////////////////////////////////////
//	FaceArea
///	Returns the area of a convex face
template <class TAAPosVRT>
UG_API 
number FaceArea(FaceVertices* f, TAAPosVRT& aaPos);

////////////////////////////////////////////////////////////////////////
//  FaceArea
///	Returns the area sum of convex faces
template <class TIterator, class TAAPosVRT>
UG_API
number FaceArea(TIterator facesBegin, TIterator facesEnd, TAAPosVRT& aaPos);

////////////////////////////////////////////////////////////////////////
///	Returns the face with the smallest area
/**
 * Make sure that TIterator::value_type equals Face* and that aaPos operates
 * on the grid from which the faces were taken.
 */
template <class TIterator, class TAAPosVRT>
UG_API 
Face* FindSmallestFace(TIterator facesBegin, TIterator facesEnd, TAAPosVRT& aaPos);

////////////////////////////////////////////////////////////////////////
//	FaceQuality
///	a simple measure for the quality of a face
/**
 * returns a value between 0 and 1, where 1 indicates a
 * good quality and 0 a bad.
 * This method checks the dot-products of edges at the corners.
 * The worst one (closest to 1 or -1) determines the quality.
 * \sa TriangleQuality
 */
UG_API 
number FaceQuality(FaceVertices* f, Grid::VertexAttachmentAccessor<APosition> aaPos);

////////////////////////////////////////////////////////////////////////
//	AreaFaceQuality
///	returns a value between 0 (bad) and 1 (good) that describes the quality of the area. 
/**
 * returns the worst FaceQuality of the faces between facesBegin and FacesEnd.
 * TIterator has to be an iterator with a value-type compatible to Face*.
 */
template <class TIterator>
UG_API 
number AreaFaceQuality(TIterator facesBegin, TIterator facesEnd,
					   Grid::VertexAttachmentAccessor<APosition>& aaPos);

////////////////////////////////////////////////////////////////////////
//	TriangleQuality
///	a simple measure for the quality of a triangle
/**
 * returns a value between 0 and 1, where 1 indicates a
 * good quality and 0 a bad.
 * This method checks the dot-products of edges at the corners.
 * The worst one (closest to 1 or -1) determines the quality.
 * \sa FaceQuality
 */
UG_API number TriangleQuality(vector3& v1, vector3& v2, vector3& v3);

////////////////////////////////////////////////////////////////////////
//	Triangulate
///	removes the quadrilateral and replaces it by two triangles.
/**
 * if paaPos is set to NULL, the quadrilateral will be splitted
 * along the edge between the first and the third vertex.
 * If paaPos points to a position-attachment-accessor,
 * then the new edge will be chosen so that the worst triangle-quality
 * is better.
 */
UG_API 
void Triangulate(Grid& grid, Quadrilateral* q,
				Grid::VertexAttachmentAccessor<APosition>* paaPos = NULL);

////////////////////////////////////////////////////////////////////////
//	Triangulate
///	replaces all specified quadrilaterals by triangles.
/**
 * if paaPos is set to NULL, the quadrilaterals will be splitted
 * along the edge between their first and their third vertex.
 * If paaPos points to a position-attachment-accessor,
 * then the new edge will be chosen so that the worst triangle-quality
 * is better.
 */
UG_API 
void Triangulate(Grid& grid,
				QuadrilateralIterator iterBegin,
				QuadrilateralIterator iterEnd,
				Grid::VertexAttachmentAccessor<APosition>* paaPos = NULL);

UG_API 
inline void Triangulate(Grid& grid,
						Grid::VertexAttachmentAccessor<APosition>* paaPos = NULL);

////////////////////////////////////////////////////////////////////////
//	GetNeighbours
///	collects neighbours of the given side of a face.
/**
 *	collects all faces that are adjacent to the given side of f.
 */
UG_API 
void GetNeighbours(std::vector<Face*>& vFacesOut, Grid& grid, Face* f,
					int side, bool clearContainer = true);
					

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	template methods
////////////////////////////////////////////////////////////////////////
//	CalculateFaceCenter
///	calculates the center of a face.
/**
 * TVertexPositionAttachmentAccessor has to be an AttachmentAccessor,
 * where AttachmentAccessor::ValueType is a vector-type compatible to
 * the lgmath vector descriptor.
 * The accessor has to access an attachment of the vertices,
 * to which f refers.
 *
 * \{
 */
template<class TVertexPositionAttachmentAccessor>
UG_API 
typename TVertexPositionAttachmentAccessor::ValueType
CalculateFaceCenter(const Face* f, TVertexPositionAttachmentAccessor& aaPosVRT);

////////////////////////////////////////////////////////////////////////
template<class TVertexPositionAttachmentAccessor>
UG_API 
typename TVertexPositionAttachmentAccessor::ValueType
CalculateCenter(const FaceVertices* f, TVertexPositionAttachmentAccessor& aaPosVRT);

/** \} */


////////////////////////////////////////////////////////////////////////
///	returns the weighted center of the vertices of the given face
/** TAAWeightVRT has to be an attachment to the vertices of the grid in which
 * f is contained, with ValueType number (or compatible).
 */
template<class TAAPosVRT, class TAAWeightVRT>
UG_API
typename TAAPosVRT::ValueType
CalculateCenter(const FaceVertices* f, TAAPosVRT& aaPos, TAAWeightVRT& aaWeight);


////////////////////////////////////////////////////////////////////////
///	Returns true if the given point lies inside the given face.
/**	\note	The method only works properly, if the point and the face are located
 * 			in the same x-y-plane.
 */
template <class vector_t, class TAAPos>
UG_API bool
ContainsPoint(const FaceVertices* f, const vector_t& p, TAAPos aaPos);


////////////////////////////////////////////////////////////////////////
//	project points to surface 
template <class TTriangleIterator, class TAAPosVRT>
UG_API 
bool ProjectPointToSurface(vector3& vOut, const vector3& v, const vector3& n,
						   TTriangleIterator trisBegin, TTriangleIterator trisEnd,
						   TAAPosVRT& aaPos, bool compareNormals = false);

////////////////////////////////////////////////////////////////////////
/**
 * returns 1 if a point lies in front of a face,
 * 0 if it lies on the face and -1 if it lies behind the face.
 * TAAPosVRT has to be an AttachmentAccessor compatible type that
 * operates on vector3.
 */
template <class TAAPosVRT>
UG_API 
int PointFaceTest(vector3& v, Face* f, TAAPosVRT& aaPos);

////////////////////////////////////////////////////////////////////////
///	returns true if the given face is degenerated.
/**	Faces are degenerated if at least one edge is shorter than the given threshold.*/
template <class TAAPosVRT>
UG_API 
bool IsDegenerated(Face* f, TAAPosVRT& aaPos, number threshold = SMALL);


////////////////////////////////////////////////////////////////////////
///	Refines the face by connecting its sides with the new center.
/**	Make sure that the specified vertex is belongs to the specified grid, and
 * that its position lies inside the specified face (self-intersections would
 * occur if it would lie outside).
 * The original face may optionally be deleted.
 */
UG_API
void InsertCenterVertex(Grid& g, Face* f, Vertex* vrt, bool eraseOldFace);

/// @}

}//	end of namespace

////////////////////////////////
//	include implementation
#include "face_util_impl.hpp"

#endif
