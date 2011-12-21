//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m01 d15

#ifndef __H__LIB_GRID__FACE_UTIL__
#define __H__LIB_GRID__FACE_UTIL__

#include "lib_grid/lg_base.h"

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
int GetFaceIndex(Volume* vol, Face* f);


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
void CalculateNormal(vector3& vNormOut, FaceVertices* face,
					Grid::AttachmentAccessor<VertexBase, APosition>& aaPos);

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
void CalculateNormalNoNormalize(vector3& vNormOut, FaceVertices* face,
						Grid::AttachmentAccessor<VertexBase, APosition>& aaPos);
					
////////////////////////////////////////////////////////////////////////
//	CalculateFaceNormals
///	calculates the normal of each face. Presumes that all faces are flat.
/**
 * aPos has to be attached to the vertices of the grid.
 * aPos should contain the position data.
 * Normals will be written to aNorm (face attachment).
 */
void CalculateFaceNormals(Grid& grid, const FaceIterator& facesBegin,
						const FaceIterator& facesEnd,
						AVector3& aPos, AVector3& aNorm);

////////////////////////////////////////////////////////////////////////
//	IsVolumeBoundaryFace
///	returns true if the given face is a boundary face.
/**	A face is regarded as a boundary face if it is adjacent
 *	to exactly one volume.
 *
 *	Please note that overloads of this function for Constrained- and
 *	ConstrainingFaces exist.*/
bool IsVolumeBoundaryFace(Grid& grid, Face* f);

////////////////////////////////////////////////////////////////////////
//	IsVolumeBoundaryFace
///	returns true if the given face is a boundary face.
/**	Overload for ConstrainedFace.
 *	A face is regarded as a boundary face if it is adjacent
 *	to exactly one volume.*/
bool IsVolumeBoundaryFace(Grid& grid, ConstrainedFace* f);

////////////////////////////////////////////////////////////////////////
//	IsVolumeBoundaryFace
///	returns true if the given face is a boundary face.
/**	Overload for ConstrainedFace.
 *	A face is regarded as a boundary face if it is adjacent
 *	to exactly one volume.*/
bool IsVolumeBoundaryFace(Grid& grid, ConstrainingFace* f);

////////////////////////////////////////////////////////////////////////
///	A wrapper for IsVolumeBoundaryFace
template <class TFace>
bool IsBoundaryFace3D(Grid& grid, TFace* f)
{return IsVolumeBoundaryFace(grid, f);}

////////////////////////////////////////////////////////////////////////
///	A wrapper for IsVolumeBoundaryFace
template <class TFace>
bool LiesOnBoundary(Grid& grid, TFace* f)
{return IsBoundaryFace3D(grid, f);}

////////////////////////////////////////////////////////////////////////
//	FaceArea
///	Returns the area of a convex face
template <class TAAPosVRT>
number FaceArea(Face* f, TAAPosVRT& aaPos);

////////////////////////////////////////////////////////////////////////
///	Returns the face with the smallest area
/**
 * Make sure that TIterator::value_type equals Face* and that aaPos operates
 * on the grid from which the faces were taken.
 */
template <class TIterator, class TAAPosVRT>
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
number FaceQuality(Face* f,
				Grid::VertexAttachmentAccessor<APosition> aaPos);

////////////////////////////////////////////////////////////////////////
//	AreaFaceQuality
///	returns a value between 0 (bad) and 1 (good) that describes the quality of the area. 
/**
 * returns the worst FaceQuality of the faces between facesBegin and FacesEnd.
 * TIterator has to be an iterator with a value-type compatible to Face*.
 */
template <class TIterator>
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
number TriangleQuality(vector3& v1,
						vector3& v2,
						vector3& v3);

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
void Triangulate(Grid& grid,
				QuadrilateralIterator iterBegin,
				QuadrilateralIterator iterEnd,
				Grid::VertexAttachmentAccessor<APosition>* paaPos = NULL);

inline void Triangulate(Grid& grid,
						Grid::VertexAttachmentAccessor<APosition>* paaPos = NULL);

////////////////////////////////////////////////////////////////////////
//	GetNeighbours
///	collects neighbours of the given side of a face.
/**
 *	collects all faces that are adjacent to the given side of f.
 */
void GetNeighbours(std::vector<Face*>& vFacesOut, Grid& grid, Face* f,
					int side, bool clearContainer = true);
					
////////////////////////////////////////////////////////////////////////
//	EdgeOrientationMatches
///	checks if the edge-orientation of the edge and the face matches.
/**
 * the match is positive if the face contains the vertices of ed
 * in the same order as ed.
 * please note: if the edge is contained by two faces and both
 * faces have the same edge-orientation as ed, then the face-orientation
 * of the faces differ.
 */
bool EdgeOrientationMatches(EdgeDescriptor& ed, Face* f);

////////////////////////////////////////////////////////////////////////
//	FixOrientation
///	creates uniform orientation of neighboured faces.
/** This algorithm uses Grid::mark
 *
 * swaps orientation of faces so that all neighboured
 * faces share the same.
 *
 * Value type of TFaceIterator has to be compatible with Face*.
 *
 * Note that all faces between faceBegin and facesEnd have to be members
 * of the specified grid.
 */
template <class TFaceIterator>
void FixFaceOrientation(Grid& grid, TFaceIterator facesBegin,
						TFaceIterator facesEnd);

////////////////////////////////////////////////////////////////////////
//	InvertOrientation
///	inverts the orientation of all faces between facesBegin and facesEnd
/**	Make sure that TFaceIterator::value_type is castable to Face*.
 *	TFaceIterator has to be compatible with stl-iterators.*/
template <class TFaceIterator>
void InvertOrientation(Grid& grid, TFaceIterator facesBegin,
					   TFaceIterator facesEnd);

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
 */
template<class TVertexPositionAttachmentAccessor>
typename TVertexPositionAttachmentAccessor::ValueType
CalculateFaceCenter(Face* f, TVertexPositionAttachmentAccessor& aaPosVRT);

////////////////////////////////////////////////////////////////////////
template<class TVertexPositionAttachmentAccessor>
typename TVertexPositionAttachmentAccessor::ValueType
CalculateCenter(FaceVertices* f, TVertexPositionAttachmentAccessor& aaPosVRT);


////////////////////////////////////////////////////////////////////////
//	FindFaceByCoordinate
///	returns the face whose center is closest to coord.
/**
 * This method does not necessarily return the face that contains the given coordinate.
 * Instead it will simply search for the face whose center is closest to the specified
 * coordinate.
 * TVertexPositionAttachmentAccessor has to be an AttachmentAccessor,
 * where AttachmentAccessor::ValueType is a vector-type compatible to
 * the lgmath vector descriptor.
 * The Accessor has to access an attachment of the vertices,
 * to which the faces between iterBegin and iterEnd refer.
 */
template<class TVertexPositionAttachmentAccessor>
Face* FindFaceByCoordinate(const typename TVertexPositionAttachmentAccessor::ValueType& coord,
							FaceIterator iterBegin, FaceIterator iterEnd,
							TVertexPositionAttachmentAccessor& aaPosVRT);

////////////////////////////////////////////////////////////////////////
//	project points to surface 
template <class TTriangleIterator, class TAAPosVRT>
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
int PointFaceTest(vector3& v, Face* f, TAAPosVRT& aaPos);

////////////////////////////////////////////////////////////////////////
///	returns true if the given face is degenerated.
/**	Faces are degenerated if at least one edge is shorter than the given threshold.*/
template <class TAAPosVRT>
bool IsDegenerated(Face* f, TAAPosVRT& aaPos, number threshold = SMALL);

/// @}

}//	end of namespace

////////////////////////////////
//	include implementation
#include "face_util_impl.hpp"

#endif
