//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m01 d15

#ifndef __H__LIB_GRID__FACE_UTIL__
#define __H__LIB_GRID__FACE_UTIL__

#include "lib_grid/lg_base.h"

namespace ug
{

/** \defgroup faceUtil Face Util
 * @{
 */


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
 */
void CalculateNormal(vector3& vNormOut, Face* face,
					Grid::VertexAttachmentAccessor<APosition>& aaPos);

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
 *	to exactly one volume.*/
bool IsVolumeBoundaryFace(Grid& grid, Face* f);

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


////////////////////////////////////////////////////////////////////////
//	GetNeighbours
///	collects neighbours of the given side of a face.
/**
 *	collects all triangles that are adjacent to the given side of f.
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
/**
 * swaps orientation of faces so that all neighboured
 * faces share the same.
 */
void FixOrientation(Grid& grid, FaceIterator facesBegin, FaceIterator facesEnd);
	
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
CalculateFaceCenter(Face* f, TVertexPositionAttachmentAccessor& aaPosVRT)
{
	uint numVrts = f->num_vertices();
	typename TVertexPositionAttachmentAccessor::ValueType v;
//	init v with 0.
	VecSet(v, 0);

//	sum up
	for(uint i = 0; i < numVrts; ++i)
	{
		VecAdd(v, v, aaPosVRT[f->vertex(i)]);
	}

//	average
	if(numVrts > 0)
		VecScale(v, v, 1./(number)numVrts);

	return v;
}

template<class TVertexPositionAttachmentAccessor>
typename TVertexPositionAttachmentAccessor::ValueType
CalculateCenter(Face* f, TVertexPositionAttachmentAccessor& aaPosVRT)
{
	uint numVrts = f->num_vertices();
	typename TVertexPositionAttachmentAccessor::ValueType v;
//	init v with 0.
	VecSet(v, 0);

//	sum up
	for(uint i = 0; i < numVrts; ++i)
	{
		VecAdd(v, v, aaPosVRT[f->vertex(i)]);
	}

//	average
	if(numVrts > 0)
		VecScale(v, v, 1./(number)numVrts);

	return v;
}


////////////////////////////////////////////////////////////////////////
//	FindFaceByCoordinate
///	returns the face whose center is closest to coord.
/**
 * This method does not necessarily return the face that contains the given coordinate.
 * Instrad it will simply search for the face whose center is closest to the specified
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
							TVertexPositionAttachmentAccessor& aaPosVRT)
{
	if(iterBegin == iterEnd)
		return NULL;

	FaceIterator iter = iterBegin;
	Face* bestFace = *iter;
	number bestDistSq = VecDistanceSq(coord, CalculateFaceCenter(bestFace, aaPosVRT));
	iter++;

	while(iter != iterEnd)
	{
		number distSq = VecDistanceSq(coord, CalculateFaceCenter(*iter, aaPosVRT));
		if(distSq < bestDistSq)
		{
			bestDistSq = distSq;
			bestFace = *iter;
		}

		++iter;
	}

	return bestFace;
}

/**@}*/ // end of doxygen defgroup command

}//	end of namespace

#endif
