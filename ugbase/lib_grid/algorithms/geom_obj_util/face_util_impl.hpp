//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m02 d05

#ifndef __H__LIB_GRID__FACE_UTIL_IMPL__
#define __H__LIB_GRID__FACE_UTIL_IMPL__

#include "face_util.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
template <class TIterator>
number AreaFaceQuality(TIterator facesBegin, TIterator facesEnd,
					   Grid::VertexAttachmentAccessor<APosition>& aaPos)
{
//	if the area is empty return 0 (bad)
	if(facesBegin == facesEnd)
		return 0;

//	get the first
	number q = FaceQuality(*facesBegin, aaPos);
	++facesBegin;

//	iterate over the others and find a worse one
	for(; facesBegin != facesEnd; ++facesBegin){
		number tq = FaceQuality(*facesBegin, aaPos);
		if(tq < q)
			q = tq;
	}

//	return the quality
	return q;
}

////////////////////////////////////////////////////////////////////////
inline void Triangulate(Grid& grid,
						Grid::VertexAttachmentAccessor<APosition>* paaPos)
{
	Triangulate(grid, grid.begin<Quadrilateral>(),
				grid.end<Quadrilateral>(), paaPos);
}

////////////////////////////////////////////////////////////////////////
//	CalculateFaceCenter
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

////////////////////////////////////////////////////////////////////////
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

////////////////////////////////////////////////////////////////////////
//	project points to surface 
template <class TTriangleIterator, class TAAPosVRT>
bool ProjectPointToSurface(vector3& vOut, const vector3& v, const vector3& n,
						   TTriangleIterator trisBegin, TTriangleIterator trisEnd,
						   TAAPosVRT& aaPos, bool compareNormals)
{
	vector3 vInter;
	bool gotOne = false;
	number b1, b2, t, tBest;

//	iterate through all triangles and find the closest intersection
	for(TTriangleIterator iter = trisBegin; iter != trisEnd; ++iter)
	{
		Triangle* tri = *iter;
		if(RayTriangleIntersection(vInter, b1, b2, t, aaPos[tri->vertex(0)],
								aaPos[tri->vertex(1)], aaPos[tri->vertex(2)], v, n))
		{
		//	check the normal
			vector3 tn = n;
			if(compareNormals){
				CalculateTriangleNormal(tn, aaPos[tri->vertex(0)], aaPos[tri->vertex(1)],
										aaPos[tri->vertex(2)]);
			}
		//	the triangle normal and the point - normal should match
			if(VecDot(tn, n) > 0){
				if(gotOne){
					if(fabs(t) < tBest){
						vOut = vInter;
						tBest = fabs(t);
					}
				}
				else{
					gotOne = true;
					vOut = vInter;
					tBest = fabs(t);
				}
			}
		}
	}

	return gotOne;
}

}//	end of namespace

#endif
