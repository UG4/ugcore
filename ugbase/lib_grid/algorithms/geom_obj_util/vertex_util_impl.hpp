//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m02 d05

#ifndef __H__LIB_GRID__VERTEX_UTIL_IMPL__
#define __H__LIB_GRID__VERTEX_UTIL_IMPL__

#include "vertex_util.h"
#include "face_util.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
template <class TAAPosVRT>
void CalculateVertexNormal(vector3& nOut, Grid& grid, VertexBase* vrt, TAAPosVRT& aaPos)
{
//	set all normal to zero
	nOut = vector3(0, 0, 0);

//	loop through all associated faces, calculate their normal and add them to thee normal
	Grid::AssociatedFaceIterator iterEnd = grid.associated_faces_end(vrt);
	for(Grid::AssociatedFaceIterator iter = grid.associated_faces_begin(vrt);
		iter != iterEnd; iter++)
	{
		vector3 vN;
		CalculateNormal(vN, *iter, aaPos);
		VecAdd(nOut, nOut, vN);
	}

	VecNormalize(nOut, nOut);
}

////////////////////////////////////////////////////////////////////////
template <class TIterator, class AAPosVRT>
void LaplacianSmooth(Grid& grid, TIterator vrtsBegin,
					TIterator vrtsEnd, AAPosVRT& aaPos,
					number alpha, int numIterations)
{
	for(int iteration = 0; iteration < numIterations; ++iteration){
	//	iterate through all vertices
		for(TIterator iter = vrtsBegin; iter != vrtsEnd; ++iter){
		//	smooth each one
			VertexBase* vrt = *iter;
			vector3 v(0, 0, 0);
			int num = 0;
			
			Grid::AssociatedEdgeIterator edgesEnd = grid.associated_edges_end(vrt);
			for(Grid::AssociatedEdgeIterator eIter = grid.associated_edges_begin(vrt);
				eIter != edgesEnd; ++eIter)
			{
				VecAdd(v, v, aaPos[GetConnectedVertex(*eIter, vrt)]);
				++num;
			}
			
			if(num > 0){
				VecScale(v, v, 1. / (number)num);
				VecSubtract(v, v, aaPos[vrt]);
				VecScale(v, v, alpha);
				VecAdd(aaPos[vrt], aaPos[vrt], v);
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////
template<class TVertexPositionAttachmentAccessor>
inline
typename TVertexPositionAttachmentAccessor::ValueType
CalculateCenter(VertexBase* v, TVertexPositionAttachmentAccessor& aaPosVRT)
{
	return aaPosVRT[v];
}

}//	end of namespace

#endif
