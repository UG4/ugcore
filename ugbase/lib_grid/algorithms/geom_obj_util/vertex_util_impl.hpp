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
	FaceIterator iterEnd = grid.associated_faces_end(vrt);
	for(FaceIterator iter = grid.associated_faces_begin(vrt);
		iter != iterEnd; iter++)
	{
		vector3 vN;
		CalculateNormal(vN, *iter, aaPos);
		VecAdd(nOut, nOut, vN);
	}

	VecNormalize(nOut, nOut);
}

}//	end of namespace

#endif
