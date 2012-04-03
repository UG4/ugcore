//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m02 d19

#ifndef __H__LIB_GRID__VOLUME_UTIL_IMPL__
#define __H__LIB_GRID__VOLUME_UTIL_IMPL__

#include "lib_grid/lg_base.h"

namespace ug
{
////////////////////////////////////////////////////////////////////////
//	PointIsInsideTetrahedron
inline bool
PointIsInsideTetrahedron(const vector3& v, Tetrahedron* tet,
						 Grid::VertexAttachmentAccessor<APosition>& aaPos)
{
	return PointIsInsideTetrahedron(v, aaPos[tet->vertex(0)], aaPos[tet->vertex(1)],
									aaPos[tet->vertex(2)], aaPos[tet->vertex(3)]);
}

////////////////////////////////////////////////////////////////////////
//	CalculateCenter
template<class TVertexPositionAttachmentAccessor>
typename TVertexPositionAttachmentAccessor::ValueType
CalculateCenter(const VolumeVertices* vol, TVertexPositionAttachmentAccessor& aaPosVRT)
{
	uint numVrts = vol->num_vertices();
	typename TVertexPositionAttachmentAccessor::ValueType v;
//	init v with 0.
	VecSet(v, 0);

//	sum up
	for(uint i = 0; i < numVrts; ++i)
	{
		VecAdd(v, v, aaPosVRT[vol->vertex(i)]);
	}

//	average
	if(numVrts > 0)
		VecScale(v, v, 1./(number)numVrts);

	return v;
}

////////////////////////////////////////////////////////////////////////
//	CheckOrientation
template<class TAAPosVRT>
bool
CheckOrientation(Volume* vol, TAAPosVRT& aaPosVRT)
{
//	some typedefs
	typedef typename TAAPosVRT::ValueType vector_t;
	
//	First calculate the center of the volume
	vector_t volCenter = CalculateCenter(vol, aaPosVRT);
	
//	now check for each side whether it points away from the center.
	size_t numFaces = vol->num_faces();
	FaceDescriptor fd;
	vector_t normal;
	for(size_t i = 0; i < numFaces; ++i){
		vol->face_desc(i, fd);
		CalculateNormal(normal, &fd, aaPosVRT);
		
	//	in order to best approximate quadrilateral faces, we'll calculate the
	//	center of the face and compare that to the volCenter.
		vector_t faceCenter = CalculateCenter(&fd, aaPosVRT);
		
	//	now compare normal and center
		vector_t dir;
		VecSubtract(dir, faceCenter, volCenter);
		if(VecDot(dir, normal) < 0)
			return false;
	}
	
//	all center / normal checks succeeded. Orientation is fine.
	return true;
}

template<class TVolIterator, class TAAPosVRT>
int
FixOrientation(Grid& grid, TVolIterator volsBegin, TVolIterator volsEnd,
			   TAAPosVRT& aaPosVRT)
{
	int numFlips = 0;
//	iterate through all volumes
	for(VolumeIterator iter = volsBegin; iter != volsEnd; ++iter){
	//	check whether the orientation is fine
		if(!CheckOrientation(*iter, aaPosVRT)){
			grid.flip_orientation(*iter);
			++numFlips;
		}
	}
	
	return numFlips;
}

}//	end of namespace

#endif
