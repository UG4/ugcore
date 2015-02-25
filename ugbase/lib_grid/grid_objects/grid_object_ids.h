// created by Sebastian Reiter
// s.b.reiter@gmail.com

#ifndef __H__UG_grid_object_ids
#define __H__UG_grid_object_ids

namespace ug{

enum GridObjectID{
	GOID_VERTEX	= 0,
	GOID_EDGE = 1,
	GOID_TRIANGLE = 2,
	GOID_QUADRILATERAL = 3,
	GOID_TETRAHEDRON = 4,
	GOID_PYRAMID = 5,
	GOID_PRISM = 6,
	GOID_OCTAHEDRON = 7,
	GOID_HEXAHEDRON = 8
};

}//	end of namespace

#endif	//__H__UG_grid_object_ids
