// created by Sebastian Reiter
// s.b.reiter@gmail.com

#ifndef __H__UG_pyramid_rules_impl
#define __H__UG_pyramid_rules_impl

#include "pyramid_rules.h"
#include "grid_object_ids.h"

namespace ug{
namespace pyra_rules{

template <class TCmp>
int ConvertToTetrahedra(int* newIndsOut, TCmp cmp)
{
//	The idea and implementation follows:
//	Dompierre et al., "How to Subdivide Pyramids, Prisms and Hexahedra into Tetrahedra"

//	find the smallest index of the base quadrilateral
	int smallest = 0;
	for(int i = 1; i < 4; ++i){
		if(cmp(i, smallest))
			smallest = i;
	}

	if((smallest == 0) || (smallest == 2)){
		int i = 0;
		newIndsOut[i++] = GOID_TETRAHEDRON;
		newIndsOut[i++] = 0; newIndsOut[i++] = 1;
		newIndsOut[i++] = 2; newIndsOut[i++] = TOP_VERTEX;

		newIndsOut[i++] = GOID_TETRAHEDRON;
		newIndsOut[i++] = 0; newIndsOut[i++] = 2;
		newIndsOut[i++] = 3; newIndsOut[i++] = TOP_VERTEX;
		return i;
	}
	else{
		int i = 0;
		newIndsOut[i++] = GOID_TETRAHEDRON;
		newIndsOut[i++] = 0; newIndsOut[i++] = 1;
		newIndsOut[i++] = 3; newIndsOut[i++] = TOP_VERTEX;

		newIndsOut[i++] = GOID_TETRAHEDRON;
		newIndsOut[i++] = 1; newIndsOut[i++] = 2;
		newIndsOut[i++] = 3; newIndsOut[i++] = TOP_VERTEX;
		return i;
	}
}

}//	end of namespace
}//	end of namespace

#endif	//__H__UG_pyramid_rules_impl
