// created by Sebastian Reiter
// s.b.reiter@gmail.com

#ifndef __H__UG_prism_rules_impl
#define __H__UG_prism_rules_impl

#include "prism_rules.h"
#include "grid_object_ids.h"

namespace ug{
namespace prism_rules{

template <class TCmp>
int ConvertToTetrahedra(int* newIndsOut, TCmp cmp)
{
//	The idea and implementation follows:
//	Dompierre et al., "How to Subdivide Pyramids, Prisms and Hexahedra into Tetrahedra"

//	find the smallest index
	int smallest = 0;
	for(int i = 1; i < 6; ++i){
		if(cmp(i, smallest))
			smallest = i;
	}

//	create a local ordering
	static const int permutations [6][6] = {{0, 1, 2, 3, 4, 5},
											{1, 2, 0, 4, 5, 3},
											{2, 0, 1, 5, 3, 4},
											{3, 5, 4, 0, 2, 1},
											{4, 3, 5, 1, 0, 2},
											{5, 4, 3, 2, 1, 0}};

	const int* I = permutations[smallest];

	int t0 = cmp(I[1], I[5])	?	I[1] : I[5];
	int t1 = cmp(I[2], I[4])	?	I[2] : I[4];

	if(cmp(t0, t1)){
		int i = 0;
		newIndsOut[i++] = GOID_TETRAHEDRON;
		newIndsOut[i++] = I[0]; newIndsOut[i++] = I[1];
		newIndsOut[i++] = I[2]; newIndsOut[i++] = I[5];

		newIndsOut[i++] = GOID_TETRAHEDRON;
		newIndsOut[i++] = I[0]; newIndsOut[i++] = I[1];
		newIndsOut[i++] = I[5]; newIndsOut[i++] = I[4];

		newIndsOut[i++] = GOID_TETRAHEDRON;
		newIndsOut[i++] = I[0]; newIndsOut[i++] = I[4];
		newIndsOut[i++] = I[5]; newIndsOut[i++] = I[3];
		return i;
	}
	else{
		int i = 0;
		newIndsOut[i++] = GOID_TETRAHEDRON;
		newIndsOut[i++] = I[0]; newIndsOut[i++] = I[1];
		newIndsOut[i++] = I[2]; newIndsOut[i++] = I[4];

		newIndsOut[i++] = GOID_TETRAHEDRON;
		newIndsOut[i++] = I[0]; newIndsOut[i++] = I[4];
		newIndsOut[i++] = I[2]; newIndsOut[i++] = I[5];

		newIndsOut[i++] = GOID_TETRAHEDRON;
		newIndsOut[i++] = I[0]; newIndsOut[i++] = I[4];
		newIndsOut[i++] = I[5]; newIndsOut[i++] = I[3];
		return i;
	}
}

}//	end of namespace	
}//	end of namespace

#endif	//__H__UG_prism_rules_impl
