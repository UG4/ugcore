/*
 * Copyright (c) 2015:  G-CSC, Goethe University Frankfurt
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

#endif