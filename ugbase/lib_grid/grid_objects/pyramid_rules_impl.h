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

#ifndef __H__UG_pyramid_rules_impl
#define __H__UG_pyramid_rules_impl

#include "pyramid_rules.h"
#include "grid_object_ids.h"

namespace ug{
namespace pyra_rules{

template <typename TCmp>
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
		newIndsOut[i++] = GridObjectID::GOID_TETRAHEDRON;
		newIndsOut[i++] = 0; newIndsOut[i++] = 1;
		newIndsOut[i++] = 2; newIndsOut[i++] = TOP_VERTEX;

		newIndsOut[i++] = GridObjectID::GOID_TETRAHEDRON;
		newIndsOut[i++] = 0; newIndsOut[i++] = 2;
		newIndsOut[i++] = 3; newIndsOut[i++] = TOP_VERTEX;
		return i;
	}
	else{
		int i = 0;
		newIndsOut[i++] = GridObjectID::GOID_TETRAHEDRON;
		newIndsOut[i++] = 0; newIndsOut[i++] = 1;
		newIndsOut[i++] = 3; newIndsOut[i++] = TOP_VERTEX;

		newIndsOut[i++] = GridObjectID::GOID_TETRAHEDRON;
		newIndsOut[i++] = 1; newIndsOut[i++] = 2;
		newIndsOut[i++] = 3; newIndsOut[i++] = TOP_VERTEX;
		return i;
	}
}

}//	end of namespace
}//	end of namespace

#endif