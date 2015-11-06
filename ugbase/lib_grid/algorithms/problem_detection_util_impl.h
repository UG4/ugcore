/*
 * Copyright (c) 2014-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG_problem_detection_util_impl
#define __H__UG_problem_detection_util_impl

#include "problem_detection_util.h"
#include "lib_grid/grid/grid_base_objects.h"

namespace ug{

template <class TIter, class TAAPos>
size_t FindSlivers(std::vector<typename TIter::value_type>& sliversOut,
				 TIter elemsBegin, TIter elemsEnd, number thresholdRatio,
				 TAAPos aaPos, bool clearContainer)
{
	typedef typename TIter::value_type elem_ptr_t;

	if(clearContainer)
		sliversOut.clear();

	size_t numInitialSlivers = sliversOut.size();
	for(TIter iter = elemsBegin; iter != elemsEnd; ++iter){
		elem_ptr_t e = *iter;
		Volume::ConstVertexArray v = e->vertices();
		if((e->reference_object_id() == ROID_TETRAHEDRON) &&
			(IsSliver(aaPos[v[0]], aaPos[v[1]], aaPos[v[2]], aaPos[v[3]], thresholdRatio) != -1))
		{
			sliversOut.push_back(e);
		}
	}

	return sliversOut.size() - numInitialSlivers;
}

}//	end of namespace

#endif	//__H__problem_detection_util_impl
