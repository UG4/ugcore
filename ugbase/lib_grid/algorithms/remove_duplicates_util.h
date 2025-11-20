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

#ifndef __H__UG_remove_duplicates_util
#define __H__UG_remove_duplicates_util

#include "lib_grid/grid/grid.h"
#include "lib_grid/grid/grid_util.h"

namespace ug{

////////////////////////////////////////////////////////////////////////
///	Removes elements which share the same set of vertices with other elements in the given range
template <typename TElemIter>
void RemoveDuplicates(Grid& grid, TElemIter elemsBegin, TElemIter elemsEnd)
{
	using TElemPtr = typename TElemIter::value_type;
	using TElem = typename PtrToValueType<TElemPtr>::base_type;

	grid.begin_marking();

//	the first is the double, the second the original
	std::vector<std::pair<TElem*, TElem*> > doubles;
	typename Grid::traits<TElem>::secure_container	elems;

	for(TElemIter iter = elemsBegin; iter != elemsEnd; ++iter){
		TElem* e = *iter;
		if(!grid.is_marked(e)){
			typename TElem::ConstVertexArray vrts = e->vertices();
			const size_t numVrts = e->num_vertices();

			bool allMarked = true;
			for(size_t i = 0; i < numVrts; ++i){
				if(!grid.is_marked(vrts[i])){
					allMarked = false;
					break;
				}
			}

			bool isDouble = false;

			if(allMarked){
			//	a necessary condition is met. However not yet sufficient.
			//	find marked elems wich connect the same vertices.
				grid.associated_elements(elems, vrts[0]);
				for(size_t i = 0; i < elems.size(); ++i){
					TElem* te = elems[i];
					if(grid.is_marked(te) && CompareVertices(e, te)){
					//	e is a double
						isDouble = true;
						doubles.push_back(std::make_pair(e, te));
						break;
					}
				}
			}

		//	finally mark e and its vertices (every processed non-double element is marked).
			if(!isDouble){
				grid.mark(e);
				for(size_t i = 0; i < numVrts; ++i)
					grid.mark(vrts[i]);
			}
		}
	}

	grid.end_marking();

//	now erase all doubles
	for(size_t i = 0; i < doubles.size(); ++i){
	//	this allows listeners to take actions
		grid.objects_will_be_merged(doubles[i].second, doubles[i].second,
									doubles[i].first);
		grid.erase(doubles[i].first);
	}
}

}//	end of namespace

#endif