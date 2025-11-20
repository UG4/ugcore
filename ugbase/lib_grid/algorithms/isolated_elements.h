/*
 * Copyright (c) 2016:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG_isolated_elements
#define __H__UG_isolated_elements

#include <vector>
#include "lib_grid/grid/grid.h"

namespace ug{

/**	Writes all elements between 'begin' and 'end' which are not sides of
 * elements in 'grid' to 'elemsOut'. This method only makes sense if called
 * on a sequence of vertices, edges, or faces.*/
template <typename TSideIterator>
size_t CollectUnconnectedSides (
			std::vector<typename TSideIterator::value_type>& elemsOut,
			Grid& grid,
			TSideIterator begin,
			TSideIterator end)
{
	using side_t = typename PtrToValueType<typename TSideIterator::value_type>::base_type;
	using sideof_t = typename side_t::sideof;

	elemsOut.clear();

	typename Grid::traits<sideof_t>::secure_container	con;

	for(TSideIterator i = begin; i != end; ++i) {
		grid.associated_elements(con, *i);
		if(con.size() == 0)
			elemsOut.push_back(*i);
	}

	return elemsOut.size();
}

}//	end of namespace

#endif