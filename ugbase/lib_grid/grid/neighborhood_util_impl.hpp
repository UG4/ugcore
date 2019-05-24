/*!
 * Copyright (c) 2010-2019:  G-CSC, Goethe University Frankfurt
 * Author: Stephan Grein
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

#include "neighborhood_util.h"
#include "grid.h"
#include "lib_grid/selector.h"
#include "lib_grid/algorithms/selection_util.h"

namespace ug {
	////////////////////////////////////////////////////////////////////////
	template <typename TElem>
	void GetNeighborhood
	(
		Grid& grid,
		size_t extSize,
		TElem* elem,
		typename geometry_traits<TElem>::const_iterator& begin,
		typename geometry_traits<TElem>::const_iterator& end
	) {
		Selector sel(grid);
		typedef typename geometry_traits<TElem>::const_iterator Iter;
		sel.template select<TElem>(elem);
		grid.begin_marking();
		for (size_t extIters = 0; extIters < extSize; ++extIters)
		{
			SelectAssociatedGridObjects(sel, 1);
			for(size_t lvl = 0; lvl < sel.num_levels(); ++lvl){
				for(Iter iter = sel.template begin<TElem>(lvl);
						iter != sel.template end<TElem>(lvl); ++iter)
				{
					TElem* el = *iter;
					if(!grid.is_marked(el)){
						grid.mark(el);
						typename Grid::traits<TElem>::secure_container elemsOut;
						grid.template associated_elements<TElem>(elemsOut, el);
						size_t size = elemsOut.size();
						for (size_t i = 0; i < size; i++) {
							sel.template select<TElem>(elemsOut[i], 1);
						}
					}
				}
			}
		}
		grid.end_marking();
		begin = sel.template begin<TElem>();
		end = sel.template end<TElem>();
	}
} // namespace ug
