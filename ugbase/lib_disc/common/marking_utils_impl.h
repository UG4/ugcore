/*
 * Copyright (c) 2011-2018:  G-CSC, Goethe University Frankfurt
 * Author: Markus Breit
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

#ifndef IG_UGBASE_LIB_DISC_COMMON_MARKING_UTILS_IMPL_H
#define IG_UGBASE_LIB_DISC_COMMON_MARKING_UTILS_IMPL_H

#include "marking_utils.h"

#include "lib_disc/common/multi_index.h"  // for DoFIndex
#include "lib_disc/domain_traits.h"  // for domain_traits
#include "lib_grid/refinement/refiner_interface.h"  // for IRefiner


namespace ug {


// //////////////////////////////////////////////////////////////////////////////
// 	Check values are within bounds
// //////////////////////////////////////////////////////////////////////////////
template <typename TGridFunction, typename TBaseElem>
unsigned long MarkOutOfRangeElems
(
	SmartPtr<IRefiner> refiner,
	ConstSmartPtr<TGridFunction> u,
	size_t cmp,
	number lowerBnd,
	number upperBnd
)
{
	using elem_it = typename TGridFunction::template traits<TBaseElem>::const_iterator;
	using elem_type = typename domain_traits<TGridFunction::domain_type::dim>::element_type;
	using elem_list = typename Grid::traits<elem_type>::secure_container;

	Grid& grid = *refiner->grid();

	unsigned long nMarked = 0;
	elem_list el;

	// loop the elements in the subset
	std::vector<DoFIndex> vDI;
	elem_it it = u->template begin<TBaseElem>();
	elem_it itEnd = u->template end<TBaseElem>();
	for (; it != itEnd; ++it)
	{
		TBaseElem* elem = *it;

		// loop indices at this element
		const size_t nInd = u->inner_dof_indices(elem, cmp, vDI, true);
		for (size_t i = 0; i < nInd; ++i)
		{
			const number& val = DoFRef(*u, vDI[i]);
			if (val < lowerBnd || val > upperBnd)
			{
				// mark neighbors for refinement
				grid.associated_elements(el, elem);
				const size_t elSz = el.size();
				for (size_t e = 0; e < elSz; ++e)
				{
					if (refiner->get_mark(el[e]) == RM_NONE)
					{
						refiner->mark(el[e], RM_FULL);
						++nMarked;
					}
				}
			}
		}
	}

	return nMarked;
}

template <typename TGridFunction>
void MarkOutOfRangeElems
(
	SmartPtr<IRefiner> refiner,
	ConstSmartPtr<TGridFunction> u,
	size_t cmp,
	number lowerBnd,
	number upperBnd
)
{
	unsigned long nMarked = 0;
	if (u->max_fct_dofs(cmp, 0))
		nMarked += MarkOutOfRangeElems<TGridFunction, Vertex>(refiner, u, cmp, lowerBnd, upperBnd);
	if (u->max_fct_dofs(cmp, 1))
		nMarked += MarkOutOfRangeElems<TGridFunction, Edge>(refiner, u, cmp, lowerBnd, upperBnd);
	if (u->max_fct_dofs(cmp, 2))
		nMarked += MarkOutOfRangeElems<TGridFunction, Face>(refiner, u, cmp, lowerBnd, upperBnd);
	if (u->max_fct_dofs(cmp, 3))
		nMarked += MarkOutOfRangeElems<TGridFunction, Volume>(refiner, u, cmp, lowerBnd, upperBnd);

#ifdef UG_PARALLEL
	if (pcl::NumProcs() > 1)
	{
		pcl::ProcessCommunicator pc;
		nMarked = pc.allreduce(nMarked, PCL_RO_SUM);
	}
#endif


	if (nMarked)
		UG_LOGN("  +++ Marked for refinement: " << nMarked << " elements.");
}

} // namespace ug

#endif