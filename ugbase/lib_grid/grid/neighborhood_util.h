/*!
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__LIB_GRID__NEIGHBORHOOD_UTIL__
#define __H__LIB_GRID__NEIGHBORHOOD_UTIL__

#include "grid.h"

namespace ug {

/*!
 * \brief Finds the neighbor connected through a side.
 * \param[in] g
 * \param[in] face
 * \param[in] elem
 * If such a neighbor does not exist, nullptr is returned.
 */
template <typename TBaseElem>
TBaseElem* GetConnectedNeighbor(Grid& g, typename TBaseElem::side* face, TBaseElem* elem);

/*!
 * \brief Finds the neighborhood of a given size for specified element and type
 * \param[in] grid
 * \param[in] extSize size of neighborhood
 * \param[in] elem start element
 * \param[out] begin iterator
 * \param[out] end iterator
 * Returns iterators begin and end for accessing elements of type
 */
template <typename TElem>
void GetNeighborhood
(
	Grid& grid,
	size_t extSize,
	TElem* elem,
	typename geometry_traits<TElem>::const_iterator& begin,
	typename geometry_traits<TElem>::const_iterator& end
);

} // namespace ug

#include "neighborhood_util_impl.hpp"

#endif // __H__LIB_GRID__NEIGHBORHOOD_UTIL__
