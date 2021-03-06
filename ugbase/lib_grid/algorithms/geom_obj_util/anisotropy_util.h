/*
 * Copyright (c) 2009-2015:  G-CSC, Goethe University Frankfurt
 * Author: Markus Breit
 * Date: 2018-05-25
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

#ifndef LIB_GRID__ALGORITHMS__GEOM_OBJ_UTIL__ANISOTROPY_UTIL_H
#define LIB_GRID__ALGORITHMS__GEOM_OBJ_UTIL__ANISOTROPY_UTIL_H

#include "common/types.h"
#include "lib_grid/grid/grid_base_objects.h"
#include "lib_grid/multi_grid.h"

#include <cstddef>
#include <vector>

namespace ug {


enum AnisotropyState
{
	ISOTROPIC = 0,
	QUAD_SHORTX,
	QUAD_SHORTY,
	PRISM_FLAT,
	PRISM_LONG,
	HEX_SHORTX,
	HEX_SHORTY,
	HEX_SHORTZ,
	HEX_SHORTXY,
	HEX_SHORTXZ,
	HEX_SHORTYZ
};


template <typename TAAPos>
AnisotropyState is_anisotropic
(
	Edge* elem,
	const TAAPos& aaPos,
	number thresholdRatio
);


template <typename TAAPos>
AnisotropyState is_anisotropic
(
	Face* elem,
	const TAAPos& aaPos,
	number thresholdRatio
);


template <typename TAAPos>
AnisotropyState is_anisotropic
(
	Volume* elem,
	const TAAPos& aaPos,
	number thresholdRatio
);



template <typename TAAPos>
AnisotropyState close_sides_of_anisotropic_elem
(
	Edge* elem,
	Grid& grid,
	const TAAPos& aaPos,
	number thresholdRatio,
	std::vector<Vertex*>& sidesOut
);


template <typename TAAPos>
AnisotropyState close_sides_of_anisotropic_elem
(
	Face* elem,
	Grid& grid,
	const TAAPos& aaPos,
	number thresholdRatio,
	std::vector<Edge*>& sidesOut
);


template <typename TAAPos>
AnisotropyState close_sides_of_anisotropic_elem
(
	Volume* elem,
	Grid& grid,
	const TAAPos& aaPos,
	number thresholdRatio,
	std::vector<Face*>& sidesOut
);



template <typename TAAPos>
AnisotropyState long_edges_of_anisotropic_elem
(
	Edge* elem,
	Grid& grid,
	const TAAPos& aaPos,
	number thresholdRatio,
	std::vector<Edge*>& longEdges
);


template <typename TAAPos>
AnisotropyState long_edges_of_anisotropic_elem
(
	Face* elem,
	Grid& grid,
	const TAAPos& aaPos,
	number thresholdRatio,
	std::vector<Edge*>& longEdges
);


template <typename TAAPos>
AnisotropyState long_edges_of_anisotropic_elem
(
	Volume* elem,
	Grid& grid,
	const TAAPos& aaPos,
	number thresholdRatio,
	std::vector<Edge*>& longEdges
);



} // namespace ug

#include "anisotropy_util_impl.h"

#endif // LIB_GRID__ALGORITHMS__GEOM_OBJ_UTIL__ANISOTROPY_UTIL_H
