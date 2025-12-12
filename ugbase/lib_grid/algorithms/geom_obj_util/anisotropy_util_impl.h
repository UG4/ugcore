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

#ifndef LIB_GRID_ALGORITHM_ANISOTROPY_UTIL_H
#define LIB_GRID_ALGORITHM_ANISOTROPY_UTIL_H

#include "anisotropy_util.h"

#include "lib_grid/algorithms/element_side_util.h"        // for GetOpposingSide
//#include "lib_grid/grid/grid_util.h"                      // for CompareVertices


namespace ug {


template <typename TAAPos>
AnisotropyState is_anisotropic
(
	Edge* elem,
	const TAAPos& aaPos,
	number thresholdRatio
)
{
	return AnisotropyState::ISOTROPIC;
}



template <typename TAAPos>
static AnisotropyState is_anisotropic
(
	Quadrilateral* q,
	const TAAPos& aaPos,
	number thresholdRatio
)
{
	// check whether elem is anisotropic
	number sideLength02 = VertexDistance(q->vertex(0), q->vertex(1), aaPos)
		+ VertexDistance(q->vertex(2), q->vertex(3), aaPos);

	number sideLength13 = VertexDistance(q->vertex(1), q->vertex(2), aaPos)
			+ VertexDistance(q->vertex(3), q->vertex(0), aaPos);

	if (sideLength02 < thresholdRatio * sideLength13)
		return AnisotropyState::QUAD_SHORTX;
	if (sideLength13 < thresholdRatio * sideLength02)
		return AnisotropyState::QUAD_SHORTY;

	return AnisotropyState::ISOTROPIC;
}



template <typename TAAPos>
AnisotropyState is_anisotropic
(
	Face* elem,
	const TAAPos& aaPos,
	number thresholdRatio
)
{
	auto q = dynamic_cast<Quadrilateral*>(elem);
	if (!q)
		return AnisotropyState::ISOTROPIC;

	return is_anisotropic(q, aaPos, thresholdRatio);
}




template <typename TAAPos>
static AnisotropyState is_anisotropic
(
	Prism* p,
	const TAAPos& aaPos,
	number thresholdRatio
)
{
	// check whether elem is anisotropic
	number length = VertexDistance(p->vertex(0), p->vertex(3), aaPos)
					+ VertexDistance(p->vertex(1), p->vertex(4), aaPos)
					+ VertexDistance(p->vertex(2), p->vertex(5), aaPos);
	number width = VertexDistance(p->vertex(0), p->vertex(1), aaPos)
					+ VertexDistance(p->vertex(1), p->vertex(2), aaPos)
					+ VertexDistance(p->vertex(2), p->vertex(0), aaPos)
					+ VertexDistance(p->vertex(3), p->vertex(4), aaPos)
					+ VertexDistance(p->vertex(4), p->vertex(5), aaPos)
					+ VertexDistance(p->vertex(5), p->vertex(3), aaPos);

	number ratio = width ? 2.0*length/width : std::numeric_limits<number>::max();

	// flat case
	if (ratio < thresholdRatio)
		return AnisotropyState::PRISM_FLAT;

	// long case
	if (ratio*thresholdRatio > 1)
		return AnisotropyState::PRISM_LONG;

	return AnisotropyState::ISOTROPIC;
}



template <typename TAAPos>
static AnisotropyState is_anisotropic
(
	Hexahedron* hex,
	const TAAPos& aaPos,
	number thresholdRatio
)
{
	number length1 = VertexDistance(hex->vertex(0), hex->vertex(1), aaPos)
					+ VertexDistance(hex->vertex(2), hex->vertex(3), aaPos)
					+ VertexDistance(hex->vertex(4), hex->vertex(5), aaPos)
					+ VertexDistance(hex->vertex(6), hex->vertex(7), aaPos);
	number length2 = VertexDistance(hex->vertex(0), hex->vertex(3), aaPos)
					+ VertexDistance(hex->vertex(1), hex->vertex(2), aaPos)
					+ VertexDistance(hex->vertex(4), hex->vertex(7), aaPos)
					+ VertexDistance(hex->vertex(5), hex->vertex(6), aaPos);
	number length3 = VertexDistance(hex->vertex(0), hex->vertex(4), aaPos)
					+ VertexDistance(hex->vertex(1), hex->vertex(5), aaPos)
					+ VertexDistance(hex->vertex(2), hex->vertex(6), aaPos)
					+ VertexDistance(hex->vertex(3), hex->vertex(7), aaPos);

	bool shortx = false;
	bool shorty = false;
	bool shortz = false;

	if (length1 < thresholdRatio * length2 || length1 < thresholdRatio * length3)
		shortx = true;

	if (length2 < thresholdRatio * length1 || length2 < thresholdRatio * length3)
		shorty = true;

	if (length3 < thresholdRatio * length1 || length3 < thresholdRatio * length2)
		shortz = true;

	if (shortx)
	{
		if (shorty)
			return AnisotropyState::HEX_SHORTXY;
		if (shortz)
			return AnisotropyState::HEX_SHORTXZ;
		return AnisotropyState::HEX_SHORTX;
	}

	if (shorty)
	{
		if (shortz)
			return AnisotropyState::HEX_SHORTYZ;
		return AnisotropyState::HEX_SHORTY;
	}

	if (shortz)
		return AnisotropyState::HEX_SHORTZ;

	return AnisotropyState::ISOTROPIC;
}



template <typename TAAPos>
AnisotropyState is_anisotropic
(
	Volume* elem,
	const TAAPos& aaPos,
	number thresholdRatio
)
{
	// treat prism case
	auto prism = dynamic_cast<Prism*>(elem);
	if (prism)
		return is_anisotropic(prism, aaPos, thresholdRatio);

	// treat hexahedron case
	auto hex = dynamic_cast<Hexahedron*>(elem);
	if (hex)
		return is_anisotropic(hex, aaPos, thresholdRatio);

	// other cases do not exist
	return AnisotropyState::ISOTROPIC;
}



template <typename TAAPos>
AnisotropyState close_sides_of_anisotropic_elem
(
	Edge* elem,
	Grid& grid,
	const TAAPos& aaPos,
	number thresholdRatio,
	std::vector<Vertex*>& sidesOut
)
{
	return AnisotropyState::ISOTROPIC;
}


template <typename TAAPos>
AnisotropyState close_sides_of_anisotropic_elem
(
	Face* elem,
	Grid& grid,
	const TAAPos& aaPos,
	number thresholdRatio,
	std::vector<Edge*>& sidesOut
)
{
	// check whether this is a quadrilateral
	auto q = dynamic_cast<Quadrilateral*>(elem);
	if (!q)
		return AnisotropyState::ISOTROPIC;

	// check whether element is anisotropic (and which case)
	AnisotropyState state = is_anisotropic(q, aaPos, thresholdRatio);
	if (state == AnisotropyState::ISOTROPIC)
		return state;

	if (state == AnisotropyState::QUAD_SHORTX)
	{
		Grid::SecureEdgeContainer assEd;
		grid.associated_elements_sorted(assEd, q);
		sidesOut.push_back(assEd[1]);
		sidesOut.push_back(assEd[3]);

		return state;
	}

	if (state == AnisotropyState::QUAD_SHORTY)
	{
		Grid::SecureEdgeContainer assEd;
		grid.associated_elements_sorted(assEd, q);
		sidesOut.push_back(assEd[0]);
		sidesOut.push_back(assEd[2]);

		return state;
	}

	return state;
}



template <typename TAAPos>
AnisotropyState close_sides_of_anisotropic_elem
(
	Volume* elem,
	Grid& grid,
	const TAAPos& aaPos,
	number thresholdRatio,
	std::vector<Face*>& sidesOut
)
{
	// treat prism case
	auto prism = dynamic_cast<Prism*>(elem);
	if (prism)
	{
		AnisotropyState state = is_anisotropic(prism, aaPos, thresholdRatio);
		if (state == AnisotropyState::ISOTROPIC)
			return state;

		// flat case
		if (state == AnisotropyState::PRISM_FLAT)
		{
			Face* side = grid.get_side(prism, 0);
			if (side)
				sidesOut.push_back(side);

			side = grid.get_side(prism, 4);
			if (side)
				sidesOut.push_back(side);
		}

		// long case
		if (state == AnisotropyState::PRISM_LONG)
		{
			Face* side = grid.get_side(prism, 1);
			if (side)
				sidesOut.push_back(side);

			side = grid.get_side(prism, 2);
			if (side)
				sidesOut.push_back(side);

			side = grid.get_side(prism, 3);
			if (side)
				sidesOut.push_back(side);
		}

		return state;
	}


	// treat heaxahedron case
	auto hex = dynamic_cast<Hexahedron*>(elem);
	if (hex)
	{
		AnisotropyState state = is_anisotropic(hex, aaPos, thresholdRatio);
		if (state == AnisotropyState::ISOTROPIC)
			return state;


		// short in x direction
		if (state == AnisotropyState::HEX_SHORTX || state == AnisotropyState::HEX_SHORTXY || state == AnisotropyState::HEX_SHORTXZ)
		{
			Face* side = grid.get_side(hex, 2);
			if (side)
				sidesOut.push_back(side);

			side = grid.get_side(hex, 4);
			if (side)
				sidesOut.push_back(side);
		}

		// short in y direction
		if (state == AnisotropyState::HEX_SHORTY || state == AnisotropyState::HEX_SHORTXY || state == AnisotropyState::HEX_SHORTYZ)
		{
			Face* side = grid.get_side(hex, 1);
			if (side)
				sidesOut.push_back(side);

			side = grid.get_side(hex, 3);
			if (side)
				sidesOut.push_back(side);
		}

		// short in z direction
		if (state == AnisotropyState::HEX_SHORTZ || state == AnisotropyState::HEX_SHORTXZ || state == AnisotropyState::HEX_SHORTYZ)
		{
			Face* side = grid.get_side(hex, 0);
			if (side)
				sidesOut.push_back(side);

			side = grid.get_side(hex, 5);
			if (side)
				sidesOut.push_back(side);
		}

		return state;
	}


	// other elements cannot be anisotropic
	return AnisotropyState::ISOTROPIC;
}





template <typename TAAPos>
AnisotropyState long_edges_of_anisotropic_elem
(
	Edge* elem,
	Grid& grid,
	const TAAPos& aaPos,
	number thresholdRatio,
	std::vector<Edge*>& longEdges
)
{
	return AnisotropyState::ISOTROPIC;
}


template <typename TAAPos>
AnisotropyState long_edges_of_anisotropic_elem
(
	Face* elem,
	Grid& grid,
	const TAAPos& aaPos,
	number thresholdRatio,
	std::vector<Edge*>& longEdges
)
{
	return close_sides_of_anisotropic_elem(elem, grid, aaPos, thresholdRatio, longEdges);
}


template <typename TAAPos>
AnisotropyState long_edges_of_anisotropic_elem
(
	Volume* elem,
	Grid& grid,
	const TAAPos& aaPos,
	number thresholdRatio,
	std::vector<Edge*>& longEdges
)
{
	// treat prism case
	auto prism = dynamic_cast<Prism*>(elem);
	if (prism)
	{
		AnisotropyState state = is_anisotropic(prism, aaPos, thresholdRatio);
		if (state == AnisotropyState::ISOTROPIC)
			return state;

		// flat case
		if (state == AnisotropyState::PRISM_FLAT)
		{
			Grid::SecureEdgeContainer assEd;
			grid.associated_elements_sorted(assEd, prism);

			UG_COND_THROW(assEd.size() < 9, "Prism needs to have 9 edges, but only "
				<< assEd.size() << " found.");

			longEdges.reserve(longEdges.size() + 6);
			longEdges.push_back(assEd[0]);
			longEdges.push_back(assEd[1]);
			longEdges.push_back(assEd[2]);
			longEdges.push_back(assEd[6]);
			longEdges.push_back(assEd[7]);
			longEdges.push_back(assEd[8]);

			return state;
		}

		// long case
		if (state == AnisotropyState::PRISM_LONG)
		{
			Grid::SecureEdgeContainer assEd;
			grid.associated_elements_sorted(assEd, prism);

			UG_COND_THROW(assEd.size() < 9, "Prism needs to have 9 edges, but only "
				<< assEd.size() << " found.");

			longEdges.reserve(longEdges.size() + 3);
			longEdges.push_back(assEd[3]);
			longEdges.push_back(assEd[4]);
			longEdges.push_back(assEd[5]);

			return state;
		}

		return state;
	}


	// treat heaxahedron case
	auto hex = dynamic_cast<Hexahedron*>(elem);
	if (hex)
	{
		AnisotropyState state = is_anisotropic(hex, aaPos, thresholdRatio);
		if (state == AnisotropyState::ISOTROPIC)
			return state;


		// short in x direction
		if (state == AnisotropyState::HEX_SHORTX)
		{
			Grid::SecureEdgeContainer assEd;
			grid.associated_elements_sorted(assEd, hex);

			UG_COND_THROW(assEd.size() < 12, "Hexahedron needs to have 12 edges, but only "
				<< assEd.size() << " found.");

			longEdges.reserve(longEdges.size() + 8);
			longEdges.push_back(assEd[1]);
			longEdges.push_back(assEd[3]);
			longEdges.push_back(assEd[4]);
			longEdges.push_back(assEd[5]);
			longEdges.push_back(assEd[6]);
			longEdges.push_back(assEd[7]);
			longEdges.push_back(assEd[8]);
			longEdges.push_back(assEd[11]);

			return state;
		}

		// short in y direction
		if (state == AnisotropyState::HEX_SHORTY)
		{
			Grid::SecureEdgeContainer assEd;
			grid.associated_elements_sorted(assEd, hex);

			UG_COND_THROW(assEd.size() < 12, "Hexahedron needs to have 12 edges, but only "
				<< assEd.size() << " found.");

			longEdges.reserve(longEdges.size() + 8);
			longEdges.push_back(assEd[0]);
			longEdges.push_back(assEd[2]);
			longEdges.push_back(assEd[4]);
			longEdges.push_back(assEd[5]);
			longEdges.push_back(assEd[6]);
			longEdges.push_back(assEd[7]);
			longEdges.push_back(assEd[8]);
			longEdges.push_back(assEd[10]);

			return state;
		}

		// short in z direction
		if (state == AnisotropyState::HEX_SHORTZ)
		{
			Grid::SecureEdgeContainer assEd;
			grid.associated_elements_sorted(assEd, hex);

			UG_COND_THROW(assEd.size() < 12, "Hexahedron needs to have 12 edges, but only "
				<< assEd.size() << " found.");

			longEdges.reserve(longEdges.size() + 8);
			longEdges.push_back(assEd[0]);
			longEdges.push_back(assEd[1]);
			longEdges.push_back(assEd[2]);
			longEdges.push_back(assEd[3]);
			longEdges.push_back(assEd[8]);
			longEdges.push_back(assEd[9]);
			longEdges.push_back(assEd[10]);
			longEdges.push_back(assEd[11]);

			return state;
		}

		// short in x and y direction
		if (state == AnisotropyState::HEX_SHORTXY)
		{
			Grid::SecureEdgeContainer assEd;
			grid.associated_elements_sorted(assEd, hex);

			UG_COND_THROW(assEd.size() < 12, "Hexahedron needs to have 12 edges, but only "
				<< assEd.size() << " found.");

			longEdges.reserve(longEdges.size() + 4);
			longEdges.push_back(assEd[4]);
			longEdges.push_back(assEd[5]);
			longEdges.push_back(assEd[6]);
			longEdges.push_back(assEd[7]);

			return state;
		}

		// short in x and z direction
		if (state == AnisotropyState::HEX_SHORTXZ)
		{
			Grid::SecureEdgeContainer assEd;
			grid.associated_elements_sorted(assEd, hex);

			UG_COND_THROW(assEd.size() < 12, "Hexahedron needs to have 12 edges, but only "
				<< assEd.size() << " found.");

			longEdges.reserve(longEdges.size() + 4);
			longEdges.push_back(assEd[1]);
			longEdges.push_back(assEd[3]);
			longEdges.push_back(assEd[9]);
			longEdges.push_back(assEd[11]);

			return state;
		}

		// short in y and z direction
		if (state == AnisotropyState::HEX_SHORTYZ)
		{
			Grid::SecureEdgeContainer assEd;
			grid.associated_elements_sorted(assEd, hex);

			UG_COND_THROW(assEd.size() < 12, "Hexahedron needs to have 12 edges, but only "
				<< assEd.size() << " found.");

			longEdges.reserve(longEdges.size() + 4);
			longEdges.push_back(assEd[0]);
			longEdges.push_back(assEd[2]);
			longEdges.push_back(assEd[8]);
			longEdges.push_back(assEd[10]);

			return state;
		}

		return state;
	}


	// other elements cannot be anisotropic
	return AnisotropyState::ISOTROPIC;
}

} // namespace ug

#endif