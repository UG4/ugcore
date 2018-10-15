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

#include "lib_grid/algorithms/element_side_util.h"        // for GetOpposingSide
#include "lib_grid/grid/grid_util.h"                      // for CompareVertices


namespace ug {


template <typename TAAPos>
bool is_anisotropic
(
	Edge* elem,
	Grid& grid,
	const TAAPos& aaPos,
	number thresholdRatio,
	std::vector<Vertex*>* nb
)
{
	return false;
}



template <typename TAAPos>
bool is_anisotropic
(
	Face* elem,
	Grid& grid,
	const TAAPos& aaPos,
	number thresholdRatio,
	std::vector<Edge*>* nb
)
{
	Quadrilateral* q = dynamic_cast<Quadrilateral*>(elem);
	if (!q)
		return false;

	// check whether elem is anisotropic
	number sideLength02 = VertexDistance(q->vertex(0), q->vertex(1), aaPos)
		+ VertexDistance(q->vertex(2), q->vertex(3), aaPos);

	number sideLength13 = VertexDistance(q->vertex(1), q->vertex(2), aaPos)
			+ VertexDistance(q->vertex(3), q->vertex(0), aaPos);

	number ratio = 0.0;
	if (sideLength02 > sideLength13)
		ratio = sideLength13 / sideLength02;
	else if (sideLength13 > sideLength02)
		ratio = sideLength02 / sideLength13;
	else
		ratio = 1.0;

	if (ratio >= thresholdRatio)
		return false;


	// get elem neighbors at larger sides
	typedef typename Grid::traits<Face>::secure_container elem_list;
	typedef typename Grid::traits<Edge>::secure_container side_list;

	side_list sl;
	EdgeDescriptor largeSide;
	if (sideLength02 > sideLength13)
		largeSide = elem->edge_desc(0);
	else
		largeSide = elem->edge_desc(1);

	grid.associated_elements(sl, elem);
	size_t slSz = sl.size();
	for (size_t s = 0; s < slSz; ++s)
	{
		if (CompareVertices(sl[s], &largeSide))
		{
			Edge* opp = GetOpposingSide(grid, elem, sl[s]);
			if (nb)
			{
				nb->push_back(sl[s]);
				nb->push_back(opp);
			}

			return true;
		}
	}

	return true;
}



template <typename TAAPos>
bool is_anisotropic
(
	Volume* elem,
	Grid& grid,
	const TAAPos& aaPos,
	number thresholdRatio,
	std::vector<Face*>* nb
)
{
	Prism* prism = dynamic_cast<Prism*>(elem);

	// treat prism case
	if (prism)
	{
		// check whether elem is anisotropic
		number length = VertexDistance(prism->vertex(0), prism->vertex(3), aaPos)
						+ VertexDistance(prism->vertex(1), prism->vertex(4), aaPos)
						+ VertexDistance(prism->vertex(2), prism->vertex(5), aaPos);
		number width = VertexDistance(prism->vertex(0), prism->vertex(1), aaPos)
						+ VertexDistance(prism->vertex(1), prism->vertex(2), aaPos)
						+ VertexDistance(prism->vertex(2), prism->vertex(0), aaPos)
						+ VertexDistance(prism->vertex(3), prism->vertex(4), aaPos)
						+ VertexDistance(prism->vertex(4), prism->vertex(5), aaPos)
						+ VertexDistance(prism->vertex(5), prism->vertex(3), aaPos);

		number ratio = width ? 2.0*length/width : std::numeric_limits<number>::max();

		if (ratio >= thresholdRatio)
			return false;


		// get elem neighbors at larger sides
		typedef typename Grid::traits<Volume>::secure_container elem_list;
		typedef typename Grid::traits<Face>::secure_container side_list;

		side_list sl;
		FaceDescriptor largeSide = prism->face_desc(0);

		grid.associated_elements(sl, prism);
		size_t slSz = sl.size();
		for (size_t s = 0; s < slSz; ++s)
		{
			if (CompareVertices(sl[s], &largeSide))
			{
				Face* opp = GetOpposingSide(grid, elem, sl[s]);
				if (nb)
				{
					nb->push_back(sl[s]);
					nb->push_back(opp);
				}

				return true;
			}
		}

		return true;
	}

	Hexahedron* hex = dynamic_cast<Hexahedron*>(elem);
	if (!hex)
		return false;

	// treat hexahedron case

	// check whether elem is anisotropic
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

	number ratio;
	number largestLength = std::max(length1, std::max(length2, length3));
	number smallestLength = std::min(length1, std::min(length2, length3));
	if (largestLength > smallestLength)
		ratio = smallestLength / largestLength;
	else
		ratio = 1.0;

	if (ratio >= thresholdRatio)
		return false;

	// get elem neighbors at larger sides
	typedef typename Grid::traits<Volume>::secure_container elem_list;
	typedef typename Grid::traits<Face>::secure_container side_list;

	size_t smallestDim = length1 < length2 ? (length1 < length3 ? 1 : 3) : (length2 < length3 ? 2 : 3);
	FaceDescriptor largeSide = hex->face_desc(3-smallestDim);

	side_list sl;
	grid.associated_elements(sl, hex);
	size_t slSz = sl.size();
	for (size_t s = 0; s < slSz; ++s)
	{
		if (CompareVertices(sl[s], &largeSide))
		{
			Face* opp = GetOpposingSide(grid, elem, sl[s]);
			if (nb)
			{
				nb->push_back(sl[s]);
				nb->push_back(opp);
			}

			return true;
		}
	}

	return true;
}

} // namespace ug

