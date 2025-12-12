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

#include "ray_element_intersection_util.h"

#include <algorithm>

#include "common/math/misc/math_util.h"
//ø #include "lib_grid/iterators/lg_for_each.h"

using namespace std;

namespace ug {

///	utility method for the full-dimensional RayElementIntersection implementation
template <typename TElem, typename vector_t>
static bool
RayElementIntersectionImpl(
		number& sminOut,
		number& smaxOut,
		const vector_t& from,
		const vector_t& dir,
		TElem* e,
		Grid& g,
		Grid::VertexAttachmentAccessor<Attachment<vector_t> > aaPos,
		number sml)
{
	using side_t = typename TElem::side;

	int numIntersections = 0;
	number smin = 0, smax = 0;

	typename Grid::traits<side_t>::secure_container sides;
	g.associated_elements(sides, e);

	for(size_t _vfeI = 0; _vfeI < sides.size(); ++_vfeI){ side_t* s = sides[_vfeI];{
		number tsmin, tsmax;
		if(RayElementIntersection(tsmin, tsmax, from, dir, s, g, aaPos, sml)){
			if(numIntersections == 0)
				smin = smax = tsmin;
			else{
				smin = min(smin, tsmin);
				smax = max(smax, tsmax);
			}
			++numIntersections;
		}
	}};

	sminOut = smin;
	smaxOut = smax;
	return numIntersections > 0;
}


//	2d edge intersection
bool RayElementIntersection(
		number& sminOut,
		number& smaxOut,
		const vector2& from,
		const vector2& dir,
		Edge* e,
		Grid&,
		Grid::VertexAttachmentAccessor<AVector2> aaPos,
		number sml)
{
	vector2 v;
	number t;
	bool ret = RayLineIntersection2d(v, t, sminOut, aaPos[e->vertex(0)],
									 aaPos[e->vertex(1)], from, dir, sml);
	smaxOut = sminOut;
	return ret;
}

//	2d face intersection
bool RayElementIntersection(
		number& sminOut,
		number& smaxOut,
		const vector2& from,
		const vector2& dir,
		Face* f,
		Grid& g,
		Grid::VertexAttachmentAccessor<AVector2> aaPos,
		number sml)
{
	return RayElementIntersectionImpl(sminOut, smaxOut, from, dir, f, g, aaPos, sml);
}


//	3d face intersection
bool RayElementIntersection(
		number& sminOut,
		number& smaxOut,
		const vector3& from,
		const vector3& dir,
		Face* f,
		Grid&,
		Grid::VertexAttachmentAccessor<AVector3> aaPos,
		number sml)
{
		vector3 v;
		number t0, t1;
		bool ret = false;
		size_t i = 0;

		do{
			ret = RayTriangleIntersection(
					v, t0, t1, sminOut,
					aaPos[f->vertex(0)],
					aaPos[f->vertex(i + 1)],
					aaPos[f->vertex(i + 2)],
					from, dir, sml);
			++i;
		}while(!ret && (i + 2 < f->num_vertices()));

		smaxOut = sminOut;
		return ret;
}

//	3d volume intersection
bool RayElementIntersection(
		number& sminOut,
		number& smaxOut,
		const vector3& from,
		const vector3& dir,
		Volume* v,
		Grid& g,
		Grid::VertexAttachmentAccessor<AVector3> aaPos,
		number sml)
{
	return RayElementIntersectionImpl(sminOut, smaxOut, from, dir, v, g, aaPos, sml);
}

//  3d edge intersection
bool RayElementIntersection(
		number& sminOut,
		number& smaxOut,
		const vector3& from,
		const vector3& dir,
		Edge* e,
		Grid& g,
		Grid::VertexAttachmentAccessor<AVector3> aaPos,
		number sml)
{
	return false;
}

}//	end of namespace
