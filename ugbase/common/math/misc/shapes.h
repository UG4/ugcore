/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG_shapes__
#define __H__UG_shapes__

#include "math_util.h"

namespace ug{

template <class vector_t>
class Sphere{
	public:
		Sphere()	{}
		Sphere(const vector_t& nCenter, number nRadius) :
			center(nCenter), radius(nRadius)	{}

	///	returns true if the point lies inside or on the bounds of the rectangle
		bool contains_point(const vector_t& point) const;

	///	returns true if the specified sphere touches or intersects this sphere.
		bool intersects(const Sphere& sphere) const;

		vector_t	center;
		number		radius;
};


template <class vector_t>
struct AABox{
	AABox()	{}
	AABox(const vector_t& nMin, const vector_t& nMax) : min(nMin), max(nMax)	{}

///	calculates the bounding box to the given set of points.
	AABox(const vector_t* points, size_t numPoints);

///	calculates the bounding box of two given bounding boxes
	AABox(const AABox& b1, const AABox& b2);

///	calculates the bounding box of the given sphere
	AABox(const Sphere<vector_t>& s);

///	calculates the bounding box of a given box and an additional point
	AABox(const AABox& b, const vector_t& v);

///	returns the center of the box
	vector_t center() const;

///	returns the extension (width/height/depth) of the box
	vector_t extension() const;

///	returns true if the given point lies in the box or on its boundary
	bool contains_point(const vector_t& point) const;

/// return true if the given line (segment) and the box overlap
	bool overlaps_line(const vector_t& point1, const vector_t& point2) const;

	vector_t min;
	vector_t max;
};

}//	end of namespace


#include "shapes_impl.h"

#endif
