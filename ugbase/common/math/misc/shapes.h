// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 27.07.2012 (d.m.y())

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

	vector_t min;
	vector_t max;
};

}//	end of namespace


#include "shapes_impl.h"

#endif
