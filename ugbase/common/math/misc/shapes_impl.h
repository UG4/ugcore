#ifndef __H__UG__shapes_impl__
#define __H__UG__shapes_impl__

#include "shapes.h"

namespace ug{

template <class vector_t>
bool Sphere<vector_t>::contains_point(const vector_t& p) const
{
	return VecDistanceSq(center, p) <= sq(radius);
}

template <class vector_t>
bool Sphere<vector_t>::intersects(const Sphere& sphere) const
{
	return VecDistanceSq(center, sphere.center)
		   <= sq(radius + sphere.radius);
}



template <class vector_t>
AABox<vector_t>::AABox(const vector_t* points, size_t numPoints)
{
	if(numPoints < 1){
		VecSet(min, 0);
		VecSet(max, 0);
		return;
	}

	min = max = points[0];

	for(size_t i = 1; i < numPoints; ++i){
		const vector_t& p = points[i];
		for(size_t j = 0; j < vector_t::Size; ++j){
			min[j] = std::min(min[j], p[j]);
			max[j] = std::max(max[j], p[j]);
		}
	}
}

template <class vector_t>
AABox<vector_t>::AABox(const AABox& b1, const AABox& b2)
{
	for(size_t j = 0; j < vector_t::Size; ++j){
		min[j] = std::min(b1.min[j], b2.min[j]);
		max[j] = std::max(b1.max[j], b2.max[j]);
	}
}

template <class vector_t>
AABox<vector_t>::AABox(const Sphere<vector_t>& s)
{
	vector_t vrad;
	SetVec(vrad, s.radius);
	VecSubtract(min, s.center, vrad);
	VecAdd(max, s.center, vrad);
}

template <class vector_t>
AABox<vector_t>::AABox(const AABox& b, const vector_t& v)
{
	for(size_t j = 0; j < vector_t::Size; ++j){
		min[j] = std::min(b.min[j], v[j]);
		max[j] = std::max(b.max[j], v[j]);
	}
}

template <class vector_t>
vector_t AABox<vector_t>::center() const
{
	vector_t v;
	VecAdd(v, max, min);
	VecScale(v, v, 0.5);
	return v;
}

template <class vector_t>
vector_t AABox<vector_t>::extension() const
{
	vector_t v;
	VecSubtract(v, max, min);
	return v;
}

template <class vector_t>
bool AABox<vector_t>::contains_point(const vector_t& point) const
{
	return BoxBoundProbe(point, min, max);
}

template <class vector_t>
bool AABox<vector_t>::overlaps_line(const vector_t& point1, const vector_t& point2) const
{
	vector_t lmin, lmax;
	lmin = lmax = point1;
	for (size_t i = 0; i < vector_t::Size; i++) {
		lmin[i] = std::min(lmin[i], point2[i]);
		lmax[i] = std::max(lmax[i], point2[i]);
	}
	return BoxBoxIntersection(lmin, lmax, min, max);
}

}// end of namespace

#endif
