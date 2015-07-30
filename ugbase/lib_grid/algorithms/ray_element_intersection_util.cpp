// created by Sebastian Reiter
// s.b.reiter@gmail.com

#include <algorithm>
#include "ray_element_intersection_util.h"
#include "common/math/misc/math_util.h"
#include "lib_grid/iterators/lg_for_each.h"

using namespace std;

namespace ug{

///	utility method for the full-dimensional RayElementIntersection implementation
template <class TElem, class vector_t>
static bool
RayElementIntersectionImpl(
		number& sminOut,
		number& smaxOut,
		const vector_t& from,
		const vector_t& dir,
		TElem* e,
		Grid& g,
		Grid::VertexAttachmentAccessor<Attachment<vector_t> > aaPos)
{
	typedef typename TElem::side	side_t;

	int numIntersections = 0;
	number smin = 0, smax = 0;

	typename Grid::traits<side_t>::secure_container sides;
	g.associated_elements(sides, e);

	for_each_in_vec(side_t* s, sides){
		number tsmin, tsmax;
		if(RayElementIntersection(tsmin, tsmax, from, dir, s, g, aaPos)){
			if(numIntersections == 0)
				smin = smax = tsmin;
			else{
				smin = min(smin, tsmin);
				smax = max(smax, tsmax);
			}
			++numIntersections;
		}
	}end_for;

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
		Grid::VertexAttachmentAccessor<AVector2> aaPos)
{
	vector2 v;
	number t;
	bool ret = RayLineIntersection2d(v, t, sminOut, aaPos[e->vertex(0)],
									 aaPos[e->vertex(1)], from, dir);
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
		Grid::VertexAttachmentAccessor<AVector2> aaPos)
{
	return RayElementIntersectionImpl(sminOut, smaxOut, from, dir, f, g, aaPos);
}


//	3d face intersection
bool RayElementIntersection(
		number& sminOut,
		number& smaxOut,
		const vector3& from,
		const vector3& dir,
		Face* f,
		Grid&,
		Grid::VertexAttachmentAccessor<AVector3> aaPos)
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
					from, dir);
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
		Grid::VertexAttachmentAccessor<AVector3> aaPos)
{
	return RayElementIntersectionImpl(sminOut, smaxOut, from, dir, v, g, aaPos);
}

//  3d edge intersection
bool RayElementIntersection(
		number& sminOut,
		number& smaxOut,
		const vector3& from,
		const vector3& dir,
		Edge* e,
		Grid& g,
		Grid::VertexAttachmentAccessor<AVector3> aaPos)
{
	return false;
}

}//	end of namespace
