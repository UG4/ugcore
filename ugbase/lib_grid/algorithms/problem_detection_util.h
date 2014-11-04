// created by Sebastian Reiter
// s.b.reiter@gmail.com

#ifndef __H__UG_problem_detection_util
#define __H__UG_problem_detection_util

#include <vector>
#include "common/math/ugmath_types.h"

namespace ug{

/**	Checks for all edges of the given tetrahedron whether the distance to the
 *	opposite edge is smaller than the given ratio of the length of the
 * longest edge in the element.
 * \returns	integer to the first edge which is closer to its opposed edge
 *			than the given threshold allows. -1 if no such edge exists.
 *			Use tet_rules::OPPOSED_EDGE to retrieve the associated opposed edge.*/
UG_API int IsSliver(const vector3& v0, const vector3& v1, const vector3& v2,
			  		 const vector3& v3, number thresholdRatio);

///	Searchs for slivers in the given list of elements.
/**	Slivers are flat tetrahedrons. Only tetrahedral elements will be regarded!*/
template <class TIter, class TAAPos>
size_t FindSlivers(std::vector<typename TIter::value_type>& sliversOut,
				 TIter elemsBegin, TIter elemsEnd, number thresholdRatio,
				 TAAPos aaPos, bool clearContainer = true);

}//	end of namespace

////////////////////////////////////////
//	include implementation
#include "problem_detection_util_impl.h"

#endif	//__H__problem_detection_util
