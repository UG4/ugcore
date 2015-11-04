#ifndef __H__UG_problem_detection_util_impl
#define __H__UG_problem_detection_util_impl

#include "problem_detection_util.h"
#include "lib_grid/grid/grid_base_objects.h"

namespace ug{

template <class TIter, class TAAPos>
size_t FindSlivers(std::vector<typename TIter::value_type>& sliversOut,
				 TIter elemsBegin, TIter elemsEnd, number thresholdRatio,
				 TAAPos aaPos, bool clearContainer)
{
	typedef typename TIter::value_type elem_ptr_t;

	if(clearContainer)
		sliversOut.clear();

	size_t numInitialSlivers = sliversOut.size();
	for(TIter iter = elemsBegin; iter != elemsEnd; ++iter){
		elem_ptr_t e = *iter;
		Volume::ConstVertexArray v = e->vertices();
		if((e->reference_object_id() == ROID_TETRAHEDRON) &&
			(IsSliver(aaPos[v[0]], aaPos[v[1]], aaPos[v[2]], aaPos[v[3]], thresholdRatio) != -1))
		{
			sliversOut.push_back(e);
		}
	}

	return sliversOut.size() - numInitialSlivers;
}

}//	end of namespace

#endif	//__H__problem_detection_util_impl
