// created by Sebastian Reiter
// s.b.reiter@gmail.com

#include "problem_detection_util.h"
#include "lib_grid/grid_objects/tetrahedron_rules.h"
#include "common/math/misc/math_util.h"

namespace ug{

int IsSliver(const vector3& v0, const vector3& v1, const vector3& v2,
			  const vector3& v3, number thresholdRatio)
{
	using namespace tet_rules;
	vector3 v[] = {v0, v1, v2, v3};

	number maxLenSq = 0;
	for(int iedge = 0; iedge < NUM_EDGES; ++iedge){
		maxLenSq = VecDistanceSq(v[EDGE_VRT_INDS[iedge][0]], v[EDGE_VRT_INDS[iedge][1]]);
	}

	number thresholdDist = sqrt(maxLenSq) * thresholdRatio;

	for(int iedge = 0; iedge + 1 < NUM_EDGES; ++iedge){
		int iop = OPPOSED_EDGE[iedge];
		if(iop > iedge){
			number dist = DistanceLineToLine(v[EDGE_VRT_INDS[iedge][0]], v[EDGE_VRT_INDS[iedge][1]],
											 v[EDGE_VRT_INDS[iop][0]], v[EDGE_VRT_INDS[iop][1]]);
			if(dist < thresholdDist)
				return iedge;
		}
	}

	return -1;
}

}//	end of namespace
