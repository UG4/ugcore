
#include "common/common.h"
#include "common/cuthill_mckee.h"
#include "cuthill_mckee.h"
#include <algorithm>
#include <vector>
#include <queue>
#include "common/profiler/profiler.h"
#include "lib_disc/domain.h"

namespace ug{

void OrderCuthillMcKee(DoFDistribution& dofDistr, bool bReverse)
{
	PROFILE_FUNC();
//	get adjacency graph
	std::vector<std::vector<size_t> > vvConnection;
	try{
		dofDistr.get_connections(vvConnection);
	}
	UG_CATCH_THROW("OrderCuthillMcKee: No adjacency graph available.");

//	get mapping for cuthill-mckee order
	std::vector<size_t> vNewIndex;
	ComputeCuthillMcKeeOrder(vNewIndex, vvConnection, bReverse);

//	reorder indices
	dofDistr.permute_indices(vNewIndex);
}

template <typename TDomain>
void OrderCuthillMcKee(ApproximationSpace<TDomain>& approxSpace, bool bReverse)
{
	std::vector<SmartPtr<DoFDistribution> > vDD = approxSpace.dof_distributions();

	for(size_t i = 0; i < vDD.size(); ++i)
		OrderCuthillMcKee(*vDD[i], bReverse);
}

#ifdef UG_DIM_1
template void OrderCuthillMcKee<Domain1d>(ApproximationSpace<Domain1d>& approxSpace, bool bReverse);
#endif
#ifdef UG_DIM_2
template void OrderCuthillMcKee<Domain2d>(ApproximationSpace<Domain2d>& approxSpace, bool bReverse);
#endif
#ifdef UG_DIM_3
template void OrderCuthillMcKee<Domain3d>(ApproximationSpace<Domain3d>& approxSpace, bool bReverse);
#endif

}
