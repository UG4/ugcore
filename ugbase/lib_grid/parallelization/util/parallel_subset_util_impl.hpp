//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m07 d14

#ifndef __H__LIB_GRID__PARALLELL_SUBSET_UTIL_IMPL__
#define __H__LIB_GRID__PARALLELL_SUBSET_UTIL_IMPL__

#include "../distributed_grid.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
//	CreateSurfaceView
template <class TIterator>
void CreateSurfaceView(SubsetHandler& shSurfaceViewOut,
						DistributedGridManager& distGridMgr,
						ISubsetHandler& sh, TIterator iterBegin,
						TIterator iterEnd)
{
	MultiGrid* pMG = dynamic_cast<MultiGrid*>(distGridMgr.get_assigned_grid());
	if(!pMG){
		UG_LOG("  Can't create surface-view. A Multigrid is required.\n");
		return;
	}

	while(iterBegin != iterEnd)
	{
		if(!distGridMgr.is_ghost(*iterBegin)){
			if(!pMG->has_children(*iterBegin)){
				shSurfaceViewOut.assign_subset(*iterBegin,
										sh.get_subset_index(*iterBegin));
			}
		}
		++iterBegin;
	}
}

}//	end of namespace

#endif
