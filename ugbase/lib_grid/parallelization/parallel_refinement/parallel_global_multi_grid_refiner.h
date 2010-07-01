//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m04 d14

#ifndef __H__LIB_GRID__PARALLEL_GLOBAL_MULTI_GRID_REFINER__
#define __H__LIB_GRID__PARALLEL_GLOBAL_MULTI_GRID_REFINER__

#include <vector>
#include "lib_grid/lg_base.h"
#include "lib_grid/multi_grid.h"
#include "lib_grid/algorithms/refinement/global_multi_grid_refiner.h"
#include "../distributed_grid.h"

namespace ug
{

class ParallelGlobalMultiGridRefiner : public GlobalMultiGridRefiner
{
	public:
		//ParallelMultiGridRefiner();
		ParallelGlobalMultiGridRefiner(DistributedGridManager& distGridMgr);
		virtual ~ParallelGlobalMultiGridRefiner();

	protected:
		virtual bool refinement_is_allowed(VertexBase* elem);
		virtual bool refinement_is_allowed(EdgeBase* elem);
		virtual bool refinement_is_allowed(Face* elem);
		virtual bool refinement_is_allowed(Volume* elem);
		
		virtual void refinement_step_begins();
		virtual void refinement_step_ends();

	protected:
		DistributedGridManager& m_distGridMgr;
};

}//	end of namespace

#endif
