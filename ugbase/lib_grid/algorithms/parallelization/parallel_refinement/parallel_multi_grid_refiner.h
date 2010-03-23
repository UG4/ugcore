//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m03 d23

#ifndef __H__LIB_GRID__PARALLEL_MULTI_GRID_REFINER__
#define __H__LIB_GRID__PARALLEL_MULTI_GRID_REFINER__

#include "lib_grid/lg_base.h"
#include "lib_grid/multi_grid.h"
#include "lib_grid/algorithms/refinement/multi_grid_refiner.h"
#include "../grid_distribution.h"

namespace ug
{

class ParallelMultiGridRefiner : public MultiGridRefiner
{
	public:
		ParallelMultiGridRefiner();
		ParallelMultiGridRefiner(MultiGrid& mg, GridLayoutMap& layoutMap);
		~ParallelMultiGridRefiner();

	protected:
		virtual void collect_objects_for_refine();
		
	protected:
		GridLayoutMap*	m_pLayoutMap;
};

}//	end of namespace

#endif
