//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m03 d23

#ifndef __H__LIB_GRID__PARALLEL_MULTI_GRID_REFINER__
#define __H__LIB_GRID__PARALLEL_MULTI_GRID_REFINER__

#include <vector>
#include "lib_grid/lg_base.h"
#include "lib_grid/multi_grid.h"
#include "lib_grid/algorithms/refinement/multi_grid_refiner.h"
#include "../distributed_grid.h"

namespace ug
{

class ParallelMultiGridRefiner : public MultiGridRefiner
{
	public:
		//ParallelMultiGridRefiner();
		ParallelMultiGridRefiner(DistributedGridManager& distGridMgr);
		~ParallelMultiGridRefiner();

	protected:
		virtual void collect_objects_for_refine();

		virtual void refinement_step_begins();
		virtual void refinement_step_ends();
		
		virtual void set_rule(VertexBase* e, RefinementMark mark);
		virtual void set_rule(EdgeBase* e, RefinementMark mark);
		virtual void set_rule(Face* e, RefinementMark mark);
		virtual void set_rule(Volume* e, RefinementMark mark);

	protected:
		DistributedGridManager& m_distGridMgr;
		
		std::vector<VertexBase*>	m_vNewlyMarkedVertices;
		std::vector<EdgeBase*>		m_vNewlyMarkedEdges;
		std::vector<Face*>			m_vNewlyMarkedFaces;
		std::vector<Volume*>		m_vNewlyMarkedVolumes;
//TODO	we need vNewlyMarkedFaces and vNewlyMarkedVolumes, too.
};

}//	end of namespace

#endif
