//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m04 d14

#include <vector>
#include "parallel_global_multi_grid_refiner.h"
#include "pcl/pcl.h"

using namespace std;

namespace ug
{

ParallelGlobalMultiGridRefiner::
ParallelGlobalMultiGridRefiner(DistributedGridManager& distGridMgr) :
	GlobalMultiGridRefiner(*distGridMgr.get_assigned_grid()),
	m_distGridMgr(distGridMgr)
{
}

ParallelGlobalMultiGridRefiner::~ParallelGlobalMultiGridRefiner()
{
}

bool ParallelGlobalMultiGridRefiner::
refinement_is_allowed(VertexBase* elem)
{
	UG_LOG("vrt");
	return !m_distGridMgr.is_ghost(elem);
}

bool ParallelGlobalMultiGridRefiner::
refinement_is_allowed(EdgeBase* elem)
{
	return !m_distGridMgr.is_ghost(elem);
}

bool ParallelGlobalMultiGridRefiner::
refinement_is_allowed(Face* elem)
{
	UG_LOG("f");
	return !m_distGridMgr.is_ghost(elem);
}

bool ParallelGlobalMultiGridRefiner::
refinement_is_allowed(Volume* elem)
{
	return !m_distGridMgr.is_ghost(elem);
}
		
void ParallelGlobalMultiGridRefiner::
refinement_step_begins()
{
	m_distGridMgr.begin_ordered_element_insertion();
}

void ParallelGlobalMultiGridRefiner::
refinement_step_ends()
{
	m_distGridMgr.end_ordered_element_insertion();
}

}//	end of namespace
