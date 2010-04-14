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
	m_distGridMgr(distGridMgr),
	GlobalMultiGridRefiner(*distGridMgr.get_assigned_grid())
{
}

ParallelGlobalMultiGridRefiner::~ParallelGlobalMultiGridRefiner()
{
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
