//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m04 d14

#ifndef __H__LIB_GRID__PARALLEL_GLOBAL_REFINER_T_IMPL__
#define __H__LIB_GRID__PARALLEL_GLOBAL_REFINER_T_IMPL__

#include <vector>
#include "parallel_global_refiner_t.h"
#include "pcl/pcl.h"

namespace ug
{

template <class TRefiner>
TParallelGlobalRefiner<TRefiner>::
TParallelGlobalRefiner(DistributedGridManager& distGridMgr) :
	TRefiner(*distGridMgr.get_assigned_grid()),
	m_distGridMgr(distGridMgr)
{
}

template <class TRefiner>
TParallelGlobalRefiner<TRefiner>::
~TParallelGlobalRefiner()
{
}

template <class TRefiner>
bool
TParallelGlobalRefiner<TRefiner>::
refinement_is_allowed(Vertex* elem)
{
	return !m_distGridMgr.is_ghost(elem);
}

template <class TRefiner>
bool TParallelGlobalRefiner<TRefiner>::
refinement_is_allowed(EdgeBase* elem)
{
	return !m_distGridMgr.is_ghost(elem);
}

template <class TRefiner>
bool TParallelGlobalRefiner<TRefiner>::
refinement_is_allowed(Face* elem)
{
	return !m_distGridMgr.is_ghost(elem);
}

template <class TRefiner>
bool TParallelGlobalRefiner<TRefiner>::
refinement_is_allowed(Volume* elem)
{
	return !m_distGridMgr.is_ghost(elem);
}

template <class TRefiner>
void TParallelGlobalRefiner<TRefiner>::
refinement_step_begins()
{
	m_distGridMgr.begin_ordered_element_insertion();
}

template <class TRefiner>
void TParallelGlobalRefiner<TRefiner>::
refinement_step_ends()
{
	m_distGridMgr.end_ordered_element_insertion();
}

}//	end of namespace

#endif
