//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m04 d14

#include <vector>
#include "parallel_global_fractured_media_refiner.h"
#include "pcl/pcl.h"
#include "lib_grid/parallelization/util/compol_boolmarker.h"

namespace ug
{

ParallelGlobalFracturedMediaRefiner::
ParallelGlobalFracturedMediaRefiner(DistributedGridManager& distGridMgr) :
	GlobalFracturedMediaRefiner(*distGridMgr.get_assigned_grid()),
	m_distGridMgr(distGridMgr)
{
}


ParallelGlobalFracturedMediaRefiner::
~ParallelGlobalFracturedMediaRefiner()
{
}


bool
ParallelGlobalFracturedMediaRefiner::
refinement_is_allowed(Vertex* elem)
{
	return !m_distGridMgr.is_ghost(elem);
}


bool ParallelGlobalFracturedMediaRefiner::
refinement_is_allowed(Edge* elem)
{
	return !m_distGridMgr.is_ghost(elem);
}


bool ParallelGlobalFracturedMediaRefiner::
refinement_is_allowed(Face* elem)
{
	return !m_distGridMgr.is_ghost(elem);
}


bool ParallelGlobalFracturedMediaRefiner::
refinement_is_allowed(Volume* elem)
{
	return !m_distGridMgr.is_ghost(elem);
}


void ParallelGlobalFracturedMediaRefiner::
refinement_step_begins()
{
	m_distGridMgr.begin_ordered_element_insertion();
}


void ParallelGlobalFracturedMediaRefiner::
refinement_step_ends()
{
	m_distGridMgr.end_ordered_element_insertion();
}


void ParallelGlobalFracturedMediaRefiner::
communicate_marks(BoolMarker& marker)
{
	GridLayoutMap& layoutMap = m_distGridMgr.grid_layout_map();

//	we have to communicate side marks. In 3d we also have to communicate edge marks.
//	we'll simply communicate edge and face marks. This is no overhead,
//	since no face interfaces don't exist in 2d anyways. The 1d case is ignored.
	ComPol_BoolMarker_AddMarks<EdgeLayout> compolMarkerEDGE(marker);
	ComPol_BoolMarker_AddMarks<FaceLayout> compolMarkerFACE(marker);

//	SLAVE->MASTER
	m_intfComEDGE.exchange_data(layoutMap, INT_H_SLAVE, INT_H_MASTER,
								compolMarkerEDGE);
	m_intfComFACE.exchange_data(layoutMap, INT_H_SLAVE, INT_H_MASTER,
								compolMarkerFACE);

	m_intfComEDGE.communicate();
	m_intfComFACE.communicate();

//	MASTER->SLAVE
	m_intfComEDGE.exchange_data(layoutMap, INT_H_MASTER, INT_H_SLAVE,
								compolMarkerEDGE);
	m_intfComFACE.exchange_data(layoutMap, INT_H_MASTER, INT_H_SLAVE,
								compolMarkerFACE);

	m_intfComEDGE.communicate();
	m_intfComFACE.communicate();
}

}//	end of namespace
