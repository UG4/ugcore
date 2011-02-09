// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 09.02.2011 (m,d,y)
 
#include "parallel_hanging_node_refiner_multi_grid.h"

namespace ug{

ParallelHangingNodeRefiner_MultiGrid::
ParallelHangingNodeRefiner_MultiGrid(
		DistributedGridManager& distGridMgr) :
	m_distGridMgr(distGridMgr),
	m_bNewInterfaceEdgesMarked(false),
	m_bNewInterfaceFacesMarked(false),
	m_bNewInterfaceVolumesMarked(false)
{
}

ParallelHangingNodeRefiner_MultiGrid::
~ParallelHangingNodeRefiner_MultiGrid()
{

}

void ParallelHangingNodeRefiner_MultiGrid::
clear_marks()
{
	HangingNodeRefiner_MultiGrid::clear_marks();
	m_bNewInterfaceEdgesMarked = false;
	m_bNewInterfaceFacesMarked = false;
	m_bNewInterfaceVolumesMarked = false;
}

void ParallelHangingNodeRefiner_MultiGrid::
mark_for_refinement(EdgeBase* e)
{
	HangingNodeRefiner_MultiGrid::mark_for_refinement(e);
	if(m_distGridMgr.is_interface_element(e))
		m_bNewInterfaceEdgesMarked = true;
}

void ParallelHangingNodeRefiner_MultiGrid::
mark_for_refinement(Face* f)
{
	HangingNodeRefiner_MultiGrid::mark_for_refinement(f);
	if(m_distGridMgr.is_interface_element(f))
		m_bNewInterfaceFacesMarked = true;
}

void ParallelHangingNodeRefiner_MultiGrid::
mark_for_refinement(Volume* v)
{
	HangingNodeRefiner_MultiGrid::mark_for_refinement(v);
	if(m_distGridMgr.is_interface_element(v))
		m_bNewInterfaceVolumesMarked = true;
}

void ParallelHangingNodeRefiner_MultiGrid::
collect_objects_for_refine()
{
//todo: This method could be improved.
//		In its current implementation a little too much
//		serial work is done.

//	first we'll call the base implementation
	while(1){
		HangingNodeRefiner_MultiGrid::collect_objects_for_refine();

	//	we now have to inform all processes whether interface elements
	//	were marked on any process.
		int newlyMarkedElems = 0;
		if(m_bNewInterfaceEdgesMarked ||
			m_bNewInterfaceFacesMarked ||
			m_bNewInterfaceVolumesMarked)
		{
			newlyMarkedElems = 1;
		}

		int exchangeFlag;
		m_procCom.allreduce(&newlyMaredElems, &exchangeFlag, 1,
							PCL_DT_INT, PCL_RO_LOR);

		if(exchangeFlag){
		//	we have to communicate the marks
		}
		else
			break;
	}
}

void ParallelHangingNodeRefiner_MultiGrid::
pre_refine()
{
	m_distGridMgr.begin_ordered_element_insertion();
}

void ParallelHangingNodeRefiner_MultiGrid::
post_refine()
{
	m_distGridMgr.end_ordered_element_insertion();
}

void ParallelHangingNodeRefiner_MultiGrid::
set_involved_processes(pcl::ProcessCommunicator com)
{
	m_procCom = com;
}

}// end of namespace
