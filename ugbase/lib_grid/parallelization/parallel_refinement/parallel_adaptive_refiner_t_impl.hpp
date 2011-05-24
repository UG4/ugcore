// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 09.02.2011 (m,d,y)

#ifndef __H__UG__parallel_adaptive_refiner_t_impl__
#define __H__UG__parallel_adaptive_refiner_t_impl__

#include "parallel_adaptive_refiner_t.h"
#include "../util/compol_selection.h"

using namespace std;

namespace ug{

template <class TRefiner>
TParallelAdaptiveRefiner<TRefiner>::
TParallelAdaptiveRefiner(
		IRefinementCallback* refCallback) :
	BaseClass(refCallback),
	m_pDistGridMgr(NULL),
	m_pMG(NULL),
	m_bNewInterfaceVerticesMarked(false),
	m_bNewInterfaceEdgesMarked(false),
	m_bNewInterfaceFacesMarked(false),
	m_bNewInterfaceVolumesMarked(false)
{
}

template <class TRefiner>
TParallelAdaptiveRefiner<TRefiner>::
TParallelAdaptiveRefiner(
		DistributedGridManager& distGridMgr,
		IRefinementCallback* refCallback) :
	BaseClass(*distGridMgr.get_assigned_grid(), refCallback),
	m_pDistGridMgr(&distGridMgr),
	m_pMG(distGridMgr.get_assigned_grid()),
	m_bNewInterfaceVerticesMarked(false),
	m_bNewInterfaceEdgesMarked(false),
	m_bNewInterfaceFacesMarked(false),
	m_bNewInterfaceVolumesMarked(false)
{
}

template <class TRefiner>
TParallelAdaptiveRefiner<TRefiner>::
~TParallelAdaptiveRefiner()
{

}

template <class TRefiner>
void
TParallelAdaptiveRefiner<TRefiner>::
set_distributed_grid_manager(DistributedGridManager& distGridMgr)
{
	m_pDistGridMgr = &distGridMgr;
	m_pMG = distGridMgr.get_assigned_grid();
}

template <class TRefiner>
void
TParallelAdaptiveRefiner<TRefiner>::
clear_marks()
{
	BaseClass::clear_marks();
	m_bNewInterfaceVerticesMarked = false;
	m_bNewInterfaceEdgesMarked = false;
	m_bNewInterfaceFacesMarked = false;
	m_bNewInterfaceVolumesMarked = false;
}

template <class TRefiner>
void
TParallelAdaptiveRefiner<TRefiner>::
mark(VertexBase* v, RefinementMark refMark)
{
	if((!BaseClass::is_marked(v))
		&& (!m_pMG->has_children(v))
		&& m_pDistGridMgr->is_interface_element(v))
		m_bNewInterfaceVerticesMarked = true;
	BaseClass::mark(v, refMark);
}

template <class TRefiner>
void
TParallelAdaptiveRefiner<TRefiner>::
mark(EdgeBase* e, RefinementMark refMark)
{
	if((!BaseClass::is_marked(e))
		&& (!m_pMG->has_children(e))
		&& m_pDistGridMgr->is_interface_element(e))
		m_bNewInterfaceEdgesMarked = true;
	BaseClass::mark(e, refMark);
}

template <class TRefiner>
void
TParallelAdaptiveRefiner<TRefiner>::
mark(Face* f, RefinementMark refMark)
{
	if((!BaseClass::is_marked(f))
		&& (!m_pMG->has_children(f))
		&& m_pDistGridMgr->is_interface_element(f))
		m_bNewInterfaceFacesMarked = true;
	BaseClass::mark(f, refMark);
}

template <class TRefiner>
void
TParallelAdaptiveRefiner<TRefiner>::
mark(Volume* v, RefinementMark refMark)
{
	if((!BaseClass::is_marked(v))
		&& (!m_pMG->has_children(v))
		&& m_pDistGridMgr->is_interface_element(v))
		m_bNewInterfaceVolumesMarked = true;
	BaseClass::mark(v, refMark);
}

template <class TRefiner>
void
TParallelAdaptiveRefiner<TRefiner>::
collect_objects_for_refine()
{
//todo: This method could be improved.
//		In its current implementation a little too much
//		serial work is done.

//	the layoutmap is used for communication
	GridLayoutMap& layoutMap = m_pDistGridMgr->grid_layout_map();

//	first we'll call the base implementation
	while(1){
	//	we call collect_objects_for_refine in each iteration.
	//	This might be a bit of an overkill, since only a few normally
	//	have changed...
		BaseClass::collect_objects_for_refine();

	//	we now have to inform all processes whether interface elements
	//	were marked on any process.
		int newlyMarkedElems = 0;
		if(m_bNewInterfaceVerticesMarked ||
			m_bNewInterfaceEdgesMarked ||
			m_bNewInterfaceFacesMarked ||
			m_bNewInterfaceVolumesMarked)
		{
			newlyMarkedElems = 1;
		}

		int exchangeFlag;
		m_procCom.allreduce(&newlyMarkedElems, &exchangeFlag, 1,
							PCL_DT_INT, PCL_RO_LOR);

	//	before we continue we'll set all flags to false
		m_bNewInterfaceVerticesMarked = false;
		m_bNewInterfaceEdgesMarked = false;
		m_bNewInterfaceFacesMarked = false;
		m_bNewInterfaceVolumesMarked = false;

		if(exchangeFlag){

		//	we have to communicate the marks.
		//	do this by first gather selection at master nodes
		//	and then distribute them to slaves.
			ComPol_Selection<VertexLayout> compolSelVRT(BaseClass::m_selMarkedElements, true, false);
			ComPol_Selection<EdgeLayout> compolSelEDGE(BaseClass::m_selMarkedElements, true, false);
			ComPol_Selection<FaceLayout> compolSelFACE(BaseClass::m_selMarkedElements, true, false);

		//	send data SLAVE -> MASTER
			m_intfComVRT.exchange_data(layoutMap, INT_H_SLAVE, INT_H_MASTER,
										compolSelVRT);

			m_intfComEDGE.exchange_data(layoutMap, INT_H_SLAVE, INT_H_MASTER,
										compolSelEDGE);

			m_intfComFACE.exchange_data(layoutMap, INT_H_SLAVE, INT_H_MASTER,
										compolSelFACE);

			m_intfComVRT.communicate();
			m_intfComEDGE.communicate();
			m_intfComFACE.communicate();

		//	and now MASTER -> SLAVE (the selection has been adjusted on the fly)
			m_intfComVRT.exchange_data(layoutMap, INT_H_MASTER, INT_H_SLAVE,
										compolSelVRT);

			m_intfComEDGE.exchange_data(layoutMap, INT_H_MASTER, INT_H_SLAVE,
										compolSelEDGE);

			m_intfComFACE.exchange_data(layoutMap, INT_H_MASTER, INT_H_SLAVE,
										compolSelFACE);

			m_intfComVRT.communicate();
			m_intfComEDGE.communicate();
			m_intfComFACE.communicate();
		}
		else
			break;
	}
}


template <class TRefiner>
bool
TParallelAdaptiveRefiner<TRefiner>::
refinement_is_allowed(VertexBase* elem)
{
	return !m_pDistGridMgr->is_ghost(elem);
}

template <class TRefiner>
bool
TParallelAdaptiveRefiner<TRefiner>::
refinement_is_allowed(EdgeBase* elem)
{
	return !m_pDistGridMgr->is_ghost(elem);
}

template <class TRefiner>
bool
TParallelAdaptiveRefiner<TRefiner>::
refinement_is_allowed(Face* elem)
{
	return !m_pDistGridMgr->is_ghost(elem);
}

template <class TRefiner>
bool
TParallelAdaptiveRefiner<TRefiner>::
refinement_is_allowed(Volume* elem)
{
	return !m_pDistGridMgr->is_ghost(elem);
}

template <class TRefiner>
void
TParallelAdaptiveRefiner<TRefiner>::
pre_refine()
{
	m_pDistGridMgr->begin_ordered_element_insertion();
	BaseClass::pre_refine();
}

template <class TRefiner>
void
TParallelAdaptiveRefiner<TRefiner>::
post_refine()
{
	BaseClass::post_refine();
	m_pDistGridMgr->end_ordered_element_insertion();
}

template <class TRefiner>
void
TParallelAdaptiveRefiner<TRefiner>::
set_involved_processes(pcl::ProcessCommunicator com)
{
	m_procCom = com;
}

template <class TRefiner>
void
TParallelAdaptiveRefiner<TRefiner>::
refine()
{
	UG_ASSERT(m_pDistGridMgr, "a distributed grid manager has to be assigned");
	if(!m_pDistGridMgr){
		throw(UGError("No distributed grid manager assigned."));
	}

	BaseClass::refine();

/*
//	DEBUG ONLY
//	Make sure that the interfaces and layouts are fine.
	pcl::ParallelCommunicator<VertexLayout::LevelLayout> com;
	GridLayoutMap& layoutMap = m_pDistGridMgr->grid_layout_map();

	UG_LOG("\nTesting horizontal layouts...\n");
	{
		VertexLayout& masterLayout = layoutMap.get_layout<VertexBase>(INT_MASTER);
		VertexLayout& slaveLayout = layoutMap.get_layout<VertexBase>(INT_H_SLAVE);
		for(size_t i = 0; i < m_pMG->num_levels(); ++i){
			UG_LOG("Testing VertexLayout on level " << i << ":" << endl);
			pcl::TestLayout(com, masterLayout.layout_on_level(i),
					slaveLayout.layout_on_level(i), true);
		}
	}

	UG_LOG("\nTesting vertical layouts...\n");
	{
		VertexLayout& masterLayout = layoutMap.get_layout<VertexBase>(INT_V_MASTER);
		VertexLayout& slaveLayout = layoutMap.get_layout<VertexBase>(INT_V_SLAVE);
		for(size_t i = 0; i < m_pMG->num_levels(); ++i){
			UG_LOG("Testing VertexLayout on level " << i << ":" << endl);
			pcl::TestLayout(com, masterLayout.layout_on_level(i),
					slaveLayout.layout_on_level(i), true);
		}
	}

	UG_LOG("\nTesting virtual layouts...\n");
	{
		VertexLayout& masterLayout = layoutMap.get_layout<VertexBase>(INT_VIRTUAL_MASTER);
		VertexLayout& slaveLayout = layoutMap.get_layout<VertexBase>(INT_VIRTUAL_SLAVE);
		for(size_t i = 0; i < m_pMG->num_levels(); ++i){
			UG_LOG("Testing VerticalVertexLayout on level " << i << ":" << endl);
			pcl::TestLayout(com, masterLayout.layout_on_level(i),
					slaveLayout.layout_on_level(i), true);
		}
	}
*/
}

}// end of namespace

#endif
