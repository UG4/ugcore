// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 09.02.2011 (m,d,y)
 
#include "parallel_hanging_node_refiner_multi_grid.h"
#include "../util/compol_selection.h"

namespace ug{

ParallelHangingNodeRefiner_MultiGrid::
ParallelHangingNodeRefiner_MultiGrid(
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

ParallelHangingNodeRefiner_MultiGrid::
ParallelHangingNodeRefiner_MultiGrid(
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

ParallelHangingNodeRefiner_MultiGrid::
~ParallelHangingNodeRefiner_MultiGrid()
{

}

void
ParallelHangingNodeRefiner_MultiGrid::
set_distributed_grid_manager(DistributedGridManager& distGridMgr)
{
	m_pDistGridMgr = &distGridMgr;
	m_pMG = distGridMgr.get_assigned_grid();
}

void
ParallelHangingNodeRefiner_MultiGrid::
clear_marks()
{
	BaseClass::clear_marks();
	m_bNewInterfaceVerticesMarked = false;
	m_bNewInterfaceEdgesMarked = false;
	m_bNewInterfaceFacesMarked = false;
	m_bNewInterfaceVolumesMarked = false;
}

bool
ParallelHangingNodeRefiner_MultiGrid::
mark(VertexBase* v, RefinementMark refMark)
{
	RefinementMark oldMark = BaseClass::get_mark(v);
	if(BaseClass::mark(v, refMark)){
		if((refMark != oldMark)
		  && m_pDistGridMgr->is_interface_element(v))
			m_bNewInterfaceVerticesMarked = true;
		return true;
	}
	return false;
}

bool
ParallelHangingNodeRefiner_MultiGrid::
mark(EdgeBase* e, RefinementMark refMark)
{
	RefinementMark oldMark = BaseClass::get_mark(e);
	if(BaseClass::mark(e, refMark)){
		if((refMark != oldMark)
		  && m_pDistGridMgr->is_interface_element(e))
			m_bNewInterfaceEdgesMarked = true;
		return true;
	}
	return false;
}

bool
ParallelHangingNodeRefiner_MultiGrid::
mark(Face* f, RefinementMark refMark)
{
	RefinementMark oldMark = BaseClass::get_mark(f);
	if(BaseClass::mark(f, refMark)){
		if((refMark != oldMark)
		  && m_pDistGridMgr->is_interface_element(f))
			m_bNewInterfaceFacesMarked = true;
		return true;
	}
	return false;
}

bool
ParallelHangingNodeRefiner_MultiGrid::
mark(Volume* v, RefinementMark refMark)
{
	RefinementMark oldMark = BaseClass::get_mark(v);
	if(BaseClass::mark(v, refMark)){
		if((refMark != oldMark)
		  && m_pDistGridMgr->is_interface_element(v))
			m_bNewInterfaceVolumesMarked = true;
		return true;
	}
	return false;
}

void
ParallelHangingNodeRefiner_MultiGrid::
collect_objects_for_refine()
{
//todo: This method could be improved.
//		In its current implementation a little too much
//		serial work is done.
	UG_DLOG(LIB_GRID, 1, "  collect_objects_for_refine started for parallel multi-grid...\n");

//	the layoutmap is used for communication
	GridLayoutMap& layoutMap = m_pDistGridMgr->grid_layout_map();

//	first we'll call the base implementation
	while(1){
		UG_DLOG(LIB_GRID, 2, "  new parallel iteration in collect_objects_for_refine\n");
	//	we call collect_objects_for_refine in each iteration.
	//	This might be a bit of an overkill, since only a few normally
	//	have changed...
		BaseClass::collect_objects_for_refine();

	//	we now have to inform all processes whether interface elements
	//	were marked on any process.
		UG_DLOG(LIB_GRID, 2, "  exchanging flag for newly marked elements (allreduce)\n");
		bool newlyMarkedElems = m_bNewInterfaceVerticesMarked ||
								m_bNewInterfaceEdgesMarked ||
								m_bNewInterfaceFacesMarked ||
								m_bNewInterfaceVolumesMarked;

		bool exchangeFlag = pcl::OneProcTrue(newlyMarkedElems);

	//	before we continue we'll set all flags to false
		m_bNewInterfaceVerticesMarked = false;
		m_bNewInterfaceEdgesMarked = false;
		m_bNewInterfaceFacesMarked = false;
		m_bNewInterfaceVolumesMarked = false;

		if(exchangeFlag){
			UG_DLOG(LIB_GRID, 2, "  there are newly marked interface elements...\n");
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
		else{
			UG_DLOG(LIB_GRID, 2, "    there are no newly marked interface elements...\n");
			UG_DLOG(LIB_GRID, 2, "    leaving parallel collect_objects_for_refine iteration...\n");
			break;
		}
	}

	UG_DLOG(LIB_GRID, 1, "  collect_objects_for_refine done for parallel multi-grid...\n");
}

void ParallelHangingNodeRefiner_MultiGrid::
assign_hnode_marks()
{
//	call base implementation
	BaseClass::assign_hnode_marks();

//	copy the hnode mark.
//	note that we're enabling the mark, but never disable it.
//	first we enable it at the master if one of the slaves is enabled,
//	then we enable it at the slaves, if the master was enabled.

	GridLayoutMap& layoutMap = m_pDistGridMgr->grid_layout_map();

	ComPol_EnableSelectionStateBits<EdgeLayout> compolEDGE(BaseClass::m_selMarkedElements,
														 BaseClass::HNRM_CONSTRAINED);
	ComPol_EnableSelectionStateBits<FaceLayout> compolFACE(BaseClass::m_selMarkedElements,
														 BaseClass::HNRM_CONSTRAINED);

	m_intfComEDGE.exchange_data(layoutMap, INT_H_SLAVE, INT_H_MASTER,
								compolEDGE);

	m_intfComFACE.exchange_data(layoutMap, INT_H_SLAVE, INT_H_MASTER,
								compolFACE);

	m_intfComEDGE.communicate();
	m_intfComFACE.communicate();

	m_intfComEDGE.exchange_data(layoutMap, INT_H_MASTER, INT_H_SLAVE,
								compolEDGE);

	m_intfComFACE.exchange_data(layoutMap, INT_H_MASTER, INT_H_SLAVE,
								compolFACE);

	m_intfComEDGE.communicate();
	m_intfComFACE.communicate();
}

bool
ParallelHangingNodeRefiner_MultiGrid::
refinement_is_allowed(VertexBase* elem)
{
	return (!m_pDistGridMgr->is_ghost(elem))
			&& BaseClass::refinement_is_allowed(elem);
}

bool
ParallelHangingNodeRefiner_MultiGrid::
refinement_is_allowed(EdgeBase* elem)
{
	return (!m_pDistGridMgr->is_ghost(elem))
			&& BaseClass::refinement_is_allowed(elem);
}

bool
ParallelHangingNodeRefiner_MultiGrid::
refinement_is_allowed(Face* elem)
{
	return (!m_pDistGridMgr->is_ghost(elem))
			&& BaseClass::refinement_is_allowed(elem);
}

bool
ParallelHangingNodeRefiner_MultiGrid::
refinement_is_allowed(Volume* elem)
{
	return (!m_pDistGridMgr->is_ghost(elem))
			&& BaseClass::refinement_is_allowed(elem);
}

void
ParallelHangingNodeRefiner_MultiGrid::
pre_refine()
{
	m_pDistGridMgr->begin_ordered_element_insertion();
	BaseClass::pre_refine();
}

void
ParallelHangingNodeRefiner_MultiGrid::
post_refine()
{
	BaseClass::post_refine();
	m_pDistGridMgr->end_ordered_element_insertion();
}

void
ParallelHangingNodeRefiner_MultiGrid::
pre_coarsen()
{
	m_pDistGridMgr->begin_element_deletion();
	BaseClass::pre_coarsen();
}

void
ParallelHangingNodeRefiner_MultiGrid::
post_coarsen()
{
	BaseClass::post_coarsen();
	m_pDistGridMgr->end_element_deletion();
}

void
ParallelHangingNodeRefiner_MultiGrid::
set_involved_processes(pcl::ProcessCommunicator com)
{
	m_procCom = com;
}


bool
ParallelHangingNodeRefiner_MultiGrid::
continue_collect_objects_for_coarsen(bool continueRequired)
{
	return (bool)m_procCom.allreduce((int)continueRequired, PCL_RO_LOR);
}


template <class TLayout>
class ComPol_BroadcastCoarsenMarks : public pcl::ICommunicationPolicy<TLayout>
{
	public:
		typedef TLayout							Layout;
		typedef typename Layout::Type			GeomObj;
		typedef typename Layout::Element		Element;
		typedef typename Layout::Interface		Interface;
		typedef typename Interface::iterator	InterfaceIter;

		ComPol_BroadcastCoarsenMarks(Selector& sel)
			 :	m_sel(sel)
		{}

		virtual int
		get_required_buffer_size(Interface& interface)		{return interface.size() * sizeof(byte);}

	///	writes writes the selection states of the interface entries
		virtual bool
		collect(ug::BinaryBuffer& buff, Interface& interface)
		{
		//	write the entry indices of marked elements.
			for(InterfaceIter iter = interface.begin();
				iter != interface.end(); ++iter)
			{
				Element elem = interface.get_element(iter);
				byte refMark = m_sel.get_selection_status(elem);
				buff.write((char*)&refMark, sizeof(byte));
			}

			return true;
		}

	///	reads marks from the given stream
		virtual bool
		extract(ug::BinaryBuffer& buff, Interface& interface)
		{
			byte val;
			for(InterfaceIter iter = interface.begin();
				iter != interface.end(); ++iter)
			{
				Element elem = interface.get_element(iter);
				buff.read((char*)&val, sizeof(byte));

			//	check the current status and adjust the mark accordingly
				byte curVal = m_sel.get_selection_status(elem);

				if(curVal != val)
					m_sel.select(elem, ParallelHangingNodeRefiner_MultiGrid::
													HNCM_SOME_NBRS_SELECTED);
			}
			return true;
		}

		Selector& m_sel;
};

void
ParallelHangingNodeRefiner_MultiGrid::
broadcast_vertex_coarsen_marks()
{
//	the broadcast has to be performed with some special operations on the
//	marks:
//	- if one mark equals HNCM_SOME_NBRS_SELECTED then all have to be set to
//		HNCM_SOME_NBRS_SELECTED.
//	- if on contains HNCM_ALL_NBRS_SELECTED and another contains another mark,
//		then all have to be set to HNCM_SOME_NBRS_SELECTED.
//	- else the mark stays as it is.
//
//	we'll collect the marks at the master-nodes, adjust the and distribute them
//	to the associated slaves afterwards.
	ComPol_BroadcastCoarsenMarks<VertexLayout>	comPol(get_refmark_selector());

	GridLayoutMap& layoutMap = m_pDistGridMgr->grid_layout_map();

	m_intfComVRT.exchange_data(layoutMap, INT_H_SLAVE, INT_H_MASTER, comPol);
	m_intfComVRT.communicate();

	m_intfComVRT.exchange_data(layoutMap, INT_H_MASTER, INT_H_SLAVE, comPol);
	m_intfComVRT.communicate();
}

void
ParallelHangingNodeRefiner_MultiGrid::
broadcast_edge_coarsen_marks()
{
//	the broadcast has to be performed with some special operations on the
//	marks:
//	- if one mark equals HNCM_SOME_NBRS_SELECTED then all have to be set to
//		HNCM_SOME_NBRS_SELECTED.
//	- if on contains HNCM_ALL_NBRS_SELECTED and another contains another mark,
//		then all have to be set to HNCM_SOME_NBRS_SELECTED.
//	- else the mark stays as it is.
//
//	we'll collect the marks at the master-nodes, adjust the and distribute them
//	to the associated slaves afterwards.
	ComPol_BroadcastCoarsenMarks<EdgeLayout>	comPol(get_refmark_selector());

	GridLayoutMap& layoutMap = m_pDistGridMgr->grid_layout_map();

	m_intfComEDGE.exchange_data(layoutMap, INT_H_SLAVE, INT_H_MASTER, comPol);
	m_intfComEDGE.communicate();

	m_intfComEDGE.exchange_data(layoutMap, INT_H_MASTER, INT_H_SLAVE, comPol);
	m_intfComEDGE.communicate();
}

void
ParallelHangingNodeRefiner_MultiGrid::
broadcast_face_coarsen_marks()
{
//	the broadcast has to be performed with some special operations on the
//	marks:
//	- if one mark equals HNCM_SOME_NBRS_SELECTED then all have to be set to
//		HNCM_SOME_NBRS_SELECTED.
//	- if on contains HNCM_ALL_NBRS_SELECTED and another contains another mark,
//		then all have to be set to HNCM_SOME_NBRS_SELECTED.
//	- else the mark stays as it is.
//
//	we'll collect the marks at the master-nodes, adjust the and distribute them
//	to the associated slaves afterwards.
	ComPol_BroadcastCoarsenMarks<FaceLayout>	comPol(get_refmark_selector());

	GridLayoutMap& layoutMap = m_pDistGridMgr->grid_layout_map();

	m_intfComFACE.exchange_data(layoutMap, INT_H_SLAVE, INT_H_MASTER, comPol);
	m_intfComFACE.communicate();

	m_intfComFACE.exchange_data(layoutMap, INT_H_MASTER, INT_H_SLAVE, comPol);
	m_intfComFACE.communicate();
}

bool
ParallelHangingNodeRefiner_MultiGrid::
contains_edges()
{
	bool containsEdges = m_pMG->num<EdgeBase>() > 0;
	return (bool)m_procCom.allreduce((int)containsEdges, PCL_RO_LOR);
}

bool
ParallelHangingNodeRefiner_MultiGrid::
contains_faces()
{
	bool containsFaces = m_pMG->num<Face>() > 0;
	return (bool)m_procCom.allreduce((int)containsFaces, PCL_RO_LOR);
}

bool
ParallelHangingNodeRefiner_MultiGrid::
contains_volumes()
{
	bool containsVolume = m_pMG->num<Volume>() > 0;
	return (bool)m_procCom.allreduce((int)containsVolume, PCL_RO_LOR);
}

/*
//	DEBUG ONLY
//	Make sure that the interfaces and layouts are fine.
	pcl::InterfaceCommunicator<VertexLayout::LevelLayout> com;
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

}// end of namespace
