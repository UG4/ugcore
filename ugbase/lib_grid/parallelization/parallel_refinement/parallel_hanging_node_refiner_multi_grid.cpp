// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 09.02.2011 (m,d,y)
 
#include "parallel_hanging_node_refiner_multi_grid.h"
#include "../util/compol_selection.h"
#include "parallel_hnode_adjuster.h"
#include "lib_grid/algorithms/debug_util.h"

namespace ug{

ParallelHangingNodeRefiner_MultiGrid::
ParallelHangingNodeRefiner_MultiGrid(
		IRefinementCallback* refCallback) :
	BaseClass(refCallback),
	m_pDistGridMgr(NULL),
	m_pMG(NULL)
{
	add_ref_mark_adjuster(ParallelHNodeAdjuster::create());
}

ParallelHangingNodeRefiner_MultiGrid::
ParallelHangingNodeRefiner_MultiGrid(
		DistributedGridManager& distGridMgr,
		IRefinementCallback* refCallback) :
	BaseClass(*distGridMgr.get_assigned_grid(), refCallback),
	m_pDistGridMgr(&distGridMgr),
	m_pMG(distGridMgr.get_assigned_grid())
{
	add_ref_mark_adjuster(ParallelHNodeAdjuster::create());
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

bool ParallelHangingNodeRefiner_MultiGrid::
continue_collect_objects_for_refine(bool continueRequired)
{
	return pcl::OneProcTrue(continueRequired);
}


static void ReplaceByNormal(MultiGrid& mg, VertexBase* v)
{
	mg.create_and_replace<Vertex>(v);
}

static void ReplaceByConstrained(MultiGrid& mg, VertexBase* v)
{
	mg.create_and_replace<ConstrainedVertex>(v);
}

static void ReplaceByConstraining(MultiGrid&, VertexBase*)
{
	UG_THROW("Can't convert vertex to constraining-vertex! "
			"(Constraining vertices don't exist!");
}

static void ReplaceByNormal(MultiGrid& mg, EdgeBase* e)
{
	mg.create_and_replace<Edge>(e);
}

static void ReplaceByConstrained(MultiGrid& mg, EdgeBase* e)
{
	mg.create_and_replace<ConstrainedEdge>(e);
}

static void ReplaceByConstraining(MultiGrid& mg, EdgeBase* e)
{
	mg.create_and_replace<ConstrainingEdge>(e);
}

static void ReplaceByNormal(MultiGrid& mg, Face* f)
{
	if(f->num_vertices() == 3)
		mg.create_and_replace<Triangle>(f);
	else{
		UG_ASSERT(f->num_vertices() == 4,
				"Only triangles and quatrilaterals currently supported");
		mg.create_and_replace<Quadrilateral>(f);
	}
}

static void ReplaceByConstrained(MultiGrid& mg, Face* f)
{
	if(f->num_vertices() == 3)
		mg.create_and_replace<ConstrainedTriangle>(f);
	else{
		UG_ASSERT(f->num_vertices() == 4,
				"Only triangles and quatrilaterals currently supported");
		mg.create_and_replace<ConstrainedQuadrilateral>(f);
	}
}

static void ReplaceByConstraining(MultiGrid& mg, Face* f)
{
	if(f->num_vertices() == 3)
		mg.create_and_replace<ConstrainingTriangle>(f);
	else{
		UG_ASSERT(f->num_vertices() == 4,
				"Only triangles and quatrilaterals currently supported");
		mg.create_and_replace<ConstrainingQuadrilateral>(f);
	}
}

/**	During collection all elements which change their constrained/constraining type
 * push their interface index and an indicator for their new type to the buffer.
 *
 * During extraction, the local type of the element at a received interface entry
 * is checked. If the type is different, the element is immediately replaced in
 * the underlying grid. This should work fine with the distributed grid manager,
 * which updates interfaces on the fly.
 *
 * This ComPol is intended for communication from vertical slaves to vertical
 * masters only. Changes are only applied to ghost-elements.
 *
 * This ComPol may currently only be used for vertices, edges and faces.
 */
template <class TLayout>
class ComPol_AdjustType : public pcl::ICommunicationPolicy<TLayout>
{
	public:
		typedef TLayout								Layout;
		typedef typename Layout::Type				GeomObj;
		typedef typename Layout::Element			Element;
		typedef typename Layout::Interface			Interface;
		typedef typename Interface::const_iterator	InterfaceIter;

		enum ConversionTypes{
			CT_IGNORE = 0,
			CT_TO_NORMAL = 1,
			CT_TO_CONSTRAINED = 2,
			CT_TO_CONSTRAINING = 3
		};

		ComPol_AdjustType(Selector& sel, DistributedGridManager& distGridMgr)
			 :	m_sel(sel), m_distGridMgr(distGridMgr)
		{}

		virtual ~ComPol_AdjustType()		{}

		virtual int
		get_required_buffer_size(const Interface& interface)		{return -1;}

	///	writes the selection states of the interface entries
		virtual bool
		collect(ug::BinaryBuffer& buff, const Interface& interface)
		{
		//	search for entries which changed their constrained/constraining status
			UG_ASSERT(m_distGridMgr.get_assigned_grid(),
					"The distributed grid manager has to operate on a valid grid!");
			MultiGrid& mg = *m_distGridMgr.get_assigned_grid();

			int counter = 0;
			for(InterfaceIter iter = interface.begin();
				iter != interface.end(); ++iter, ++counter)
			{
				byte mark = CT_IGNORE;
				Element elem = interface.get_element(iter);

				if(m_sel.is_selected(elem)){
					if(elem->is_constraining()){
						if((m_sel.get_selection_status(elem)
							& HangingNodeRefinerBase::HNRM_CONSTRAINED) == 0)
						{
							mark = CT_TO_NORMAL;
						}
					}
					else if(elem->is_constrained()){
						if((m_sel.get_selection_status(elem)
							& HangingNodeRefinerBase::HNRM_CONSTRAINED))
						{
							mark = CT_TO_CONSTRAINING;
						}
					}
					else{
						if((m_sel.get_selection_status(elem)
							& HangingNodeRefinerBase::HNRM_CONSTRAINED))
						{
							mark = CT_TO_CONSTRAINING;
						}
					}
				}
				else{
					if(elem->is_constrained()){
						GeometricObject* go = mg.get_parent(elem);
						if(go && m_sel.is_selected(go)){
							if((m_sel.get_selection_status(go)
								& HangingNodeRefinerBase::HNRM_CONSTRAINED) == 0)
							{
								mark = CT_TO_NORMAL;
							}
						}
					}
				}

				if(mark != CT_IGNORE){
					buff.write((char*)&counter, sizeof(int));
					buff.write((char*)&mark, sizeof(byte));
				}
			}

			int eof = -1;
			buff.write((char*)&eof, sizeof(int));
			return true;
		}

	///	reads marks from the given stream
		virtual bool
		extract(ug::BinaryBuffer& buff, const Interface& interface)
		{
		//	search for entries which changed their constrained/constraining status
			UG_ASSERT(m_distGridMgr.get_assigned_grid(),
					"The distributed grid manager has to operate on a valid grid!");
			MultiGrid& mg = *m_distGridMgr.get_assigned_grid();

			int counter = 0;
			InterfaceIter iter = interface.begin();
			while(1){
				int index;
				buff.read((char*)&index, sizeof(int));
				if(index == -1)
					break;

				byte mark;
				buff.read((char*)&mark, sizeof(byte));

				while((counter < index) && (iter != interface.end())){
					++iter;
					++counter;
				}

				if(iter == interface.end()){
					UG_LOG("WARNING: Unexpectedly reached end of interface\n");
				}

				Element elem = interface.get_element(iter);

			//	if the element is also contained in a horizontal interface, we
			//	don't have (and indeed must not) convert the element during this
			//	communication step.
				if(m_distGridMgr.is_in_horizontal_interface(elem))
					continue;

				switch(mark){
					case CT_TO_NORMAL:
						if(elem->is_constraining() || elem->is_constrained()){
//							UG_DLOG(LIB_GRID, 2, "ParHNodeRef: replacing with normal element at "
//									<< GetGeometricObjectCenter(mg, elem) << ".\n");
							ReplaceByNormal(mg, elem);
						}
						break;
					case CT_TO_CONSTRAINING:
						if(!elem->is_constraining()){
//							UG_DLOG(LIB_GRID, 2, "ParHNodeRef: replacing with constraining element at "
//									<< GetGeometricObjectCenter(mg, elem) << ".\n");
							ReplaceByConstraining(mg, elem);
						}
						break;
					case CT_TO_CONSTRAINED:
						if(!elem->is_constrained()){
//							UG_DLOG(LIB_GRID, 2, "ParHNodeRef: replacing with constrained element at "
//									<< GetGeometricObjectCenter(mg, elem) << ".\n");
							ReplaceByConstrained(mg, elem);
						}
						break;
				}
			}
			return true;
		}

	private:
		Selector& 				m_sel;
		DistributedGridManager&	m_distGridMgr;
};


void ParallelHangingNodeRefiner_MultiGrid::
assign_hnode_marks()
{
	UG_DLOG(LIB_GRID, 1, "ParHNodeRef-start: assign_hnode_marks\n");
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


//	until now we only communicated hnode-marks over horizontal interfaces.
//	However, while we won't refine ghosts (vertical masters which do not lie
//	in any horizontal interface) we have to adjust their type (constrained / constraining)
//	so that they match the type of their vertical slaves. While this is not really
//	important for most uses, differing types of ghosts cause problems during
//	redistribution, which possibly causes issues in methods executed afterwards.
//	Note: we adjust types of ghosts before refinement here. This is not a problem,
//	since it only involves replacement, which can be handled by the distributed
//	grid manager even outside of begin_ordered_element_insertion/end_...

	ComPol_AdjustType<VertexLayout> compolAdjustVRT(BaseClass::m_selMarkedElements,
													*m_pDistGridMgr);
	ComPol_AdjustType<EdgeLayout> compolAdjustEDGE(BaseClass::m_selMarkedElements,
												   *m_pDistGridMgr);
	ComPol_AdjustType<FaceLayout> compolAdjustFACE(BaseClass::m_selMarkedElements,
												   *m_pDistGridMgr);

	m_intfComVRT.exchange_data(layoutMap, INT_V_SLAVE, INT_V_MASTER,
							   compolAdjustVRT);
	m_intfComEDGE.exchange_data(layoutMap, INT_V_SLAVE, INT_V_MASTER,
								compolAdjustEDGE);
	m_intfComFACE.exchange_data(layoutMap, INT_V_SLAVE, INT_V_MASTER,
								compolAdjustFACE);

	m_intfComVRT.communicate();
	m_intfComEDGE.communicate();
	m_intfComFACE.communicate();

	UG_DLOG(LIB_GRID, 1, "ParHNodeRef-stop: assign_hnode_marks\n");
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
		typedef TLayout								Layout;
		typedef typename Layout::Type				GeomObj;
		typedef typename Layout::Element			Element;
		typedef typename Layout::Interface			Interface;
		typedef typename Interface::const_iterator	InterfaceIter;

		ComPol_BroadcastCoarsenMarks(Selector& sel)
			 :	m_sel(sel)
		{}

		virtual int
		get_required_buffer_size(const Interface& interface)
		{return interface.size() * sizeof(byte);}

	///	writes writes the selection states of the interface entries
		virtual bool
		collect(ug::BinaryBuffer& buff, const Interface& interface)
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
		extract(ug::BinaryBuffer& buff, const Interface& interface)
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
