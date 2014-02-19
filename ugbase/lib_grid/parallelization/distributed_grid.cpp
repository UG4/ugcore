// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y09 m08 d17

#include "distributed_grid.h"
#include "common/serialization.h"
#include "common/util/hash.h"
#include "pcl/pcl_interface_communicator.h"
#include "lib_grid/algorithms/debug_util.h"

using namespace std;

namespace ug
{
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	constructor / destructor
DistributedGridManager::
DistributedGridManager() :
	m_aElemInfoVrt("DistributedGridManager_ElemInfoVrt", false),
	m_aElemInfoEdge("DistributedGridManager_ElemInfoEdge", false),
	m_aElemInfoFace("DistributedGridManager_ElemInfoFace", false),
	m_aElemInfoVol("DistributedGridManager_ElemInfoVol", false)
{
	m_interfaceManagementEnabled = true;
	m_bOrderedInsertionMode = false;
	m_bElementDeletionMode = false;
	m_pGrid = NULL;
}

DistributedGridManager::
DistributedGridManager(MultiGrid& grid) :
	m_aElemInfoVrt("DistributedGridManager_ElemInfoVrt", false),
	m_aElemInfoEdge("DistributedGridManager_ElemInfoEdge", false),
	m_aElemInfoFace("DistributedGridManager_ElemInfoFace", false),
	m_aElemInfoVol("DistributedGridManager_ElemInfoVol", false)
{
	m_interfaceManagementEnabled = true;
	m_bOrderedInsertionMode = false;
	m_bElementDeletionMode = false;
	m_pGrid = NULL;
	assign(grid);
}

DistributedGridManager::		
~DistributedGridManager()
{
	if(m_pGrid){
		m_pGrid->unregister_observer(this);
	//	remove attached data
		m_pGrid->detach_from_vertices(m_aElemInfoVrt);
		m_pGrid->detach_from_edges(m_aElemInfoEdge);
		m_pGrid->detach_from_faces(m_aElemInfoFace);
		m_pGrid->detach_from_volumes(m_aElemInfoVol);
	}
}

////////////////////////////////////////////////////////////////////////
//	assignment
void
DistributedGridManager::
assign(MultiGrid& grid)
{
	set_grid(&grid);
}

void DistributedGridManager::
set_grid(Grid* grid)
{
	if(m_pGrid){
		m_pGrid->unregister_observer(this);
		free_grid_data();
	}
	
	if(grid){
		m_pGrid = dynamic_cast<MultiGrid*>(grid);
		m_pGrid->register_observer(this, OT_FULL_OBSERVER);
		UG_ASSERT(m_pGrid, "Only MultiGrids are supported by the DistributedGridManager"
						   " in its current implementation.");
						   
	//	attach element infos
		grid->attach_to_vertices(m_aElemInfoVrt);
		grid->attach_to_edges(m_aElemInfoEdge);
		grid->attach_to_faces(m_aElemInfoFace);
		grid->attach_to_volumes(m_aElemInfoVol);
		
	//	access them
		m_aaElemInfoVRT.access(*grid, m_aElemInfoVrt);
		m_aaElemInfoEDGE.access(*grid, m_aElemInfoEdge);
		m_aaElemInfoFACE.access(*grid, m_aElemInfoFace);
		m_aaElemInfoVOL.access(*grid, m_aElemInfoVol);
		
	//	initialise the element statuses. This is automatically done
	//	on a call to grid_layout_changed.
		grid_layouts_changed();						   
	}
}

void DistributedGridManager::free_grid_data()
{
//	remove attached data
	if(m_pGrid){
		m_pGrid->detach_from_vertices(m_aElemInfoVrt);
		m_pGrid->detach_from_edges(m_aElemInfoEdge);
		m_pGrid->detach_from_faces(m_aElemInfoFace);
		m_pGrid->detach_from_volumes(m_aElemInfoVol);
		m_pGrid = NULL;
	}

//	clear the layout-map
	m_gridLayoutMap = GridLayoutMap();
}

void DistributedGridManager::
grid_to_be_destroyed(Grid* grid)
{
	if(m_pGrid)
		free_grid_data();
}
/*
template <class TElem>
void DistributedGridManager::set_preliminary_ghost_states()
{//	Works but is currently unused and would only be required, if
//	horizontal interfaces between vertical-masters would exist.

	typedef typename geometry_traits<TElem>::iterator iterator;

	for(iterator iter = m_pGrid->begin<TElem>(); iter != m_pGrid->end<TElem>(); ++iter)
	{
		TElem* elem = *iter;
		byte status = get_status(elem);

		bool isGhost = true;

		if((status & ES_V_MASTER)){
		//	the element is a vertical master and thus potentially is a ghost.
		//	however - if an horizontal interface exists, which connects the
		//	vertical master with its vertical slaves, then it is not regarded
		//	as a ghost.
			if(status & (ES_H_MASTER | ES_H_SLAVE))
			{
			//	it lies in vertical interfaces, too. check if the element
			//	lies in a horizontal and in a vertical interface, which point
			//	to the same proc.
				const ElementInfo<TElem>& inf = elem_info(elem);

				typedef typename ElementInfo<TElem>::ConstEntryIterator EntryIter;

				for(EntryIter iter_v = inf.entries_begin();
					iter_v != inf.entries_end(); ++iter_v)
				{
					int intfcTypeV = inf.get_interface_type(iter_v);
					if(intfcTypeV == ES_V_MASTER){
					//	search for a H_MASTER or H_SLAVE Entry with the same target proc
						int vtarget = inf.get_target_proc(iter_v);

						for(EntryIter iter_h = inf.entries_begin();
							iter_h != inf.entries_end(); ++iter_h)
						{
							int intfcTypeH = inf.get_interface_type(iter_h);
							if(intfcTypeH & (ES_H_MASTER | ES_H_SLAVE)){
							//	compare target procs
								if(inf.get_target_proc(iter_h) == vtarget){
								//	we've got a connection.
								//	The node thus isn't a ghost.
									isGhost = false;
									break;
								}
							}
						}
					}
				}
			}
		}
		else
			isGhost = false;

		if(isGhost)
			elem_info(elem).set_status(status | ES_GHOST);
		else
			elem_info(elem).set_status(status & (~ES_GHOST));
	}
}
*/
/*
////////////////////////////////////////////////////////////////////////
void DistributedGridManager::update_ghost_states()
{//	Works but is currently unused and would only be required, if
//	horizontal interfaces between vertical-masters would exist.
//	first we'll calculate a preliminary ghost state
	set_preliminary_ghost_states<VertexBase>();
	set_preliminary_ghost_states<EdgeBase>();
	set_preliminary_ghost_states<Face>();
	set_preliminary_ghost_states<Volume>();

//	now we have to convert some ghost to non-ghosts.
//	if a ghost is a part of a higher dimensional non-ghost, then
//	it has to be marked as non-ghost, too.
	vector<Face*> faces;
	vector<EdgeBase*> edges;

	for(VolumeIterator iter = m_pGrid->begin<Volume>();
		iter != m_pGrid->end<Volume>(); ++iter)
	{
		Volume* v = *iter;
		if(contains_status(v, ES_GHOST))
			continue;

		CollectAssociated(faces, *m_pGrid, v);
		for(size_t i = 0; i < faces.size(); ++i)
			elem_info(faces[i]).set_status(get_status(faces[i]) & (~ES_GHOST));

		CollectAssociated(edges, *m_pGrid, v);
		for(size_t i = 0; i < edges.size(); ++i)
			elem_info(edges[i]).set_status(get_status(edges[i]) & (~ES_GHOST));

		for(size_t i = 0; i < v->num_vertices(); ++i)
			elem_info(v->vertex(i)).set_status(get_status(v->vertex(i)) & (~ES_GHOST));
	}

	for(FaceIterator iter = m_pGrid->begin<Face>();
		iter != m_pGrid->end<Face>(); ++iter)
	{
		Face* f = *iter;
		if(contains_status(f, ES_GHOST))
			continue;

		CollectAssociated(edges, *m_pGrid, f);
		for(size_t i = 0; i < edges.size(); ++i)
			elem_info(edges[i]).set_status(get_status(edges[i]) & (~ES_GHOST));

		for(size_t i = 0; i < f->num_vertices(); ++i)
			elem_info(f->vertex(i)).set_status(get_status(f->vertex(i)) & (~ES_GHOST));
	}

	for(EdgeBaseIterator iter = m_pGrid->begin<EdgeBase>();
		iter != m_pGrid->end<EdgeBase>(); ++iter)
	{
		EdgeBase* e = *iter;
		if(contains_status(e, ES_GHOST))
			continue;

		for(size_t i = 0; i < e->num_vertices(); ++i)
			elem_info(e->vertex(i)).set_status(get_status(e->vertex(i)) & (~ES_GHOST));
	}
}
*/

////////////////////////////////////////////////////////////////////////
void DistributedGridManager::grid_layouts_changed(bool addedElemsOnly)
{
//	I don't think that addedElemsOnly is correctly implemented in the moment.
//	I thus disabled it here.
	//if(!addedElemsOnly)
	{
	//	first we have to reset all elem infos
		reset_elem_infos<VertexBase>();
		reset_elem_infos<EdgeBase>();
		reset_elem_infos<Face>();
		reset_elem_infos<Volume>();
	}
	
//	call for each layout in the grid-layout the corresponding
//	init_elem_status_from_layout function.
//	every layout has multiple levels
	update_all_elem_infos<VertexBase>();
	update_all_elem_infos<EdgeBase>();
	update_all_elem_infos<Face>();
	update_all_elem_infos<Volume>();

//	update ghost states
	//update_ghost_states();
}

////////////////////////////////////////////////////////////////////////
template <class TGeomObj>
void DistributedGridManager::update_all_elem_infos()
{
	update_elem_info<TGeomObj>(m_gridLayoutMap, INT_H_MASTER,
							   ES_IN_INTERFACE | ES_H_MASTER);

	update_elem_info<TGeomObj>(m_gridLayoutMap, INT_H_SLAVE,
							   ES_IN_INTERFACE | ES_H_SLAVE);

	update_elem_info<TGeomObj>(m_gridLayoutMap, INT_V_MASTER,
							   ES_V_MASTER, true);

	update_elem_info<TGeomObj>(m_gridLayoutMap, INT_V_SLAVE,
							   ES_V_SLAVE, true);
}
		
////////////////////////////////////////////////////////////////////////
template <class TGeomObj>
void DistributedGridManager::reset_elem_infos()
{
	for(typename geometry_traits<TGeomObj>::iterator
		iter = m_pGrid->begin<TGeomObj>();
		iter != m_pGrid->end<TGeomObj>(); ++iter)
	{
		elem_info(*iter).reset();
	}
}

////////////////////////////////////////////////////////////////////////
template <class TGeomObj, class TLayoutMap>
void DistributedGridManager::
update_elem_info(TLayoutMap& layoutMap, int nodeType, byte newStatus, bool addStatus)
{	
	typedef typename TLayoutMap::template Types<TGeomObj>::Layout Layout;
	if(layoutMap.template has_layout<TGeomObj>(nodeType)){
	//	get the layout
		 Layout& layout = layoutMap.template get_layout<TGeomObj>(nodeType);
		 
	//	iterate through the levels of the layout
		for(size_t l = 0; l < layout.num_levels(); ++l){
		//	iterate through the interfaces of the layout
			for(typename Layout::iterator iiter = layout.begin(l);
				iiter != layout.end(l); ++iiter)
			{
				typename Layout::Interface& interface = layout.interface(iiter);
				//int procID = layout.proc_id(iiter);
				
			///	iterate through the elements of the interface
				for(typename Layout::Interface::iterator iter = interface.begin();
					iter != interface.end(); ++iter)
				{
				//	set the new status
					//set_status(*iter, newStatus);
					if(addStatus)
						elem_info(interface.get_element(iter)).set_status(
							elem_info(interface.get_element(iter)).get_status() | newStatus);
					else
						elem_info(interface.get_element(iter)).set_status(newStatus);
					
				//	add the iterator to the iterator list and set the proc-id
					//if(newStatus & ES_IN_INTERFACE)
					elem_info(interface.get_element(iter)).add_entry(&interface, iter, nodeType);
				/*
					elem_info(*iter).lstProcIterPairs.push_back(
								typename ElementInfo<TGeomObj>(procID, iter));
				*/
				}
			}
		}
	}
}

byte DistributedGridManager::
get_status(GridObject* go) const
{
	int baseType = go->base_object_id();
	switch(baseType)
	{
		case VERTEX:
			return get_status(static_cast<VertexBase*>(go));
		case EDGE:
			return get_status(static_cast<EdgeBase*>(go));
		case FACE:
			return get_status(static_cast<Face*>(go));
		case VOLUME:
			return get_status(static_cast<Volume*>(go));
	}
	return 0;
}

////////////////////////////////////////////////////////////////////////
//	begin_ / end_element_creation
void
DistributedGridManager::
begin_ordered_element_insertion()
{
	assert(!m_bElementDeletionMode);
	m_bOrderedInsertionMode = true;
}
/*
template <class TElem, class TCommGrp>
void DistributedGridManager::
add_element_to_interface(TElem* pElem, TCommGrp& commGrp,
						int procID, InterfaceNodeTypes nodeType)
{
	typename TCommGrp::HNODE hNode;
	
//	if ele hasn't already been inserted into the group, then do it now.
	if(!check_status(pElem, ES_IN_COMM_GRP))
	{
		hNode = commGrp.add_node(pElem);
		commGrp.set_node_type(hNode, nodeType);
		set_status(pElem, ES_IN_COMM_GRP);
	}
	else
		hNode = commGrp.get_handle(pElem);
	
//	add ele to the interface
	commGrp.add_node_to_interface(hNode, procID);
	
//	update status
	if(nodeType == pcl::NT_MASTER)
		set_status(pElem, get_status(pElem) | ES_H_MASTER);
	else if(nodeType == pcl::NT_SLAVE)
		set_status(pElem, get_status(pElem) | ES_H_SLAVE);
		
	set_status(pElem, get_status(pElem) | ES_IN_INTERFACE);

}
*/
template <class TElem>
void DistributedGridManager::
add_element_to_interface(TElem* pElem, int procID)
{
	byte status = get_status(pElem);
	
	typename GridLayoutMap::Types<TElem>::Interface::iterator iter;
	typename GridLayoutMap::Types<TElem>::Interface* interface;
	int intfcType = ES_NONE;
	byte s = get_status(pElem);
	
	if(status & ES_H_MASTER){
		interface = &m_gridLayoutMap.get_layout<TElem>(INT_H_MASTER)
						.interface(procID, m_pGrid->get_level(pElem));
		iter = interface->push_back(pElem);
		s &= (~ES_SCHEDULED_FOR_INTERFACE);
		s |= ES_IN_INTERFACE | ES_H_MASTER;
		elem_info(pElem).set_status(s);
		intfcType = ES_H_MASTER;
	}
	else{
		UG_ASSERT(status & ES_H_SLAVE, "interface-elements have to be either master or slave!");
		interface = &m_gridLayoutMap.get_layout<TElem>(INT_H_SLAVE)
						.interface(procID, m_pGrid->get_level(pElem));
		iter = interface->push_back(pElem);

		s &= (~ES_SCHEDULED_FOR_INTERFACE);
		s |= ES_IN_INTERFACE | ES_H_SLAVE;
		elem_info(pElem).set_status(s);
		intfcType = ES_H_SLAVE;
	}
	
//	add the interface-entry to the info
	elem_info(pElem).add_entry(interface, iter, intfcType);
}
		
template <class TScheduledElemMap>
void DistributedGridManager::
perform_ordered_element_insertion(TScheduledElemMap& elemMap)
{		
	for(typename TScheduledElemMap::iterator iter = elemMap.begin();
		iter != elemMap.end(); ++iter)
	{
		ScheduledElement& schedElem = iter->second;
		int objType = schedElem.geomObj->base_object_id();
		switch(objType)
		{
			case VERTEX:
				add_element_to_interface(static_cast<VertexBase*>(schedElem.geomObj),
										schedElem.connectedProcID);
				break;
			case EDGE:
				add_element_to_interface(static_cast<EdgeBase*>(schedElem.geomObj),
										schedElem.connectedProcID);
				break;
			case FACE:
				add_element_to_interface(static_cast<Face*>(schedElem.geomObj),
										schedElem.connectedProcID);
				break;
			case VOLUME:
				add_element_to_interface(static_cast<Volume*>(schedElem.geomObj),
										schedElem.connectedProcID);
				break;
		}
	}
}
										
void
DistributedGridManager::
end_ordered_element_insertion()
{
	if(m_bOrderedInsertionMode && m_interfaceManagementEnabled){
	//TODO: support all elements
		perform_ordered_element_insertion(m_vrtMap);
		perform_ordered_element_insertion(m_edgeMap);
		perform_ordered_element_insertion(m_faceMap);
		perform_ordered_element_insertion(m_volMap);

		clear_scheduled_elements();
	}

	m_bOrderedInsertionMode = false;

//	we have to update ghost-states
	//update_ghost_states();
}

void
DistributedGridManager::
clear_scheduled_elements()
{
	m_vrtMap.clear();
	m_edgeMap.clear();
	m_faceMap.clear();
	m_volMap.clear();
}

//	adds a ScheduledElement to scheduled-element-maps, depending on
//	the interfaces in which pParent is contained.
template <class TElem, class TScheduledElemMap, class TParent>
void DistributedGridManager::
schedule_element_for_insertion(TScheduledElemMap& elemMap,
										TElem* elem,
										TParent* pParent)
{
	typedef typename ElementInfo<TParent>::EntryIterator entry_iter;
	ElementInfo<TParent>& parentInfo = elem_info(pParent);

//	schedule one element for each horizontal parent-interface
	if(parentInfo.is_interface_element()){
		for(entry_iter iter = parentInfo.entries_begin();
			iter != parentInfo.entries_end(); ++iter)
		{
			int intfcType = parentInfo.get_interface_type(iter);
			if(!(intfcType & (ES_V_MASTER | ES_V_SLAVE))){
				//UG_DLOG(LIB_GRID, 3, parentInfo.get_target_proc(iter) << ", ");
				elemMap.insert(make_pair(parentInfo.get_local_id(iter),
						ScheduledElement(elem, parentInfo.get_target_proc(iter))));
			}
		}
	
	//	set the status
		if(parentInfo.get_status() & (ES_H_MASTER))
			elem_info(elem).set_status(ES_H_MASTER | ES_SCHEDULED_FOR_INTERFACE);
		else{
			UG_ASSERT(parentInfo.get_status() & (ES_H_SLAVE), "interface-elements have to be either master or slave!");
			elem_info(elem).set_status(ES_H_SLAVE | ES_SCHEDULED_FOR_INTERFACE);
		}
	}
}

//	new elements are handled here.
template <class TElem>
void DistributedGridManager::
handle_created_element(TElem* pElem, GridObject* pParent,
						bool replacesParent)
{
	if(replacesParent){
	//	we have to replace the parent entry. To do this
	//	we'll replace all occurences of the parent in all interfaces.
	//	If a replace takes place, the parent has to be of the same type
	//	as pElem.
		TElem* parent = dynamic_cast<TElem*>(pParent);
		if(parent){
		//	copy info
			ElementInfo<TElem>& elemInfo = elem_info(pElem);
			elemInfo = elem_info(parent);

			if(elemInfo.is_interface_element()){
			//	update all pointers
				for(typename ElementInfo<TElem>::EntryIterator iter = elemInfo.entries_begin();
					iter != elemInfo.entries_end(); ++iter)
				{
					typename ElementInfo<TElem>::Entry& entry = *iter;
					entry.m_interface->get_element(entry.m_interfaceElemIter) = pElem;
				}
			}

		//	clear the parent-info.
			elem_info(parent).reset();

		//	we also have to check wheter the new element has been converted
		//	from normal to constrained while lying in a v-interface during coarsening
			if(m_bElementDeletionMode && pElem->is_constrained()
				&& (!pParent->is_constrained()) && is_in_vertical_interface(pElem))
			{
				got_new_constrained_vertical(pElem);
			}

			return;
		}
		else{
			throw(UGError("Can't replace an element with an element of another type."));
		}
	}

	elem_info(pElem).set_status(ES_NONE);
	
	if(!m_interfaceManagementEnabled)
		return;
		
//	if there is no parent, we can immediately leave.
	if(!pParent)
		return;
	
//	if the parent is not in an interface or scheduled to for an interface,
//	there is nothing to do either
	if(!(get_status(pParent) &
		(ES_IN_INTERFACE | ES_SCHEDULED_FOR_INTERFACE)))
		return;

	int parentType = pParent->base_object_id();

//	if ordered insertion mode is active, we have to insert the elements
//	into a map, instead of directly adding them to the group.
	if(m_bOrderedInsertionMode)
	{
		switch(parentType)
		{
			case VERTEX:
				//UG_DLOG(LIB_GRID, 3, "scheduling element with vertex-parent to interfaces ");
				schedule_element_for_insertion(m_vrtMap,
												pElem,
												(VertexBase*)pParent);
				//UG_DLOG(LIB_GRID, 3, endl);
				break;
				
			case EDGE:
				//UG_DLOG(LIB_GRID, 3, "scheduling element with edge-parent to interfaces ");
				schedule_element_for_insertion(m_edgeMap, pElem,
												(EdgeBase*)pParent);
				//UG_DLOG(LIB_GRID, 3, endl);
				break;

			case FACE:
				//UG_DLOG(LIB_GRID, 3, "scheduling element with face-parent to interfaces ");
				schedule_element_for_insertion(m_faceMap, pElem,
												(Face*)pParent);
				//UG_DLOG(LIB_GRID, 3, endl);
				break;

			case VOLUME:
				//UG_DLOG(LIB_GRID, 3, "scheduling element with volume-parent to interfaces ");
				schedule_element_for_insertion(m_volMap, pElem,
												(Volume*)pParent);
				//UG_DLOG(LIB_GRID, 3, endl);
				break;
		}
	}
	else
	{
	//	directly add the edge to the communication group
		assert(!"only sorted mode supported in the moment!");
	}
}

void DistributedGridManager::
vertex_created(Grid* grid, VertexBase* vrt, GridObject* pParent,
				bool replacesParent)
{
	handle_created_element(vrt, pParent, replacesParent);
}


void DistributedGridManager::
edge_created(Grid* grid, EdgeBase* e, GridObject* pParent,
			 bool replacesParent)
{
	handle_created_element(e, pParent, replacesParent);
}

void DistributedGridManager::
face_created(Grid* grid, Face* f, GridObject* pParent,
			 bool replacesParent)
{
	handle_created_element(f, pParent, replacesParent);
}

void DistributedGridManager::
volume_created(Grid* grid, Volume* v, GridObject* pParent,
			   bool replacesParent)
{
	handle_created_element(v, pParent, replacesParent);
}


////////////////////////////////////////////////////////////////////////////////
//	Element deletion
void DistributedGridManager::
begin_element_deletion()
{
	assert(!m_bOrderedInsertionMode);
	m_bElementDeletionMode = true;
	m_newConstrainedVerticalVrts.clear();
	m_newConstrainedVerticalEdges.clear();
	m_newConstrainedVerticalFaces.clear();
}


template <class TLayout>
class ComPol_NewConstrainedVerticals : public pcl::ICommunicationPolicy<TLayout>
{
	public:
		typedef TLayout								Layout;
		typedef typename Layout::Type				GeomObj;
		typedef typename Layout::Element			Element;
		typedef typename Layout::Interface			Interface;
		typedef typename Interface::const_iterator	InterfaceIter;

	/**	Note that a reference to newConstrained is passed to the constructor
	 * and that its content may be changed by some methods in this class.*/
		ComPol_NewConstrainedVerticals(DistributedGridManager* dgm,
									   std::vector<GeomObj*>& newConstrained) :
			m_newConstrained(newConstrained),
			m_dgm(dgm),
			m_hash(newConstrained.size() + 97),
			m_localHMasterCount(0),
			m_exchangeVMasterRanks(false),
			m_initialHandshake(false)
			//m_checkHOrder(false)
		{
			m_hash.reserve(size_t(newConstrained.size() * 1.5));

		//	insert each new constrained into the hash.
			for(size_t i_nc = 0; i_nc < newConstrained.size(); ++i_nc){
				GeomObj* e = newConstrained[i_nc];
				create_initial_hash_entry(e);
			}
		}

		void create_initial_hash_entry(GeomObj* e)
		{
			if(m_dgm->contains_status(e, ES_H_MASTER)){
				m_hash.insert(e, Entry(pcl::ProcRank(), m_localHMasterCount));
				++m_localHMasterCount;
			}
			else if(m_dgm->contains_status(e, ES_V_MASTER)
					&& (!m_dgm->contains_status(e, ES_H_SLAVE)))
			{
			//	find lowest connected vslave proc
				vector<pair<int, size_t> >	interfaceEntries;
				m_dgm->collect_interface_entries(interfaceEntries, e, ES_V_MASTER);
				UG_ASSERT(!interfaceEntries.empty(),
						  "Elem with type " << e->base_object_id() << " at " <<
						   GetGridObjectCenter(*m_dgm->get_assigned_grid(), e)
						   << " is marked as v-master but is not contained in a vslave interface!");
				int lp = interfaceEntries.front().first;
				for(size_t i = 1; i < interfaceEntries.size(); ++i)
					lp = min(lp, interfaceEntries[i].first);

				m_hash.insert(e, Entry(lp, -1));
			}
			else
				m_hash.insert(e, Entry(-1, -1));
		}

		virtual ~ComPol_NewConstrainedVerticals()	{}

		virtual int
		get_required_buffer_size(const Interface& interface)
		{return -1;}

		virtual bool
		collect(ug::BinaryBuffer& buff, const Interface& interface)
		{
			vector<pair<int, size_t> > vInterfaces;

			int targetProc = interface.get_target_proc();
			int counter = 0;
			for(InterfaceIter iter = interface.begin();
				iter != interface.end(); ++iter, ++counter)
			{
				bool sendVMasterRanks = false;
				Element elem = interface.get_element(iter);
				if(m_hash.has_entry(elem)){
					Serialize(buff, counter);
					Entry& entry = m_hash.get_entry(elem);
					Serialize(buff, entry.hmasterProcInfo);
					sendVMasterRanks = (entry.hmasterProcInfo.first == targetProc);
				}

				if(m_exchangeVMasterRanks){
					UG_ASSERT(!m_initialHandshake, "initial handshake and vmaster "
								 "exchange may not be active at the same time.");
					UG_ASSERT(m_dgm->contains_status(elem, ES_V_SLAVE),
							  "Only v-slaves can communicate associated v-masters.");

					if(sendVMasterRanks){
						m_dgm->collect_interface_entries(vInterfaces, elem, ES_V_SLAVE);
						int numOtherMasters = (int)vInterfaces.size() - 1;
						Serialize(buff, numOtherMasters);
						for(size_t i = 0; i < vInterfaces.size(); ++i){
							if(vInterfaces[i].first != targetProc){
								Serialize(buff, vInterfaces[i].first);
							}
						}
					}
				}
			}

			int val = -1;
			Serialize(buff, val);

			return true;
		}

		virtual bool
		extract(ug::BinaryBuffer& buff, const Interface& interface)
		{
			int index;
			Deserialize(buff, index);

			int localProc = pcl::ProcRank();
			int counter = 0;
			InterfaceIter iter = interface.begin();

			while((index != -1) && (iter != interface.end())){
				if(counter == index){
					Element elem = interface.get_element(iter);
					std::pair<int, int> val;
					Deserialize(buff, val);

//					if(m_checkHOrder){
//						if(!m_hash.has_entry(elem)){
//							UG_THROW("No matching entry on proc " << pcl::ProcRank()
//									 << " for element " << ElementDebugInfo(*m_dgm->get_assigned_grid(), elem));
//							Entry& entry = m_hash.get_entry(elem);
//							if((entry.hmasterProcInfo.first == val.first)
//								&& (entry.hmasterProcInfo.second != val.second)){
//								UG_THROW("H-Order mismatch on proc " << pcl::ProcRank()
//										<< " for element " << ElementDebugInfo(*m_dgm->get_assigned_grid(), elem));
//							}
//						}
//					}
//					else if(m_initialHandshake){
					if(m_initialHandshake){
						UG_ASSERT(!m_exchangeVMasterRanks, "initial handshake and vmaster "
								 "exchange may not be active at the same time.");
						if(!m_hash.has_entry(elem)){
							create_initial_hash_entry(elem);
							m_newConstrained.push_back(elem);
						}
					}
					else{
						UG_ASSERT(m_hash.has_entry(elem),
								  "A matching element has to exist in the local procs "
								  "new constrained list:"
								  << ElementDebugInfo(*m_dgm->get_assigned_grid(), elem));

					//	the entry whose second value is specified has the highest priority
						Entry& entry = m_hash.get_entry(elem);
						if((entry.hmasterProcInfo.second == -1) && (val.first != -1))
							entry.hmasterProcInfo = val;

						if(m_exchangeVMasterRanks){
							UG_ASSERT(m_dgm->contains_status(elem, ES_V_MASTER),
									  "Only v-masters can receive from associated v-slaves.");
							if(entry.hmasterProcInfo.first == localProc){
								int numOtherMasters;
								Deserialize(buff, numOtherMasters);
								entry.otherVMasterRanks.clear();
								entry.otherVMasterRanks.reserve(numOtherMasters);
								for(int i = 0; i < numOtherMasters; ++i){
									int om;
									Deserialize(buff, om);
									entry.otherVMasterRanks.push_back(om);
									UG_ASSERT(om != localProc, "Only other procs should arrive here!");
								}
							}
						}
					}
					Deserialize(buff, index);
				}
				++counter;
				++iter;
			}

			UG_ASSERT(index == -1, "Not all entries in the stream have been read!");
			return (index == -1);
		}

		void exchange_data()
		{
			pcl::InterfaceCommunicator<TLayout> com;
			GridLayoutMap& glm = m_dgm->grid_layout_map();

		//	we perform an initial handshake between vmasters and vslaves, to make
		//	sure, that both have matching elements in their hashes.
			UG_DLOG(LIB_GRID, 3, "  Initial handshake\n");
			m_initialHandshake = true;
			com.exchange_data(glm, INT_V_MASTER, INT_V_SLAVE, *this);
			com.exchange_data(glm, INT_V_SLAVE, INT_V_MASTER, *this);
			com.communicate();
		//	This extra communication step is required for the rare cases, in which
		//	only a vslave was registered as new-constrained. It most likely wouldn't
		//	be necessary if one would guarantee, that lower dimensional elements
		//	of vmasters are always vmasters again.
			com.exchange_data(glm, INT_V_MASTER, INT_V_SLAVE, *this);
			com.communicate();
			m_initialHandshake = false;

		//	In the setup phase the entries of all h-master nodes already received
		//	their correct values (h-masters won't change).
		//	If a v-master is not contained in a h-interface, it will assume that the
		//	new h-master will be created at the lowest v-slave. This guess may however
		//	be corrected in later communication steps.
		//	h-slaves will only receive an invalid entry, since they don't have to
		//	be adjusted at all and since they don't really contribute to the algorithm.

		//	in the first communication step we notify v-slaves about what their
		//	associated v-masters presume to be the new h-master. Most of the time
		//	this is, however, only a guess by the v-masters -- only if they are a
		//	h-master them selves, this guess should be right.
			UG_DLOG(LIB_GRID, 3, "  communicating vmasters->vslaves...\n");
			m_exchangeVMasterRanks = false;
//			if(glm.has_layout<GeomObj>(INT_V_MASTER))
//				com.send_data(glm.get_layout<GeomObj>(INT_V_MASTER), *this);
//			if(glm.has_layout<GeomObj>(INT_V_SLAVE))
//				com.receive_data(glm.get_layout<GeomObj>(INT_V_SLAVE), *this);
			com.exchange_data(glm, INT_V_MASTER, INT_V_SLAVE, *this);
			com.communicate();

		//	v-slaves now know where their associated v-masters assume that the
		//	associated h-master node lies. By checking (hmasterProcInfo.second == -1),
		//	one knows that this was only a guess.
		//	iterate over all entries. Those who contain the local proc as h-master
		//	proc have to be supplied with a local h-master index (if none is present yet)
			int localProc = pcl::ProcRank();
			for(size_t i = 0; i < m_newConstrained.size(); ++i){
				GeomObj* elem = m_newConstrained[i];
				if(!m_dgm->contains_status(elem, ES_H_SLAVE)){
					Entry& entry = m_hash.get_entry(elem);
					if((entry.hmasterProcInfo.first == localProc)
						&& (entry.hmasterProcInfo.second == -1))
					{
						entry.hmasterProcInfo.second = m_localHMasterCount++;
					}
				}
			}


			UG_DLOG(LIB_GRID, 3, "  communicating vslaves->vmasters...\n");
			m_exchangeVMasterRanks = true;
//			if(glm.has_layout<GeomObj>(INT_V_SLAVE))
//				com.send_data(glm.get_layout<GeomObj>(INT_V_SLAVE), *this);
//			if(glm.has_layout<GeomObj>(INT_V_MASTER))
//				com.receive_data(glm.get_layout<GeomObj>(INT_V_MASTER), *this);
			com.exchange_data(glm, INT_V_SLAVE, INT_V_MASTER, *this);
			com.communicate();
			m_exchangeVMasterRanks = false;

			//check_corresponding_h_order();
		}

	/**	returns a std::pair<int, int> where 'first' represents the h-master rank and
	 * where 'second' defines the order in which entries have to be added to
	 * the h-interface.*/
		std::pair<int, int> get_h_master_info(GeomObj* o)	{return m_hash.get_entry(o).hmasterProcInfo;}

	/**	returns the array of other v-master procs which build an interface to the local proc,
	 * if the local proc will be the new h-master proc and if it is also a v-master proc.*/
		std::vector<int>& other_v_masters(GeomObj* o)		{return m_hash.get_entry(o).otherVMasterRanks;}

	private:
//		void check_corresponding_h_order()
//		{
//			pcl::InterfaceCommunicator<TLayout> com;
//			GridLayoutMap& glm = m_dgm->grid_layout_map();
//
//		//	we perform an initial handshake between vmasters and vslaves, to make
//		//	sure, that both have matching elements in their hashes.
////			UG_LOG("\n DEBUG CHECK check_corresponding_h_order\n");
////			for(size_t i = 0; i < m_newConstrained.size(); ++i){
////				GeomObj* elem = m_newConstrained[i];
////				if(vector3(0.015625, -0.25, 0) == GetGridObjectCenter(*m_dgm->get_assigned_grid(), elem)){
////					UG_LOG("(0.015625, -0.25, 0) h-master: " << m_hash.get_entry(elem).hmasterProcInfo.first << "\n");
////					UG_LOG("(0.015625, -0.25, 0) h-order: " << m_hash.get_entry(elem).hmasterProcInfo.second << "\n");
////				}
////				if(vector3(0, -0.265625, 0) == GetGridObjectCenter(*m_dgm->get_assigned_grid(), elem)){
////					UG_LOG("\n");
////					UG_LOG("(0, -0.265625, 0) h-master: " << m_hash.get_entry(elem).hmasterProcInfo.first << "\n");
////					UG_LOG("(0, -0.265625, 0) h-order: " << m_hash.get_entry(elem).hmasterProcInfo.second << "\n");
////				}
////				if(vector3(0, -0.25, 0) == GetGridObjectCenter(*m_dgm->get_assigned_grid(), elem)){
////					UG_LOG("\n");
////					UG_LOG("(0, -0.25, 0) h-master: " << m_hash.get_entry(elem).hmasterProcInfo.first << "\n");
////					UG_LOG("(0, -0.25, 0) h-order: " << m_hash.get_entry(elem).hmasterProcInfo.second << "\n");
////				}
////				if(vector3(0.03125, -0.25, 0) == GetGridObjectCenter(*m_dgm->get_assigned_grid(), elem)){
////					UG_LOG("\n");
////					UG_LOG("(0.03125, -0.25, 0) h-master: " << m_hash.get_entry(elem).hmasterProcInfo.first << "\n");
////					UG_LOG("(0.03125, -0.25, 0) h-order: " << m_hash.get_entry(elem).hmasterProcInfo.second << "\n");
////				}
////				if(vector3(0.15625, -0.28125, 0) == GetGridObjectCenter(*m_dgm->get_assigned_grid(), elem)){
////					UG_LOG("\n");
////					UG_LOG("(0.15625, -0.28125, 0) h-master: " << m_hash.get_entry(elem).hmasterProcInfo.first << "\n");
////					UG_LOG("(0.15625, -0.28125, 0) h-order: " << m_hash.get_entry(elem).hmasterProcInfo.second << "\n");
////				}
////			}
//
//			m_checkHOrder = true;
//			com.exchange_data(glm, INT_V_MASTER, INT_V_SLAVE, *this);
//			com.exchange_data(glm, INT_V_SLAVE, INT_V_MASTER, *this);
//			com.communicate();
//			m_checkHOrder = false;
//		}

		struct Entry{
			Entry(int hmasterRank, int localHMasterCount) :
				hmasterProcInfo(hmasterRank, localHMasterCount)	{}

			std::pair<int, int>	hmasterProcInfo;
			std::vector<int> otherVMasterRanks;
		};

		std::vector<GeomObj*>&	m_newConstrained;
		DistributedGridManager* m_dgm;
		Hash<GeomObj*, Entry>	m_hash;
		int						m_localHMasterCount;
		bool					m_exchangeVMasterRanks;
		bool					m_initialHandshake;
		//bool					m_checkHOrder;
};

template <class TElem>
void DistributedGridManager::
create_missing_constrained_h_interfaces(vector<TElem*>& newConstrainedElems)
{
//	some notes:
//	The process on which the new hmaster lies is carefully chosen so that if an
//	element already lies in a h-interface, its hmaster entry will also be the hmaster
//	of the newly created h-interfaces.
//	h-interfaces between all v-slaves exist at this time. h-interfaces between
//	v-slaves and v-masters may exist.
//	h-interfaces are created between the new hmaster and all other associated vmasters.
//	We don't have to create h-interfaces between vslaves since those would already
//	exist if they were required.
//	Before creating an interface we make sure that no hinterface between those procs
//	exists.

	MultiGrid& mg = *m_pGrid;
	typedef typename GridLayoutMap::Types<TElem>::Layout	layout_t;

	ComPol_NewConstrainedVerticals<layout_t> compolHMasters(mg.distributed_grid_manager(),
															newConstrainedElems);

	compolHMasters.exchange_data();

	int localRank = pcl::ProcRank();

	ScheduledElemMap	scheduledElems;
	vector<pair<int, size_t> > vInterfaces;
	vector<pair<int, size_t> > hInterfaces;

	for(size_t i_nce = 0; i_nce < newConstrainedElems.size(); ++i_nce){
		TElem* e = newConstrainedElems[i_nce];
	//	nothing to do for h-slave entries, since the h-master won't change.
		if(contains_status(e, ES_H_SLAVE))
			continue;

		std::pair<int, int> hmasterInfo = compolHMasters.get_h_master_info(e);
		int hmasterRank = hmasterInfo.first;
		int hmasterOrder = hmasterInfo.second;

		UG_ASSERT(hmasterRank != -1, "A hmasterRank has to be provided!");
		UG_ASSERT(hmasterOrder != -1, "HMaster orders have not been communicated properly");

		if(hmasterRank == localRank){
		//	make the local element a hmaster and create hinterfaces to all other vmasters
			if(contains_status(e, ES_V_SLAVE)){
				collect_interface_entries(hInterfaces, e, ES_H_SLAVE);
				collect_interface_entries(hInterfaces, e, ES_H_MASTER, false);
				collect_interface_entries(vInterfaces, e, ES_V_SLAVE);

				for(size_t i_v = 0; i_v < vInterfaces.size(); ++i_v){
					int tp = vInterfaces[i_v].first;
				//	make sure that no h-interface to this process exists already!
					bool hInterfaceExists = false;
					for(size_t i_h = 0; i_h < hInterfaces.size(); ++i_h){
						if(tp == hInterfaces[i_h].first){
							hInterfaceExists = true;
							break;
						}
					}

					if(!hInterfaceExists){
						scheduledElems.insert(make_pair(hmasterOrder, ScheduledElement(e, tp)));
						elem_info(e).set_status(get_status(e) | ES_H_MASTER | ES_SCHEDULED_FOR_INTERFACE);
					}
				}
			}
			else{
				UG_ASSERT(contains_status(e, ES_V_MASTER), "Only vslaves and vmasters should be handled here!");
				collect_interface_entries(hInterfaces, e, ES_H_SLAVE);
				collect_interface_entries(hInterfaces, e, ES_H_MASTER, false);
				vector<int>& otherVMasters = compolHMasters.other_v_masters(e);
				for(size_t i_om = 0; i_om < otherVMasters.size(); ++i_om){
					int om = otherVMasters[i_om];
				//	check whether a h-interface to that process exists already
					bool hInterfaceExists = false;
					for(size_t i_h = 0; i_h < hInterfaces.size(); ++i_h){
						if(om == hInterfaces[i_h].first){
							hInterfaceExists = true;
							break;
						}
					}

					if(!hInterfaceExists){
						scheduledElems.insert(make_pair(hmasterOrder, ScheduledElement(e, om)));
						elem_info(e).set_status(get_status(e) | ES_H_MASTER | ES_SCHEDULED_FOR_INTERFACE);
					}
				}

				collect_interface_entries(vInterfaces, e, ES_V_MASTER);
				for(size_t i_vm = 0; i_vm < vInterfaces.size(); ++i_vm){
					int tp = vInterfaces[i_vm].first;
				//	check whether a h-interface to that process exists already
					bool hInterfaceExists = false;
					for(size_t i_h = 0; i_h < hInterfaces.size(); ++i_h){
						if(tp == hInterfaces[i_h].first){
							hInterfaceExists = true;
							break;
						}
					}

					if(!hInterfaceExists){
						scheduledElems.insert(make_pair(hmasterOrder, ScheduledElement(e, tp)));
						elem_info(e).set_status(get_status(e) | ES_H_MASTER | ES_SCHEDULED_FOR_INTERFACE);
					}
				}
			}
		}
		else if(contains_status(e, ES_V_MASTER)){
			UG_ASSERT(!contains_status(e, ES_H_MASTER), "This proc should not be considered as h-master proc!");
		//	create a hslaveinterface to the hmasterRank.
			collect_interface_entries(hInterfaces, e, ES_H_SLAVE);
			collect_interface_entries(hInterfaces, e, ES_H_MASTER, false);

		//	make sure that no h-interface to this process exists already!
			bool hInterfaceExists = false;
			for(size_t i_h = 0; i_h < hInterfaces.size(); ++i_h){
				if(hmasterRank == hInterfaces[i_h].first){
					hInterfaceExists = true;
					break;
				}
			}

			if(!hInterfaceExists){
				scheduledElems.insert(make_pair(hmasterOrder, ScheduledElement(e, hmasterRank)));
				elem_info(e).set_status(get_status(e) | ES_H_SLAVE | ES_SCHEDULED_FOR_INTERFACE);
			}
		}
	}

//	finally insert the scheduled elements in the actual interfaces
	perform_ordered_element_insertion(scheduledElems);
}

void DistributedGridManager::
end_element_deletion()
{
	UG_DLOG(LIB_GRID, 1, "DistributedGridManager start - end_element_deletion\n");
// todo	defragment interfaces
	m_bElementDeletionMode = false;

//	we have to make sure that h-interfaces exist between new constrained vertical
//	interface elements.
	MultiGrid& mg = *m_pGrid;

//	make sure that all associated elements are contained in the newConstrained... containers.
	std::vector<EdgeBase*> edges;
	for(size_t i = 0; i < m_newConstrainedVerticalFaces.size(); ++i){
		CollectAssociated(edges, mg, m_newConstrainedVerticalFaces[i], false);
	}

	mg.begin_marking();
	mg.mark(m_newConstrainedVerticalEdges.begin(), m_newConstrainedVerticalEdges.end());

	for(size_t i = 0; i < edges.size(); ++i){
		EdgeBase* e = edges[i];
		if(!mg.is_marked(e)){
			mg.mark(e);
			m_newConstrainedVerticalEdges.push_back(e);
		}
	}

	mg.mark(m_newConstrainedVerticalVrts.begin(), m_newConstrainedVerticalVrts.end());
	for(size_t i = 0; i < m_newConstrainedVerticalEdges.size(); ++i){
		EdgeBase* e = m_newConstrainedVerticalEdges[i];
		for(size_t j = 0; j < 2; ++j){
			if(!mg.is_marked(e->vertex(j))){
				mg.mark(e->vertex(j));
				m_newConstrainedVerticalVrts.push_back(e->vertex(j));
			}
		}
	}

	mg.end_marking();

	UG_DLOG(LIB_GRID, 2, "  creating missing constrained h interfaces for vertices...\n");
	create_missing_constrained_h_interfaces(m_newConstrainedVerticalVrts);
	UG_DLOG(LIB_GRID, 2, "  creating missing constrained h interfaces for edges...\n");
	create_missing_constrained_h_interfaces(m_newConstrainedVerticalEdges);
	UG_DLOG(LIB_GRID, 2, "  creating missing constrained h interfaces for faces...\n");
	create_missing_constrained_h_interfaces(m_newConstrainedVerticalFaces);

	UG_DLOG(LIB_GRID, 1, "DistributedGridManager stop - end_element_deletion\n");
}

template <class TElem>
void DistributedGridManager::
element_to_be_erased(TElem* elem)
{
	if(!m_interfaceManagementEnabled)
		return;

	ElementInfo<TElem>& elemInfo = elem_info(elem);

	if(elemInfo.is_interface_element()){
		UG_ASSERT(m_bElementDeletionMode, "Call begin_element_deletion() before deleting elements.");
	//	erase the element from all associated interfaces
		for(typename ElementInfo<TElem>::EntryIterator iter = elemInfo.entries_begin();
			iter != elemInfo.entries_end(); ++iter)
		{
			typename ElementInfo<TElem>::Entry& entry = *iter;
			entry.m_interface->erase(entry.m_interfaceElemIter);
		}
	}
}

void DistributedGridManager::
vertex_to_be_erased(Grid* grid, VertexBase* vrt, VertexBase* replacedBy)
{
	if(!replacedBy)
		element_to_be_erased(vrt);
}

void DistributedGridManager::
edge_to_be_erased(Grid* grid, EdgeBase* e, EdgeBase* replacedBy)
{
	if(!replacedBy)
		element_to_be_erased(e);
}

void DistributedGridManager::
face_to_be_erased(Grid* grid, Face* f, Face* replacedBy)
{
	if(!replacedBy)
		element_to_be_erased(f);
}

void DistributedGridManager::
volume_to_be_erased(Grid* grid, Volume* vol, Volume* replacedBy)
{
	if(!replacedBy)
		element_to_be_erased(vol);
}

}// end of namespace

