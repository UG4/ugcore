// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y09 m08 d17

#include "distributed_grid.h"

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
	//	remove attached data
		m_pGrid->detach_from_vertices(m_aElemInfoVrt);
		m_pGrid->detach_from_edges(m_aElemInfoEdge);
		m_pGrid->detach_from_faces(m_aElemInfoFace);
		m_pGrid->detach_from_volumes(m_aElemInfoVol);
		m_pGrid = NULL;
		
	//	clear the layout-map
		m_gridLayoutMap = GridLayoutMap();
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

void DistributedGridManager::
grid_to_be_destroyed(Grid* grid)
{
	if(m_pGrid)
		set_grid(NULL);
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
get_status(GeometricObject* go) const
{
	int baseType = go->base_object_id();
	switch(baseType)
	{
		case VERTEX:
			return get_status(static_cast<VertexBase*>(go));;
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
	
	if(status & ES_H_MASTER){
		interface = &m_gridLayoutMap.get_layout<TElem>(INT_H_MASTER)
						.interface(procID, m_pGrid->get_level(pElem));
		iter = interface->push_back(pElem);
		elem_info(pElem).set_status(ES_IN_INTERFACE | ES_H_MASTER);
		intfcType = ES_H_MASTER;
	}
	else{
		UG_ASSERT(status & ES_H_SLAVE, "interface-elements have to be either master or slave!");
		interface = &m_gridLayoutMap.get_layout<TElem>(INT_H_SLAVE)
						.interface(procID, m_pGrid->get_level(pElem));
		iter = interface->push_back(pElem);
						
		elem_info(pElem).set_status(ES_IN_INTERFACE | ES_H_SLAVE);
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
				UG_DLOG(LIB_GRID, 3, parentInfo.get_target_proc(iter) << ", ");
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
handle_created_element(TElem* pElem, GeometricObject* pParent,
						bool replacesParent)
{
	if(replacesParent){
	//	we only have to replace the parent entry. To do this
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
				UG_DLOG(LIB_GRID, 3, "scheduling element with vertex-parent to interfaces ");
				schedule_element_for_insertion(m_vrtMap,
												pElem,
												(VertexBase*)pParent);
				UG_DLOG(LIB_GRID, 3, endl);
				break;
				
			case EDGE:
				UG_DLOG(LIB_GRID, 3, "scheduling element with edge-parent to interfaces ");
				schedule_element_for_insertion(m_edgeMap, pElem,
												(EdgeBase*)pParent);
				UG_DLOG(LIB_GRID, 3, endl);
				break;

			case FACE:
				UG_DLOG(LIB_GRID, 3, "scheduling element with face-parent to interfaces ");
				schedule_element_for_insertion(m_faceMap, pElem,
												(Face*)pParent);
				UG_DLOG(LIB_GRID, 3, endl);
				break;

			case VOLUME:
				UG_DLOG(LIB_GRID, 3, "scheduling element with volume-parent to interfaces ");
				schedule_element_for_insertion(m_volMap, pElem,
												(Volume*)pParent);
				UG_DLOG(LIB_GRID, 3, endl);
				break;
		}
	}
	else
	{
	//	directly add the edge to the communication group
		assert(!"only sorted mode supported in the moment!");
	}
}
/*
// when handle_erased_element is implemented, one has to be careful with
// the local interface ids. Those won't be continuous after an erasure
// was performed. This is normally no problem. But when it comes to
// grid-redistribution, a method is required that brings those id-s
// back to continuous. Keep that in mind and have a look at the
// todos in RedistributeGrid.

template <class TElem, class TCommGrp>
void DistributedGrid::
handle_erased_element(TElem* e, TCommGrp& commGrp)
{
//	for each deleted element a check has to be performed whether it was
//	contained in the associated communication-group.
//	If so it has to be removed from the group.
//TODO:
//	If we are in creation mode, we also have to check whether the element
//	is contained in any map.
//	if so we'll store it in a hash and will check later, on insertion of new
//	edges into the interface, whether the pointer is contained in the hash.
//	If so we won't add it.
	if(check_status(e, ES_IN_COMM_GRP))
	{
	//	delete the element from the communication group
		typename TCommGrp::HNODE hNode = commGrp.get_handle(e);
		commGrp.erase_node(hNode);
	}
	
	if(check_status(e, ES_SCHEDULED_FOR_COMM_GRP))
	{
	//	todo: remove the element from the map.
		assert(0);
	}
}

template <class TElem, class TCommGrp>
void DistributedGrid::
handle_replaced_element(TElem* eOld, TElem* eNew, TCommGrp& commGrp)
{
//	if the replaced element is a memeber of the elements communication group,
//	we have to perform replace-node on the group.
	if(check_status(eOld, ES_IN_COMM_GRP))
	{
	//	replace eOld wiht eNew in the communication group
		commGrp.replace_node(eOld, eNew);
	}
	
	if(check_status(eOld, ES_SCHEDULED_FOR_COMM_GRP))
	{
	//	todo: remove the element from the map.
		assert(0);
	}
}
*/
void DistributedGridManager::
vertex_created(Grid* grid, VertexBase* vrt, GeometricObject* pParent,
				bool replacesParent)
{
	handle_created_element(vrt, pParent, replacesParent);
}


void DistributedGridManager::
edge_created(Grid* grid, EdgeBase* e, GeometricObject* pParent,
			 bool replacesParent)
{
	handle_created_element(e, pParent, replacesParent);
}

void DistributedGridManager::
face_created(Grid* grid, Face* f, GeometricObject* pParent,
			 bool replacesParent)
{
	handle_created_element(f, pParent, replacesParent);
}

void DistributedGridManager::
volume_created(Grid* grid, Volume* v, GeometricObject* pParent,
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
}

void DistributedGridManager::
end_element_deletion()
{
// todo	defragment interfaces
	m_bElementDeletionMode = false;
}

template <class TElem>
void DistributedGridManager::
element_to_be_erased(TElem* elem)
{
	if(!m_interfaceManagementEnabled)
		return;

	ElementInfo<TElem>& elemInfo = elem_info(elem);

	if(elemInfo.is_interface_element()){
		assert(m_bElementDeletionMode);
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

