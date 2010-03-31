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
DistributedGridManager()
{
	m_bOrderedInsertionMode = false;
	m_pGrid = NULL;
}

DistributedGridManager::
DistributedGridManager(MultiGrid& grid)
{
	m_bOrderedInsertionMode = false;
	m_pGrid = NULL;
	assign(grid);
}

DistributedGridManager::		
~DistributedGridManager()
{
	if(m_pGrid)
		m_pGrid->unregister_observer(this);
}

////////////////////////////////////////////////////////////////////////
//	assignment
void
DistributedGridManager::
assign(MultiGrid& grid)
{
	grid.register_observer(this);
}

void DistributedGridManager::registered_at_grid(Grid* grid)
{
	if(m_pGrid)
		m_pGrid->unregister_observer(this);

	m_pGrid = dynamic_cast<MultiGrid*>(grid);
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

void DistributedGridManager::unregistered_from_grid(Grid* grid)
{
	UG_ASSERT(m_pGrid == grid, "Grids do not match in DistributedGridManager::unregistered_from_grid");
	
//	remove attached data
	grid->detach_from_vertices(m_aElemInfoVrt);
	grid->detach_from_edges(m_aElemInfoEdge);
	grid->detach_from_faces(m_aElemInfoFace);
	grid->detach_from_volumes(m_aElemInfoVol);
	m_pGrid = NULL;
	
//	clear the layout-map
	m_gridLayoutMap = GridLayoutMap();
}

////////////////////////////////////////////////////////////////////////
void DistributedGridManager::grid_layouts_changed(bool addedElemsOnly)
{
	if(!addedElemsOnly){
	//	first we have to reset all elem infos
		reset_elem_infos<VertexBase>();
		reset_elem_infos<EdgeBase>();
		reset_elem_infos<Face>();
		reset_elem_infos<Volume>();
	}
	
//	call for each layout in the grid-layout the corresponding
//	init_elem_status_from_layout function.
//	every layout has multiple levels

//	VERTICES
	update_elem_info<VertexBase>(m_gridLayoutMap, INT_MASTER,
								ES_IN_INTERFACE | ES_MASTER);

	update_elem_info<VertexBase>(m_gridLayoutMap, INT_SLAVE,
								ES_IN_INTERFACE | ES_SLAVE);

//	EDGES						
	update_elem_info<EdgeBase>(m_gridLayoutMap, INT_MASTER,
								ES_IN_INTERFACE | ES_MASTER);

	update_elem_info<EdgeBase>(m_gridLayoutMap, INT_SLAVE,
								ES_IN_INTERFACE | ES_SLAVE);
								
//	FACES
	update_elem_info<Face>(m_gridLayoutMap, INT_MASTER,
								ES_IN_INTERFACE | ES_MASTER);

	update_elem_info<Face>(m_gridLayoutMap, INT_SLAVE,
								ES_IN_INTERFACE | ES_SLAVE);

//	VOLUMES				
	update_elem_info<Volume>(m_gridLayoutMap, INT_MASTER,
								ES_IN_INTERFACE | ES_MASTER);

	update_elem_info<Volume>(m_gridLayoutMap, INT_SLAVE,
								ES_IN_INTERFACE | ES_SLAVE);}

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
update_elem_info(TLayoutMap& layoutMap, int nodeType, byte newStatus)
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
				int procID = layout.proc_id(iiter);
				
			///	iterate through the elements of the interface
				for(typename Layout::Interface::iterator iter = interface.begin();
					iter != interface.end(); ++iter)
				{
				//	set the new status
					//set_status(*iter, newStatus);
					elem_info(interface.get_element(iter)).set_status(newStatus);
					
				//	add the iterator to the iterator list and set the proc-id
					elem_info(interface.get_element(iter)).add_entry(&interface, iter);
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
get_status(GeometricObject* go)
{
	int baseType = go->base_object_type_id();
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
}

////////////////////////////////////////////////////////////////////////
//	begin_ / end_element_creation
void
DistributedGridManager::
begin_ordered_element_insertion()
{
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
		set_status(pElem, get_status(pElem) | ES_MASTER);
	else if(nodeType == pcl::NT_SLAVE)
		set_status(pElem, get_status(pElem) | ES_SLAVE);
		
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
	
	if(status & ES_MASTER){
		interface = &m_gridLayoutMap.get_layout<TElem>(INT_MASTER)
						.interface(procID, m_pGrid->get_level(pElem));
		iter = interface->push_back(pElem);
		elem_info(pElem).set_status(ES_IN_INTERFACE | ES_MASTER);
	}
	else{
		UG_ASSERT(status & ES_SLAVE, "interface-elements have to be either master or slave!");
		interface = &m_gridLayoutMap.get_layout<TElem>(INT_SLAVE)
						.interface(procID, m_pGrid->get_level(pElem));
		iter = interface->push_back(pElem);
						
		elem_info(pElem).set_status(ES_IN_INTERFACE | ES_SLAVE);
	}
	
//	add the interface-entry to the info
	elem_info(pElem).add_entry(interface, iter);
}
		
template <class TScheduledElemMap>
void DistributedGridManager::
perform_ordered_element_insertion(TScheduledElemMap& elemMap)
{		
	for(typename TScheduledElemMap::iterator iter = elemMap.begin();
		iter != elemMap.end(); ++iter)
	{
		ScheduledElement& schedElem = iter->second;
		int objType = schedElem.geomObj->base_object_type_id();
		switch(objType)
		{
			case VERTEX:
				add_element_to_interface((VertexBase*)schedElem.geomObj,
										schedElem.connectedProcID);
				break;
			case EDGE:
				add_element_to_interface((EdgeBase*)schedElem.geomObj,
										schedElem.connectedProcID);
				break;
//TODO: add support for faces and volumes
		}
	}
}
										
void
DistributedGridManager::
end_ordered_element_insertion()
{
//TODO: support all elements
	perform_ordered_element_insertion(m_vrtMap);
	perform_ordered_element_insertion(m_edgeMap);

	clear_scheduled_elements();
	m_bOrderedInsertionMode = false;
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

//	schedule one element for each parent-interface
	for(entry_iter iter = parentInfo.entries_begin();
		iter != parentInfo.entries_end(); ++iter)
	{
		UG_DLOG(LIB_GRID, 3, parentInfo.get_target_proc(iter) << ", ");
		elemMap.insert(make_pair(parentInfo.get_local_id(iter),
				ScheduledElement(elem, parentInfo.get_target_proc(iter))));
	}
	
//	set the status
	if(parentInfo.get_status() & ES_MASTER)
		elem_info(elem).set_status(ES_MASTER | ES_SCHEDULED_FOR_INTERFACE);
	else{
		UG_ASSERT(parentInfo.get_status() & ES_SLAVE, "interface-elements have to be either master or slave!");
		elem_info(elem).set_status(ES_SLAVE | ES_SCHEDULED_FOR_INTERFACE);
	}
}

//	new elements are handled here.
template <class TElem>
void DistributedGridManager::
handle_created_element(TElem* pElem,
						GeometricObject* pParent)
{
	elem_info(pElem).set_status(ES_NONE);
		
//	if there is no parent, we can immediatly leave.
	if(!pParent)
		return;
	
//	if the parent is not in an interface or scheduled to for an interface,
//	there is nothing to do either
	if(!(get_status(pParent) &
		(ES_IN_INTERFACE | ES_SCHEDULED_FOR_INTERFACE)))
		return;
		
	int parentType = pParent->base_object_type_id();

//	if ordered insertion mode is active, we have to insert the elements
//	into a map, instead of directly adding them to the group.
	if(m_bOrderedInsertionMode)
	{
		switch(parentType)
		{
//TODO: add support for all elements
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
		}
		elem_info(pElem).set_status(get_status(pElem) | ES_SCHEDULED_FOR_INTERFACE);
	}
	else
	{
	//	directly add the edge to the communication group
		assert(!"only sorted mode supported in the moment!");
	}
}
/*
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
vertex_created(Grid* grid, VertexBase* vrt, GeometricObject* pParent)
{
	handle_created_element(vrt, pParent);
}


void DistributedGridManager::
edge_created(Grid* grid, EdgeBase* e, GeometricObject* pParent)
{
	handle_created_element(e, pParent);
}



#ifdef __OLD_IMPLEMENTATION__
////////////////////////////////////////////////////////////////////////
//	GetElemInfo
///	returns element-info for the given element.
/*
template <class TGeomObj, class TCommGrp>
uint GetElemInfo(TGeomObj* pObj, TCommGrp& commGrp)
{
	typename TCommGrp::NodeHandleIterator hNode = commGrp.get_handle(pObj);
	if(commGrp.handle_is_valid(hNode))
	{
		uint retVal = 0;
		if(commGrp.num_containing_interfaces(hNode) > 0)
			retVal |= DistributedGrid::EI_INTERFACE_NODE;
		int nt = commGrp.get_node_type(hNode);
		switch(nt)
		{
			case pcl::NT_NORMAL: retVal |= DistributedGrid::EI_NORMAL; break;
			case pcl::NT_MASTER: retVal |= DistributedGrid::EI_MASTER; break;
			case pcl::NT_SLAVE: retVal |= DistributedGrid::EI_SLAVE; break;
		}
		return retVal;
	}
	return 0;
}
*/
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	constructor / destructor
DistributedGrid::
DistributedGrid()
{
	m_bOrderedInsertionMode = false;
	m_pGrid = NULL;
}

DistributedGrid::
DistributedGrid(Grid& grid)
{
	m_bOrderedInsertionMode = false;
	assign(grid);
}

DistributedGrid::		
~DistributedGrid()
{
}

////////////////////////////////////////////////////////////////////////
//	assignment
void
DistributedGrid::
assign(Grid& grid)
{
	grid.register_observer(this);
}

////////////////////////////////////////////////////////////////////////
//	begin_ / end_element_creation
void
DistributedGrid::
begin_ordered_element_insertion()
{
	m_bOrderedInsertionMode = true;
}

template <class TElem, class TCommGrp>
void DistributedGrid::
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
		set_status(pElem, get_status(pElem) | ES_MASTER);
	else if(nodeType == pcl::NT_SLAVE)
		set_status(pElem, get_status(pElem) | ES_SLAVE);
		
	set_status(pElem, get_status(pElem) | ES_IN_INTERFACE);

}

template <class TScheduledElemMap>
void DistributedGrid::
perform_ordered_element_insertion(TScheduledElemMap& elemMap)
{		
	for(typename TScheduledElemMap::iterator iter = elemMap.begin();
		iter != elemMap.end(); ++iter)
	{
		ScheduledElement& schedElem = iter->second;
		int objType = schedElem.pObj->base_object_type_id();
		switch(objType)
		{
			case VERTEX:
				add_element_to_interface((VertexBase*)schedElem.pObj,
										m_pCommSet->vrtGroup,
										schedElem.connectedProcID,
										schedElem.nodeType);
				break;
			case EDGE:
				add_element_to_interface((EdgeBase*)schedElem.pObj,
										m_pCommSet->edgeGroup,
										schedElem.connectedProcID,
										schedElem.nodeType);
				break;
//TODO: add support for faces and volumes
		}
	}
}
										
void
DistributedGrid::
end_ordered_element_insertion()
{
//TODO: support all elements
	perform_ordered_element_insertion(m_vrtMap);
	perform_ordered_element_insertion(m_edgeMap);

	clear_scheduled_elements();
	m_bOrderedInsertionMode = false;
}

////////////////////////////////////////////////////////////////////////
//	protected methods
void
DistributedGrid::
clear_scheduled_elements()
{
	m_vrtMap.clear();
	m_edgeMap.clear();
	m_faceMap.clear();
	m_volMap.clear();
}

/*
template <class TCommGrp, class TAttachmentAccessor>
void DistributedGrid::
init_elem_status_from_comm_group(TCommGrp& commGrp,
							TAttachmentAccessor& aaStatus)
{
	for(typename TCommGrp::NodeHandleIterator iter = commGrp.handles_begin();
		iter != commGrp.handles_end(); ++iter)
	{
		typename TCommGrp::HNODE hNode = *iter;
		byte& stat = aaStatus[commGrp.get_node(hNode)];
		stat = ES_IN_COMM_GRP;
		if(commGrp.num_containing_interfaces(hNode) > 0)
			stat |= ES_IN_INTERFACE;
		if(commGrp.get_node_type(hNode) == pcl::NT_MASTER)
			stat |= ES_MASTER;
		else if(commGrp.get_node_type(hNode) == pcl::NT_SLAVE)
			stat |= ES_SLAVE;
	}
}
*/

////////////////////////////////////////////////////////////////////////
void DistributedGrid::grid_layouts_changed(bool addedElemsOnly)
{
	if(!addedElemsOnly){
	//	first we have to reset all elem infos
		reset_elem_infos<VertexBase>();
		reset_elem_infos<EdgeBase>();
		reset_elem_infos<Face>();
		reset_elem_infos<Volume>();
	}
	
//	call for each layout in the grid-layout the corresponding
//	init_elem_status_from_layout function.
//	every layout has multiple levels
//	VERTICES
	if(m_gridLayoutMap.has_layout<VertexBase>(INT_MASTER))
		set_elem_statuses(m_gridLayoutMap.get_layout<VertexBase>(INT_MASTER),
						m_aaStatusVRT, m_aaIEIterVRT, ES_IN_INTERFACE | ES_MASTER);
	
	if(parallelGridLayout.has_vertex_layout(INT_SLAVE))
		set_elem_statuses(parallelGridLayout.vertex_layout(INT_SLAVE),
						m_aaStatusVRT, ES_IN_INTERFACE | ES_SLAVE);

//	EDGES						
	if(m_gridLayoutMap.has_edge_layout(INT_MASTER))
		set_elem_statuses(m_gridLayoutMap.edge_layout(INT_MASTER),
						m_aaStatusEDGE, ES_IN_INTERFACE | ES_MASTER);
	if(m_gridLayoutMap.has_edge_layout(INT_SLAVE))
		set_elem_statuses(m_gridLayoutMap.edge_layout(INT_SLAVE),
						m_aaStatusEDGE, ES_IN_INTERFACE | ES_SLAVE);

//	FACES
	if(m_gridLayoutMap.has_face_layout(INT_MASTER))
		set_elem_statuses(m_gridLayoutMap.face_layout(INT_MASTER),
						m_aaStatusFACE, ES_IN_INTERFACE | ES_MASTER);
	if(m_gridLayoutMap.has_face_layout(INT_SLAVE))
		set_elem_statuses(m_gridLayoutMap.face_layout(INT_SLAVE),
						m_aaStatusFACE, ES_IN_INTERFACE | ES_SLAVE);

//	VOLUMES				
	if(m_gridLayoutMap.has_volume_layout(INT_MASTER))
		set_elem_statuses(m_gridLayoutMap.volume_layout(INT_MASTER),
						m_aaStatusVOL, ES_IN_INTERFACE | ES_MASTER);
	if(m_gridLayoutMap.has_volume_layout(INT_SLAVE))
		set_elem_statuses(m_gridLayoutMap.volume_layout(INT_SLAVE),
						m_aaStatusVOL, ES_IN_INTERFACE | ES_SLAVE);						
}

////////////////////////////////////////////////////////////////////////
template <class TGeomObj>
void DistributedGrid::reset_elem_infos()
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
void DistributedGrid::
update_elem_info(TLayoutMap& layoutMap, int nodeType, byte newStatus)
{
	typedef typename TLayoutMap::template Types<TGeomObj>::Layout Layout;
	if(layoutMap.template has_layout<TGeomObj>(nodeType)){
	//	get the layout
		 Layout& layout = layoutMap.template get_layout<TGeomObj>(nodeType);
		 
	//	iterate through the levels of the layout
		for(size_t l = 0; l < layout.num_levels(); ++l){
		//	iterate through the interfaces of the layout
			for(Layout::iterator iiter = layout.begin(l);
				iiter != layout.end(l); ++iiter)
			{
				Layout::Interface& interface = layout.interface(iiter);
				int procID = layout.proc_id(iiter);
				
			///	iterate through the elements of the interface
				for(Layout::Interface::iterator iter = interface.begin();
					iter != interface.end(); ++iter)
				{
				//	set the new status
					set_status(*iter, newStatus);
				//	add the iterator to the iterator list and set the proc-id
					elem_info(*iter).lstProcIterPairs.push_back(
								typename ElementInfo<TGeomObj>(procID, iter));
				}
			}
		}
	}
}
/*
template <class TParallelElemLayout>
void DistributedGrid::
set_elem_statuses(TParallelElemLayout& pel,	byte newStatus)
{
	typedef typename TParallelElemLayout::Layout 	Layout
	typedef typename Layout::iterator				LayoutIter;
	typedef typename Layout::Interface				Interface;
	typedef typename Interface::iterator			InterfaceIter;
	
	for(size_t i = 0; i < pel.num_levels(); ++i){
		Layout& layout = pel.layout(i);
	//	iterate over all interfaces and assign status
		for(LayoutIter lIter = layout.begin();
			lIter != layout.end(); ++lIter)
		{
			Interface& interface = layout.interface(lIter);
			for(InterfaceIter iter = interface.begin();
				iter != interface.end(); ++iter)
			{
			//	value type of InterfaceIter is a GeometricObject-pointer.
				set_status(*iter, newStatus);
			}
		}
	}
}
*/
byte DistributedGrid::
get_geom_obj_status(GeometricObject* go)
{
	int baseType = go->base_object_type_id();
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
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	implementation of observer-methods.
//	grid callbacks
void
DistributedGrid::
registered_at_grid(Grid* grid)
{
	if(m_pGrid != NULL)
		m_pGrid->unregister_observer(this);

	m_pGrid = grid;
	
//	attach element infos
	m_pGrid->attach_to_vertices(m_aElemInfoVrt);
	m_pGrid->attach_to_edges(m_aElemInfoEdge);
	m_pGrid->attach_to_faces(m_aElemInfoFace);
	m_pGrid->attach_to_volumes(m_aElemInfoVol);
	
//	access them
	m_aaElemInfoVRT.access(*grid, m_aElemInfoVrt);
	m_aaElemInfoEDGE.access(*grid, m_aElemInfoEdge);
	m_aaElemInfoFACE.access(*grid, m_aElemInfoFace);
	m_aaElemInfoVOL.access(*grid, m_aElemInfoVol);
	
//	initialise the element statuses. This is automatically done
//	on a call to grid_layout_changed.
	grid_layout_changed();
}

void
DistributedGrid::
unregistered_from_grid(Grid* grid)
{
	if(m_bOrderedInsertionMode)
	{
		m_bOrderedInsertionMode = false;
		clear_scheduled_elements();
	}
	
	m_pGrid->detach_from_vertices(m_aElemInfoVrt);
	m_pGrid->detach_from_edges(m_aElemInfoEdge);
	m_pGrid->detach_from_faces(m_aElemInfoFace);
	m_pGrid->detach_from_volumes(m_aElemInfoVol);
	
	m_pGrid = NULL;
}

void
DistributedGrid::
elements_to_be_cleared(Grid* grid)
{
	if(m_bOrderedInsertionMode)
		clear_scheduled_elements();
}

//	adds a ScheduledElement to scheduled-element-maps, depending on
//	the interfaces in which pParent is contained.
template <class TElem, class TScheduledElemMap, class TParent>
void ScheduleElementForCommunicationGroup(TScheduledElemMap& elemMap,
										TElem* elem,
										TParent* pParent)
{
	typedef typename ElementInfo<TElem>::ProcIterPairList TList;
	TList& lst = elem_info(elem).lstProcIterPairs;
	
	for(typename TList::iterator iter = lst.begin(); iter != lst.end(); ++iter)
	{
		elemMap.insert(make_pair(vEntryInds[i],
			DistributedGrid::ScheduledElement(pObj,
				parentCommGrp.get_node_type(hParent),
				parentCommGrp.get_connected_proc_of_interface(vInterfaceInds[i]))));
	}
}

//	new elements are handled here.
template <class TElem>
void DistributedGrid::
handle_created_element(TElem* pElem,
						GeometricObject* pParent)
{
	set_status(pElem, ES_NONE);
	
//	if there is no parent, we can immediatly leave.
	if(!pParent)
		return;
		
	int parentType = pParent->base_object_type_id();
	
//	if the parent is not in an interface or scheduled to for an interface,
//	there is nothing to do either
	if(!(get_geom_obj_status(pParent) &
		(ES_IN_INTERFACE | ES_SHEDULED_FOR_INTERFACE)))
		return;
		
//	if ordered insertion mode is active, we have to insert the elements
//	into a map, instead of directly adding them to the group.
	if(m_bOrderedInsertionMode)
	{
		switch(parentType)
		{
//TODO: add support for all elements
			case EDGE:
				ScheduleElementForCommunicationGroup(m_edgeMap,
													pElem,
													(EdgeBase*)pParent);
				break;
		}
		set_status(pElem, get_status(pElem) | ES_SCHEDULED_FOR_INTERFACE);
	}
	else
	{
	//	directly add the edge to the communication group
		assert(!"only sorted mode supported in the moment!");
	}
}

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

void DistributedGrid::
vertex_created(Grid* grid, VertexBase* vrt, GeometricObject* pParent)
{
	handle_created_element(vrt, pParent);
}

void DistributedGrid::
vertex_to_be_erased(Grid* grid, VertexBase* vrt)
{
	handle_erased_element(vrt, m_pCommSet->vrtGroup);
}

void DistributedGrid::
vertex_to_be_replaced(Grid* grid, VertexBase* vrtOld, VertexBase* vrtNew)
{
	handle_replaced_element(vrtOld, vrtNew, m_pCommSet->vrtGroup);
}

void DistributedGrid::edge_created(Grid* grid, EdgeBase* edge,
											GeometricObject* pParent)
{
	handle_created_element(edge, pParent);
}

void DistributedGrid::edge_to_be_erased(Grid* grid, EdgeBase* edge)
{
	handle_erased_element(edge, m_pCommSet->edgeGroup);
}

void DistributedGrid::edge_to_be_replaced(Grid* grid, EdgeBase* edgeOld,
													EdgeBase* edgeNew)
{
	handle_replaced_element(edgeOld, edgeNew, m_pCommSet->edgeGroup);
}

#endif __OLD_IMPLEMENTATION__
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	implementation of protected methods.

}// end of namespace

