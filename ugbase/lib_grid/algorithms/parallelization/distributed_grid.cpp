// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y09 m08 d17

#include "distributed_grid_observer.h"

using namespace std;

namespace libGrid
{

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
			retVal |= DistributedGridObserver::EI_INTERFACE_NODE;
		int nt = commGrp.get_node_type(hNode);
		switch(nt)
		{
			case pcl::NT_NORMAL: retVal |= DistributedGridObserver::EI_NORMAL; break;
			case pcl::NT_MASTER: retVal |= DistributedGridObserver::EI_MASTER; break;
			case pcl::NT_SLAVE: retVal |= DistributedGridObserver::EI_SLAVE; break;
		}
		return retVal;
	}
	return 0;
}
*/
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	constructor / destructor
DistributedGridObserver::
DistributedGridObserver()
{
	m_bOrderedInsertionMode = false;
	m_pGrid = NULL;
}

DistributedGridObserver::
DistributedGridObserver(Grid& grid)
{
	m_bOrderedInsertionMode = false;
	assign(grid);
}

DistributedGridObserver::		
~DistributedGridObserver()
{
}

////////////////////////////////////////////////////////////////////////
//	assignment
void
DistributedGridObserver::
assign(Grid& grid)
{
	grid.register_observer(this);
}

////////////////////////////////////////////////////////////////////////
//	begin_ / end_element_creation
void
DistributedGridObserver::
begin_ordered_element_insertion()
{
	m_bOrderedInsertionMode = true;
}

template <class TElem, class TCommGrp>
void DistributedGridObserver::
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
void DistributedGridObserver::
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
DistributedGridObserver::
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
DistributedGridObserver::
clear_scheduled_elements()
{
	m_vrtMap.clear();
	m_edgeMap.clear();
	m_faceMap.clear();
	m_volMap.clear();
}

/*
template <class TCommGrp, class TAttachmentAccessor>
void DistributedGridObserver::
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
void DistributedGridObserver::grid_layouts_changed(bool addedElemsOnly)
{
	if(!addedElemsOnly){
	//	first we have to clear all status-entries
		SetAttachmentValues(m_aaStatusVRT, m_pGrid->vertices_begin(),
							m_pGrid->vertices_end(), ES_NONE);
		SetAttachmentValues(m_aaStatusEDGE, m_pGrid->edges_begin(),
							m_pGrid->edges_end(), ES_NONE);
		SetAttachmentValues(m_aaStatusFACE, m_pGrid->faces_begin(),
							m_pGrid->faces_end(), ES_NONE);
		SetAttachmentValues(m_aaStatusVOL, m_pGrid->volumes_begin(),
							m_pGrid->volumes_end(), ES_NONE);
	}
	
//	call for each layout in the grid-layout the corresponding
//	init_elem_status_from_layout function.
//	every layout has multiple levels
//	VERTICES
	if(m_gridLayoutMap.has_vertex_layout(INT_MASTER))
		set_elem_statuses(m_gridLayoutMap.vertex_layout(INT_MASTER),
						m_aaStatusVRT, ES_IN_INTERFACE | ES_MASTER);
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
template <class TParallelElemLayout, class TAttachmentAccessor>
void DistributedGridObserver::
set_elem_statuses(TParallelElemLayout& pel, TAttachmentAccessor& aaStatus,
					byte newStatus)
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
				aaStatus[*iter] = newStatus;
			}
		}
	}
}

byte DistributedGridObserver::
get_geom_obj_status(GeometricObject* go)
{
	int baseType = go->base_object_type_id();
	switch(baseType)
	{
		case VERTEX:
			return m_aaStatusVRT[go];
		case EDGE:
			return m_aaStatusEDGE[go];
		case FACE:
			return m_aaStatusFACE[go];
		case VOLUME:
			return m_aaStatusVOL[go];
	}
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	implementation of observer-methods.
//	grid callbacks
void
DistributedGridObserver::
registered_at_grid(Grid* grid)
{
	if(m_pGrid != NULL)
		m_pGrid->unregister_observer(this);

	m_pGrid = grid;
	m_pGrid->attach_to_vertices(m_aStatus);
	m_pGrid->attach_to_edges(m_aStatus);
	m_pGrid->attach_to_faces(m_aStatus);
	m_pGrid->attach_to_volumes(m_aStatus);

	m_aaStatusVRT.access(*grid, m_aStatus);
	m_aaStatusEDGE.access(*grid, m_aStatus);
	m_aaStatusFACE.access(*grid, m_aStatus);
	m_aaStatusVOL.access(*grid, m_aStatus);
	
//	initialise the element statuses. This is automatically done
//	on a call to grid_layout_changed.
	grid_layout_changed();
}

void
DistributedGridObserver::
unregistered_from_grid(Grid* grid)
{
	if(m_bOrderedInsertionMode)
	{
		m_bOrderedInsertionMode = false;
		clear_scheduled_elements();
	}
	
	m_pGrid->detach_from_vertices(m_aStatus);
	m_pGrid->detach_from_edges(m_aStatus);
	m_pGrid->detach_from_faces(m_aStatus);
	m_pGrid->detach_from_volumes(m_aStatus);
	
	m_pGrid = NULL;
}

void
DistributedGridObserver::
elements_to_be_cleared(Grid* grid)
{
	if(m_bOrderedInsertionMode)
		clear_scheduled_elements();
}

//	adds a ScheduledElement to scheduled-element-maps, depending on
//	the interfaces in which pParent is contained.
template <class TScheduledElemMap, class TParent, class TCommGrp>
void ScheduleElementForCommunicationGroup(TScheduledElemMap& elemMap,
										GeometricObject* pObj,
										TParent* pParent,
										TCommGrp& parentCommGrp)
{
	typename TCommGrp::HNODE hParent = parentCommGrp.get_handle(pParent);

//	insert the element into the maps
	const vector<int>& vInterfaceInds = parentCommGrp.get_interface_indices(hParent);
	const vector<int>& vEntryInds = parentCommGrp.get_interface_entry_indices(hParent);
	
	for(int i = 0; i < vInterfaceInds.size(); ++i)
	{
		elemMap.insert(make_pair(vEntryInds[i],
			DistributedGridObserver::ScheduledElement(pObj,
				parentCommGrp.get_node_type(hParent),
				parentCommGrp.get_connected_proc_of_interface(vInterfaceInds[i]))));
	}
}

//	new elements are handled here.
template <class TElem>
void DistributedGridObserver::
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
													(EdgeBase*)pParent,
													m_pCommSet->edgeGroup);
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
void DistributedGridObserver::
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
void DistributedGridObserver::
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

void DistributedGridObserver::
vertex_created(Grid* grid, VertexBase* vrt, GeometricObject* pParent)
{
	handle_created_element(vrt, pParent);
}

void DistributedGridObserver::
vertex_to_be_erased(Grid* grid, VertexBase* vrt)
{
	handle_erased_element(vrt, m_pCommSet->vrtGroup);
}

void DistributedGridObserver::
vertex_to_be_replaced(Grid* grid, VertexBase* vrtOld, VertexBase* vrtNew)
{
	handle_replaced_element(vrtOld, vrtNew, m_pCommSet->vrtGroup);
}

void DistributedGridObserver::edge_created(Grid* grid, EdgeBase* edge,
											GeometricObject* pParent)
{
	handle_created_element(edge, pParent);
}

void DistributedGridObserver::edge_to_be_erased(Grid* grid, EdgeBase* edge)
{
	handle_erased_element(edge, m_pCommSet->edgeGroup);
}

void DistributedGridObserver::edge_to_be_replaced(Grid* grid, EdgeBase* edgeOld,
													EdgeBase* edgeNew)
{
	handle_replaced_element(edgeOld, edgeNew, m_pCommSet->edgeGroup);
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	implementation of protected methods.

}// end of namespace

