// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y10 m02 d15

#include "lib_grid/algorithms/attachment_util.h"
#include "selector_interface.h"
#ifdef UG_PARALLEL
	#include "lib_grid/parallelization/distributed_grid.h"
	#include "lib_grid/parallelization/util/compol_selection.h"
#endif


namespace ug
{
	
ISelector::ISelector(uint supportedElements) :
	m_aSelected("ISelector_IsSelected")
{
	m_pGrid = NULL;
	m_supportedElements = supportedElements;// Since no grid is available, we can't yet activate them.
	m_bAutoselectionEnabled = false;
	m_bSelectionInheritanceEnabled = true;
	m_bStrictInheritanceEnabled = false;
}

ISelector::ISelector(Grid& grid, uint supportedElements) :
	m_aSelected("ISelector_IsSelected")
{
	m_pGrid = &grid;
	m_supportedElements = SE_NONE;
	m_bAutoselectionEnabled = false;
	m_bSelectionInheritanceEnabled = true;
	m_bStrictInheritanceEnabled = false;

//	register at grid. Don't use set_grid here since it invokes virtual methods.
	if(m_pGrid){
		m_pGrid->register_observer(this, OT_GRID_OBSERVER | OT_VERTEX_OBSERVER | OT_EDGE_OBSERVER |
									OT_FACE_OBSERVER | OT_VOLUME_OBSERVER);
									
	//	initialise attachments and accessors.
		enable_element_support(supportedElements);
	}
}

ISelector::~ISelector()
{
//	unregister from grid
//	don't use set_grid here, since it invokes virtual methods.
	if(m_pGrid){
	//	disable all currently supported elements (this will remove any attachments)
		disable_element_support(m_supportedElements);
		m_pGrid->unregister_observer(this);
		m_pGrid = NULL;
	}
}

void ISelector::
set_supported_elements(uint shElements)
{
//	do this in two steps:
//	1: disable the element-support that is no longer required.
//	2: enable the element-support that was not already enabled.
//	disable the elements that shall be disabled.

//	(the ones which shall not be set, but are currently active.)
	disable_element_support((~shElements) & m_supportedElements);

//	enable the elements that are not already enabled
	enable_element_support(shElements & (~m_supportedElements));
}

void ISelector::enable_element_support(uint shElements)
{
//	if no grid is assigned, we can't do anything.
	if(m_pGrid){
	//	check for each option whether it should be enabled.
	//	to reduce unnecessary operations, we have to make sure that
	//	that option hasn't already been enabled.

		if((shElements & SE_VERTEX) &&
			(!elements_are_supported(SE_VERTEX))){
//LOG("enabling vertex support\n");
		//	enable vertex-support.
			m_pGrid->attach_to_vertices(m_aSelected);
			m_aaSelVRT.access(*m_pGrid, m_aSelected);
			SetAttachmentValues(m_aaSelVRT, m_pGrid->begin<VertexBase>(),
								m_pGrid->end<VertexBase>(), 0);
			m_supportedElements |= SE_VERTEX;
		}

		if((shElements & SE_EDGE) &&
			(!elements_are_supported(SE_EDGE))){
//LOG("enabling edge support\n");
		//	enable edge support
			m_pGrid->attach_to_edges(m_aSelected);
			m_aaSelEDGE.access(*m_pGrid, m_aSelected);
			SetAttachmentValues(m_aaSelEDGE, m_pGrid->begin<EdgeBase>(),
								m_pGrid->end<EdgeBase>(), 0);
			m_supportedElements |= SE_EDGE;
		}

		if((shElements & SE_FACE) &&
			(!elements_are_supported(SE_FACE))){
//LOG("enabling face support\n");
		//	enable face support
			m_pGrid->attach_to_faces(m_aSelected);
			m_aaSelFACE.access(*m_pGrid, m_aSelected);
			SetAttachmentValues(m_aaSelFACE, m_pGrid->begin<Face>(),
								m_pGrid->end<Face>(), 0);
			m_supportedElements |= SE_FACE;
		}

		if((shElements & SE_VOLUME) &&
			(!elements_are_supported(SE_VOLUME))){
//LOG("enabling volume support\n");
		//	enable volume support
			m_pGrid->attach_to_volumes(m_aSelected);
			m_aaSelVOL.access(*m_pGrid, m_aSelected);
			SetAttachmentValues(m_aaSelVOL, m_pGrid->begin<Volume>(),
								m_pGrid->end<Volume>(), 0);
			m_supportedElements |= SE_VOLUME;
		}
	}
	else{
		m_supportedElements |= shElements;
	}
}

void ISelector::disable_element_support(uint shElements)
{
//	if no grid is assigned, we can't do anything.
	if(m_pGrid){
	//	check for each option whether it should be disabled.
	//	to reduce unnecessary operations, we have to make sure that
	//	that option hasn't already been disabled.

		if((shElements & SE_VERTEX) && elements_are_supported(SE_VERTEX)){
//LOG("disabling vertex support\n");
			m_pGrid->detach_from_vertices(m_aSelected);
		}

		if((shElements & SE_EDGE) && elements_are_supported(SE_EDGE)){
//LOG("disabling edge support\n");
			m_pGrid->detach_from_edges(m_aSelected);
		}

		if((shElements & SE_FACE) && elements_are_supported(SE_FACE)){
//LOG("disabling face support\n");
			m_pGrid->detach_from_faces(m_aSelected);
		}

		if((shElements & SE_VOLUME) && elements_are_supported(SE_VOLUME)){
//LOG("disabling volume support\n");
			m_pGrid->detach_from_volumes(m_aSelected);
		}
	}

//	remove the disabled elements from the set of currently supported elements.
	m_supportedElements &= (~shElements);
}

void ISelector::enable_autoselection(bool bEnable)
{
	m_bAutoselectionEnabled = bEnable;
}

void ISelector::enable_selection_inheritance(bool bEnable)
{
	m_bSelectionInheritanceEnabled = bEnable;
}

void ISelector::enable_strict_inheritance(bool bEnable)
{
	m_bStrictInheritanceEnabled = bEnable;
}


void ISelector::set_grid(Grid* grid)
{
//	if we're already registered at this grid then return
	if(m_pGrid == grid)
		return;
		
//	if we're already registered at a grid unregister first.
	if(m_pGrid){
	//	disable all currently supported elements (this will remove any attachments)
		clear();
		disable_element_support(m_supportedElements);
		m_pGrid->unregister_observer(this);
		m_pGrid = NULL;
	}

//	if the new grid is not empty, we'll initialise and register
	if(grid){
		grid->register_observer(this, OT_GRID_OBSERVER | OT_VERTEX_OBSERVER | OT_EDGE_OBSERVER |
									OT_FACE_OBSERVER | OT_VOLUME_OBSERVER);
		m_pGrid = grid;

	//	initialise attachments and accessors.
	//	do this whith a little trick:
	//	set the supported-element-options to SE_NONE,
	//	then call enable for all element-types that should be supported.
		uint tmpOpts = m_supportedElements;
		m_supportedElements = SE_NONE;
		enable_element_support(tmpOpts);
	}
}


#ifdef UG_PARALLEL
template <class TIntfcCom>
void ISelector::broadcast_selection_states(bool deselect,
										   bool includeGhosts,
										   TIntfcCom& icom)
{
	DistributedGridManager& dgm = *m_pGrid->distributed_grid_manager();
	GridLayoutMap& glm = dgm.grid_layout_map();

	if(deselect){
		ComPol_Selection<typename TIntfcCom::Layout> compol(*this, false, true);
		if(includeGhosts){
		//	if ghosts shall be included, we will first have to make sure that
		//	copies in h-interfaces know about the selection state of associated
		//	ghosts
			icom.exchange_data(glm, INT_V_MASTER, INT_V_SLAVE, compol);
			icom.communicate();
		}

	//	gather selection state at h-master
		icom.exchange_data(glm, INT_H_SLAVE, INT_H_MASTER, compol);
		icom.communicate();

	//	copy it to all h-slaves
		icom.exchange_data(glm, INT_H_MASTER, INT_H_SLAVE, compol);
		icom.communicate();

		if(includeGhosts){
		//	if ghosts shall be included, we will now have to copy the values back
		//	to those ghosts
			icom.exchange_data(glm, INT_V_SLAVE, INT_V_MASTER, compol);
			icom.communicate();
		}
	}
	else{
		ComPol_EnableSelectionStateBits<typename TIntfcCom::Layout> compol(*this, 0xFF);
		if(includeGhosts){
		//	if ghosts shall be included, we will first have to make sure that
		//	copies in h-interfaces know about the selection state of associated
		//	ghosts
			icom.exchange_data(glm, INT_V_MASTER, INT_V_SLAVE, compol);
			icom.communicate();
		}

	//	gather selection state at h-master
		icom.exchange_data(glm, INT_H_SLAVE, INT_H_MASTER, compol);
		icom.communicate();

	//	copy it to all h-slaves
		icom.exchange_data(glm, INT_H_MASTER, INT_H_SLAVE, compol);
		icom.communicate();

		if(includeGhosts){
		//	if ghosts shall be included, we will now have to copy the values back
		//	to those ghosts
			icom.exchange_data(glm, INT_V_SLAVE, INT_V_MASTER, compol);
			icom.communicate();
		}
	}
}
#endif

void ISelector::broadcast_selection_states(bool deselect,
										   bool includeGhosts)
{
#ifdef UG_PARALLEL
	if(!m_pGrid)
		return;

	broadcast_selection_states(deselect, includeGhosts, m_icomVRT);
	broadcast_selection_states(deselect, includeGhosts, m_icomEDGE);
	broadcast_selection_states(deselect, includeGhosts, m_icomFACE);
	broadcast_selection_states(deselect, includeGhosts, m_icomVOL);
#endif
}

////////////////////////////////////////////////////////////////////////
//	grid callbacks
/*
void ISelector::registered_at_grid(Grid* grid)
{
//	if we're already registered at this grid then return
	if(m_pGrid == grid)
		return;

//	if we're already registered at a grid, then unregister first
	if(m_pGrid)
		m_pGrid->unregister_observer(this);

//	assign grid
	m_pGrid = grid;

//	initialise attachments and accessors.
//	do this whith a little trick:
//	set the supported-element-options to SE_NONE,
//	then call enable for all element-types that should be supported.
	uint tmpOpts = m_supportedElements;
	m_supportedElements = SE_NONE;
	enable_element_support(tmpOpts);
}

void ISelector::unregistered_from_grid(Grid* grid)
{
	assert(m_pGrid == grid && "grids do not match!");

	if(m_pGrid == grid){
	//	disable all currently supported elements (this will remove any attachments)
		disable_element_support(m_supportedElements);
		m_pGrid = NULL;
	}
}
*/
void ISelector::grid_to_be_destroyed(Grid* grid)
{
	assert(m_pGrid == grid && "grids do not match!");

	if(m_pGrid == grid){
		set_grid(NULL);
	}
}

void ISelector::elements_to_be_cleared(Grid* grid)
{
	clear();
}

//	vertex callbacks
void ISelector::vertex_created(Grid* grid, VertexBase* vrt,
								GeometricObject* pParent,
								bool replacesParent)
{
	assert((m_pGrid == grid) && "grids do not match.");
	
//TODO: this if could be removed if the subset-handler was only registered for
//		the elements that it supports. Note that a dynamic register/unregister
//		would be required...
	if(elements_are_supported(SE_VERTEX)){
	//	init the element
		mark_deselected(vrt);
		if(autoselection_enabled())
			select(vrt);
		else if((pParent != NULL) && selection_inheritance_enabled()){
			if(m_bStrictInheritanceEnabled){
				if(pParent->base_object_id() == VERTEX){
					select(vrt, get_selection_status(static_cast<VertexBase*>(pParent)));
				}
			}
			else
				select(vrt, get_selection_status(pParent));
		}
		else if(replacesParent){
			UG_ASSERT(pParent, "A parent has to exist if it shall be replaced");
			UG_ASSERT(dynamic_cast<VertexBase*>(pParent), "Only parents of the same type may be replaced.");
			select(vrt, get_selection_status(static_cast<VertexBase*>(pParent)));
		}
	}
}

void ISelector::vertex_to_be_erased(Grid* grid, VertexBase* vrt,
									 VertexBase* replacedBy)
{
	assert((m_pGrid == grid) && "grids do not match.");
	
//TODO: this if could be removed if the subset-handler was only registered for
//		the elements that it supports. Note that a dynamic register/unregister
//		would be required...
	if(elements_are_supported(SE_VERTEX)){
		deselect(vrt);
	}
}

//	edge callbacks
void ISelector::edge_created(Grid* grid, EdgeBase* edge,
							GeometricObject* pParent,
							bool replacesParent)
{
	assert((m_pGrid == grid) && "grids do not match.");
	
//TODO: this if could be removed if the subset-handler was only registered for
//		the elements that it supports. Note that a dynamic register/unregister
//		would be required...
	if(elements_are_supported(SE_EDGE)){
	//	init the element
		mark_deselected(edge);
		if(autoselection_enabled())
			select(edge);
		else if((pParent != NULL) && selection_inheritance_enabled()){
			if(m_bStrictInheritanceEnabled){
				if(pParent->base_object_id() == EDGE){
					select(edge, get_selection_status(static_cast<EdgeBase*>(pParent)));
				}
			}
			else
				select(edge, get_selection_status(pParent));
		}
		else if(replacesParent){
			UG_ASSERT(pParent, "A parent has to exist if it shall be replaced");
			UG_ASSERT(dynamic_cast<EdgeBase*>(pParent), "Only parents of the same type may be replaced.");
			select(edge, get_selection_status(static_cast<EdgeBase*>(pParent)));
		}
	}
}

void ISelector::edge_to_be_erased(Grid* grid, EdgeBase* edge,
									EdgeBase* replacedBy)
{
	assert((m_pGrid == grid) && "grids do not match.");
	
//TODO: this if could be removed if the subset-handler was only registered for
//		the elements that it supports. Note that a dynamic register/unregister
//		would be required...
	if(elements_are_supported(SE_EDGE)){
		deselect(edge);
	}
}

//	face callbacks
void ISelector::face_created(Grid* grid, Face* face,
							GeometricObject* pParent,
							bool replacesParent)
{
	assert((m_pGrid == grid) && "grids do not match.");
	
//TODO: this if could be removed if the subset-handler was only registered for
//		the elements that it supports. Note that a dynamic register/unregister
//		would be required...
	if(elements_are_supported(SE_FACE)){
	//	init the element
		mark_deselected(face);
		if(autoselection_enabled())
			select(face);
		else if((pParent != NULL) && selection_inheritance_enabled()){
			if(m_bStrictInheritanceEnabled){
				if(pParent->base_object_id() == FACE){
					select(face, get_selection_status(static_cast<Face*>(pParent)));
				}
			}
			else
				select(face, get_selection_status(pParent));
		}
		else if(replacesParent){
			UG_ASSERT(pParent, "A parent has to exist if it shall be replaced");
			UG_ASSERT(dynamic_cast<Face*>(pParent), "Only parents of the same type may be replaced.");
			select(face, get_selection_status(static_cast<Face*>(pParent)));
		}
	}
}

void ISelector::face_to_be_erased(Grid* grid, Face* face,
								 Face* replacedBy)
{
	assert((m_pGrid == grid) && "grids do not match.");
	
//TODO: this if could be removed if the subset-handler was only registered for
//		the elements that it supports. Note that a dynamic register/unregister
//		would be required...
	if(elements_are_supported(SE_FACE)){
		deselect(face);
	}
}

//	volume callbacks
void ISelector::volume_created(Grid* grid, Volume* vol,
								GeometricObject* pParent,
								bool replacesParent)
{
	assert((m_pGrid == grid) && "grids do not match.");
	
//TODO: this if could be removed if the subset-handler was only registered for
//		the elements that it supports. Note that a dynamic register/unregister
//		would be required...
	if(elements_are_supported(SE_VOLUME)){
	//	init the element
		mark_deselected(vol);
		if(autoselection_enabled())
			select(vol);
		else if((pParent != NULL) && selection_inheritance_enabled()){
			if(m_bStrictInheritanceEnabled){
				if(pParent->base_object_id() == VOLUME){
					select(vol, get_selection_status(static_cast<Volume*>(pParent)));
				}
			}
			else
				select(vol, get_selection_status(pParent));
		}
		else if(replacesParent){
			UG_ASSERT(pParent, "A parent has to exist if it shall be replaced");
			UG_ASSERT(dynamic_cast<Volume*>(pParent), "Only parents of the same type may be replaced.");
			select(vol, get_selection_status(static_cast<Volume*>(pParent)));
		}
	}
}

void ISelector::volume_to_be_erased(Grid* grid, Volume* vol,
									 Volume* replacedBy)
{
	assert((m_pGrid == grid) && "grids do not match.");
	
//TODO: this if could be removed if the subset-handler was only registered for
//		the elements that it supports. Note that a dynamic register/unregister
//		would be required...
	if(elements_are_supported(SE_VOLUME)){
		deselect(vol);
	}
}


template <class TElem>
void ISelector::
elems_to_be_merged(Grid* grid, TElem* target,
					TElem* elem1, TElem* elem2)
{
//	if at least one was selected, we'll select the new one too.
	if(is_selected(elem1) || is_selected(elem2))
		select(target);
}

void ISelector::
vertices_to_be_merged(Grid* grid, VertexBase* target,
					 VertexBase* elem1, VertexBase* elem2)
{
	elems_to_be_merged(grid, target, elem1, elem2);
}

void ISelector::
edges_to_be_merged(Grid* grid, EdgeBase* target,
				  EdgeBase* elem1, EdgeBase* elem2)
{
	elems_to_be_merged(grid, target, elem1, elem2);
}

void ISelector::
faces_to_be_merged(Grid* grid, Face* target,
					Face* elem1, Face* elem2)
{
	elems_to_be_merged(grid, target, elem1, elem2);
}

void ISelector::
volumes_to_be_merged(Grid* grid, Volume* target,
					Volume* elem1, Volume* elem2)
{
	elems_to_be_merged(grid, target, elem1, elem2);
}

}//	end of namespace
