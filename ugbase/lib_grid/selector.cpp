// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y08 m11 d20

#include <cassert>
#include "selector.h"
#include "common/common.h"
#include "multi_grid.h"

using namespace std;

namespace ug
{
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	implementation of GenericElementSelector
template <class TElem, class SelectionPolicy>
GenericElementSelector<TElem, SelectionPolicy>::
GenericElementSelector()
{
	m_bAutoselectionEnabled = true;
	m_pGrid = NULL;
}

template <class TElem, class SelectionPolicy>
GenericElementSelector<TElem, SelectionPolicy>::
GenericElementSelector(TGridRef grid)
{
	m_bAutoselectionEnabled = true;
	m_pGrid = NULL;
	assign_grid(grid);
}

template <class TElem, class SelectionPolicy>
GenericElementSelector<TElem, SelectionPolicy>::
GenericElementSelector(const GenericElementSelector<TElem, SelectionPolicy>& gSel)
{
	assert(!"WARNING in GenericElementSelector::GenericElementSelector(const GenericElementSelector& sel): Copy-Constructor not yet implemented!");
	LOG("WARNING in GenericElementSelector::GenericElementSelector(const GenericElementSelector& sel): Copy-Constructor not yet implemented! Expect unexpected behaviour!");
}

template <class TElem, class SelectionPolicy>
GenericElementSelector<TElem, SelectionPolicy>::~GenericElementSelector()
{
	if(m_pGrid)
	{
	//	unregister the observer.
		m_pGrid->unregister_observer(this);
	}
}

template <class TElem, class SelectionPolicy>
void
GenericElementSelector<TElem, SelectionPolicy>::
assign_grid(TGridRef grid)
{
//	get the observer type
	uint observerType = OT_GRID_OBSERVER;

	int elemType = geometry_traits<TElem>::BASE_OBJECT_TYPE_ID;
	switch(elemType)
	{
		case VERTEX:	observerType |= OT_VERTEX_OBSERVER; break;
		case EDGE:		observerType |= OT_EDGE_OBSERVER; break;
		case FACE:		observerType |= OT_FACE_OBSERVER; break;
		case VOLUME:	observerType |= OT_VOLUME_OBSERVER; break;
	}

	grid.register_observer(this, observerType);
}

template <class TElem, class SelectionPolicy>
typename GenericElementSelector<TElem, SelectionPolicy>::TGridPtr
GenericElementSelector<TElem, SelectionPolicy>::
get_assigned_grid()
{
	return m_pGrid;
}


template <class TElem, class SelectionPolicy>
void
GenericElementSelector<TElem, SelectionPolicy>::
enable_autoselection(bool bEnable)
{
	m_bAutoselectionEnabled = bEnable;
}

template <class TElem, class SelectionPolicy>
void
GenericElementSelector<TElem, SelectionPolicy>::
enable_selection_inheritance(bool bEnable)
{
	m_bSelectionInheritanceEnabled = bEnable;
}
		
template <class TElem, class SelectionPolicy>
void
GenericElementSelector<TElem, SelectionPolicy>::
elem_created(Grid* grid, TElem* elem, GeometricObject* pParent)
{
	assert((m_pGrid == grid) && "ERROR in GenericElementSelector::elem_created(...): selectors grid and elements grid do not match!");
//	init deselected.
	SelectionPolicy::init_element(elem);

	if(autoselection_enabled() == true)
		SelectionPolicy::select(elem);
	else if(selection_inheritance_enabled() == true)
	{
		if(pParent != NULL)
		{
			if(SelectionPolicy::is_selected(pParent))
				SelectionPolicy::select(elem);
		}
	}
}

template <class TElem, class SelectionPolicy>
void
GenericElementSelector<TElem, SelectionPolicy>::
elem_to_be_erased(Grid* grid, TElem* elem)
{
	assert((m_pGrid == grid) && "ERROR in GenericElementSelector::elem_to_be_erased(...): selectors grid and elements grid do not match!");

//TODO: It is sufficient to erase the element from the list, if it was selected.
//		A function that handles this only task could be added to the SelectionPolicy.
//		This should lead to a slightly improved performance.
	SelectionPolicy::deselect(elem);
}

//	grid callbacks
template <class TElem, class SelectionPolicy>
void
GenericElementSelector<TElem, SelectionPolicy>::
registered_at_grid(Grid* grid)
{
	if(m_pGrid != NULL)
		m_pGrid->unregister_observer(this);

	m_pGrid = dynamic_cast<TGridPtr>(grid);

	SelectionPolicy::new_grid(m_pGrid);
}

template <class TElem, class SelectionPolicy>
void
GenericElementSelector<TElem, SelectionPolicy>::
unregistered_from_grid(Grid* grid)
{
	SelectionPolicy::new_grid(NULL);
}

template <class TElem, class SelectionPolicy>
void
GenericElementSelector<TElem, SelectionPolicy>::
elements_to_be_cleared(Grid* grid)
{
	assert((grid == m_pGrid) && "ERROR in GenericElementSelector::elements_to_be_cleared(...): grids do not match!");
	clear_selection();
}

//	vertex callbacks
template <class TElem, class SelectionPolicy>
void
GenericElementSelector<TElem, SelectionPolicy>::
vertex_created(Grid* grid, VertexBase* vrt, GeometricObject* pParent)
{
//	this method is only called if TElem == VertexBase...
	elem_created(grid, reinterpret_cast<TElem*>(vrt), pParent);
}

template <class TElem, class SelectionPolicy>
void
GenericElementSelector<TElem, SelectionPolicy>::
vertex_to_be_erased(Grid* grid, VertexBase* vrt)
{
//	this method is only called if TElem == VertexBase...
	elem_to_be_erased(grid, reinterpret_cast<TElem*>(vrt));
}

//	edge callbacks
template <class TElem, class SelectionPolicy>
void
GenericElementSelector<TElem, SelectionPolicy>::
edge_created(Grid* grid, EdgeBase* edge, GeometricObject* pParent)
{
//	this method is only called if TElem == EdgeBase...
	elem_created(grid,reinterpret_cast<TElem*>(edge), pParent);
}

template <class TElem, class SelectionPolicy>
void
GenericElementSelector<TElem, SelectionPolicy>::
edge_to_be_erased(Grid* grid, EdgeBase* edge)
{
//	this method is only called if TElem == EdgeBase...
	elem_to_be_erased(grid,reinterpret_cast<TElem*>(edge));
}

//	face callbacks
template <class TElem, class SelectionPolicy>
void
GenericElementSelector<TElem, SelectionPolicy>::
face_created(Grid* grid, Face* face, GeometricObject* pParent)
{
//	this method is only called if TElem == Face...
	elem_created(grid,reinterpret_cast<TElem*>(face), pParent);
}

template <class TElem, class SelectionPolicy>
void
GenericElementSelector<TElem, SelectionPolicy>::
face_to_be_erased(Grid* grid, Face* face)
{
//	this method is only called if TElem == Face...
	elem_to_be_erased(grid,reinterpret_cast<TElem*>(face));
}

//	volume callbacks
template <class TElem, class SelectionPolicy>
void
GenericElementSelector<TElem, SelectionPolicy>::
volume_created(Grid* grid, Volume* vol, GeometricObject* pParent)
{
//	this method is only called if TElem == Volume...
	elem_created(grid,reinterpret_cast<TElem*>(vol), pParent);
}

template <class TElem, class SelectionPolicy>
void
GenericElementSelector<TElem, SelectionPolicy>::
volume_to_be_erased(Grid* grid, Volume* vol)
{
//	this method is only called if TElem == Volume...
	elem_to_be_erased(grid,reinterpret_cast<TElem*>(vol));
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	explicit instantiation
template class GenericElementSelector<VertexBase, GridSelectionPolicy<VertexBase> >;
template class GenericElementSelector<EdgeBase, GridSelectionPolicy<EdgeBase> >;
template class GenericElementSelector<Face, GridSelectionPolicy<Face> >;
template class GenericElementSelector<Volume, GridSelectionPolicy<Volume> >;

template class GenericElementSelector<VertexBase, MultiGridSelectionPolicy<VertexBase> >;
template class GenericElementSelector<EdgeBase, MultiGridSelectionPolicy<EdgeBase> >;
template class GenericElementSelector<Face, MultiGridSelectionPolicy<Face> >;
template class GenericElementSelector<Volume, MultiGridSelectionPolicy<Volume> >;


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	implementation of selector
template <class TElementSelectors>
GenericSelector<TElementSelectors>::
GenericSelector()
{
	m_bAutoselectionEnabled = false;
	m_pGrid = NULL;
}

template <class TElementSelectors>
GenericSelector<TElementSelectors>::
GenericSelector(TGridRef grid)
{
	m_bAutoselectionEnabled = false;
	m_pGrid = NULL;
	assign_grid(grid);
}

template <class TElementSelectors>
GenericSelector<TElementSelectors>::
GenericSelector(const ClassType& sel)
{
	assert(!"WARNING in GenericSelector::GenericSelector(const GenericSelector& sel): Copy-Constructor not yet implemented!");
	LOG("WARNING in GenericSelector::GenericSelector(const GenericSelector& sel): Copy-Constructor not yet implemented! Expect unexpected behavior." << endl);
}

template <class TElementSelectors>
GenericSelector<TElementSelectors>::
~GenericSelector()
{
	if(m_pGrid)
	{
	//	unregister observer
		m_pGrid->unregister_observer(this);
	}
}

template <class TElementSelectors>
void
GenericSelector<TElementSelectors>::
assign_grid(TGridRef grid)
{
	grid.register_observer(this, OT_GRID_OBSERVER);
}

template <class TElementSelectors>
typename GenericSelector<TElementSelectors>::TGridPtr
GenericSelector<TElementSelectors>::
get_assigned_grid()
{
	return m_pGrid;
}

template <class TElementSelectors>
void
GenericSelector<TElementSelectors>::
enable_autoselection(bool bEnable)
{
	m_bAutoselectionEnabled = bEnable;
	m_vertexSelector.enable_autoselection(bEnable);
	m_edgeSelector.enable_autoselection(bEnable);
	m_faceSelector.enable_autoselection(bEnable);
	m_volumeSelector.enable_autoselection(bEnable);
}

template <class TElementSelectors>
void
GenericSelector<TElementSelectors>::
enable_selection_inheritance(bool bEnable)
{
	m_bSelectionInheritanceEnabled = bEnable;
	m_vertexSelector.enable_selection_inheritance(bEnable);
	m_edgeSelector.enable_selection_inheritance(bEnable);
	m_faceSelector.enable_selection_inheritance(bEnable);
	m_volumeSelector.enable_selection_inheritance(bEnable);
}

template <class TElementSelectors>
void
GenericSelector<TElementSelectors>::
select(GeometricObject* obj)
{
	int type = obj->base_object_type_id();
	switch(type)
	{
		case VERTEX:	m_vertexSelector.select(static_cast<VertexBase*>(obj)); return;
		case EDGE:		m_edgeSelector.select(static_cast<EdgeBase*>(obj)); return;
		case FACE:		m_faceSelector.select(static_cast<Face*>(obj)); return;
		case VOLUME:	m_volumeSelector.select(static_cast<Volume*>(obj)); return;
	}
}

template <class TElementSelectors>
void
GenericSelector<TElementSelectors>::
deselect(GeometricObject* obj)
{
	int type = obj->base_object_type_id();
	switch(type)
	{
		case VERTEX:	m_vertexSelector.deselect(static_cast<VertexBase*>(obj)); return;
		case EDGE:		m_edgeSelector.deselect(static_cast<EdgeBase*>(obj)); return;
		case FACE:		m_faceSelector.deselect(static_cast<Face*>(obj)); return;
		case VOLUME:	m_volumeSelector.deselect(static_cast<Volume*>(obj)); return;
	}
}

template <class TElementSelectors>
bool
GenericSelector<TElementSelectors>::
is_selected(GeometricObject* obj)
{
	int type = obj->base_object_type_id();
	switch(type)
	{
		case VERTEX:	return m_vertexSelector.is_selected(static_cast<VertexBase*>(obj));
		case EDGE:		return m_edgeSelector.is_selected(static_cast<EdgeBase*>(obj));
		case FACE:		return m_faceSelector.is_selected(static_cast<Face*>(obj));
		case VOLUME:	return m_volumeSelector.is_selected(static_cast<Volume*>(obj));
	}
	return false;
}



//	grid callbacks
template <class TElementSelectors>
void
GenericSelector<TElementSelectors>::
registered_at_grid(Grid* grid)
{
	if(m_pGrid != NULL)
		m_pGrid->unregister_observer(this);

	m_pGrid = dynamic_cast<TGridPtr>(grid);
	if(m_pGrid)
	{
		m_vertexSelector.assign_grid(*m_pGrid);
		m_edgeSelector.assign_grid(*m_pGrid);
		m_faceSelector.assign_grid(*m_pGrid);
		m_volumeSelector.assign_grid(*m_pGrid);
	}
}

template <class TElementSelectors>
void
GenericSelector<TElementSelectors>::
unregistered_from_grid(Grid* grid)
{
	if(m_pGrid)
	{
		assert(dynamic_cast<TGridPtr>(grid) && "grid can not be casted to TGridPtr");
		assert((dynamic_cast<TGridPtr>(grid) == m_pGrid) && "ERROR in Selector::unregistered from grid(...): grids do not match!");
		m_pGrid->unregister_observer(&m_vertexSelector);
		m_pGrid->unregister_observer(&m_edgeSelector);
		m_pGrid->unregister_observer(&m_faceSelector);
		m_pGrid->unregister_observer(&m_volumeSelector);
		m_pGrid = NULL;
	}
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	explicit instantiation
template class GenericSelector<ElementSelectors>;
template class GenericSelector<MGElementSelectors>;



////////////////////////////////////////////////////////////////////////
//	implementation of Selector
//	geometric-object-collection
GeometricObjectCollection
Selector::
get_geometric_object_collection()
{
//TODO: ugly casts! GenericElementSelector should store its selected elements
//		in a GeometricObjectSectionContainer!
	return GeometricObjectCollection(
			(GeometricObjectSectionContainer*)&m_vertexSelector.m_selectedElements,
			(GeometricObjectSectionContainer*)&m_edgeSelector.m_selectedElements,
			(GeometricObjectSectionContainer*)&m_faceSelector.m_selectedElements,
			(GeometricObjectSectionContainer*)&m_volumeSelector.m_selectedElements);
}

////////////////////////////////////////////////////////////////////////
//	implementation of MGSelector
//	geometric-object-collection
MultiLevelGeometricObjectCollection
MGSelector::
get_multi_level_geometric_object_collection()
{
//TODO: ugly casts! GenericElementSelector should store its selected elements
//		in a GeometricObjectSectionContainer!
	uint numLevels = num_levels();
	MultiLevelGeometricObjectCollection mgoc(numLevels);
	
	for(uint i = 0; i < numLevels; ++i)
	{
		mgoc.add_level(
			(GeometricObjectSectionContainer*)&m_vertexSelector.get_section(i),
			(GeometricObjectSectionContainer*)&m_edgeSelector.get_section(i),
			(GeometricObjectSectionContainer*)&m_faceSelector.get_section(i),
			(GeometricObjectSectionContainer*)&m_volumeSelector.get_section(i));
	}
	
	return mgoc;
}

}//	end of namespace
