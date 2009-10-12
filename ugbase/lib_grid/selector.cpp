// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y08 m11 d20

#include <cassert>
#include "selector.h"
#include "common/common.h"

using namespace std;

namespace ug
{
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	implementation of GenericSelector
template <class TElem>
GenericSelector<TElem>::
GenericSelector() : m_baseObjectType(geometry_traits<TElem>::BASE_OBJECT_TYPE_ID),
					m_aElemIterator(false)
{
	m_bAutoselectionEnabled = false;
	m_pGrid = NULL;
	m_invalidContainer.push_back(NULL);
}

template <class TElem>
GenericSelector<TElem>::
GenericSelector(Grid& grid) : m_baseObjectType(geometry_traits<TElem>::BASE_OBJECT_TYPE_ID)
{
	m_bAutoselectionEnabled = false;
	m_pGrid = NULL;
	m_invalidContainer.push_back(NULL);
	assign_grid(grid);
}

template <class TElem>
GenericSelector<TElem>::
GenericSelector(const GenericSelector<TElem>& gSel) : m_baseObjectType(geometry_traits<TElem>::BASE_OBJECT_TYPE_ID)
{
	assert(!"WARNING in GenericSelector::GenericSelector(const GenericSelector& sel): Copy-Constructor not yet implemented!");
	LOG("WARNING in GenericSelector::GenericSelector(const GenericSelector& sel): Copy-Constructor not yet implemented! Expect unexpected behaviour!");
}

template <class TElem>
GenericSelector<TElem>::~GenericSelector()
{
	if(m_pGrid)
	{
	//	unregister the observer.
		m_pGrid->unregister_observer(this);
	}
}

template <class TElem>
void
GenericSelector<TElem>::
assign_grid(Grid& grid)
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

template <class TElem>
Grid*
GenericSelector<TElem>::
get_assigned_grid()
{
	return m_pGrid;
}


template <class TElem>
void
GenericSelector<TElem>::
enable_autoselection(bool bEnable)
{
	m_bAutoselectionEnabled = bEnable;
}

template <class TElem>
void
GenericSelector<TElem>::
enable_selection_inheritance(bool bEnable)
{
	m_bSelectionInheritanceEnabled = bEnable;
}

template <class TElem>
void
GenericSelector<TElem>::
select(TElem* elem)
{
	assert(m_pGrid && "ERROR in GenericSelector::select(...): selector not registered at any grid!");

	if(!is_selected(elem))
	{
		m_aaElemIterator[elem] = m_selectedElements.insert(elem, elem->shared_pipe_section());
	}
}

template <class TElem>
void
GenericSelector<TElem>::
deselect(TElem* elem)
{
	assert(m_pGrid && "ERROR in GenericSelector::deselect(...): selector not registered at any grid!");
	if(is_selected(elem))
	{
		m_selectedElements.erase(m_aaElemIterator[elem], elem->shared_pipe_section());
		m_aaElemIterator[elem] = m_invalidContainer.begin();
	}
}

template <class TElem>
bool
GenericSelector<TElem>::
is_selected(TElem* elem)
{
	assert(m_pGrid && "ERROR in GenericSelector::is_selected(...): selector not registered at any grid!");
	if(m_aaElemIterator[elem] == (m_invalidContainer.begin()))
		return false;
	return true;
}

template <class TElem>
bool
GenericSelector<TElem>::
is_selected(GeometricObject* elem)
{
	assert(m_pGrid && "ERROR in GenericSelector::is_selected(...): selector not registered at any grid!");
	TElem* pElem = dynamic_cast<TElem*>(elem);
	if(pElem != NULL)
	{
		if(m_aaElemIterator[elem] == (m_invalidContainer.begin()))
			return false;
		return true;
	}
	return false;
}
		
template <class TElem>
void
GenericSelector<TElem>::
elem_created(Grid* grid, TElem* elem, GeometricObject* pParent)
{
	assert((m_pGrid == grid) && "ERROR in GenericSelector::elem_created(...): selectors grid and elements grid do not match!");
//	init deselected.
	m_aaElemIterator[elem] = m_invalidContainer.begin();

	if(autoselection_enabled() == true)
		select(elem);
	else if(selection_inheritance_enabled() == true)
	{
		if(pParent != NULL)
		{
			if(is_selected(pParent))
				select(elem);
		}
	}
}

template <class TElem>
void
GenericSelector<TElem>::
elem_to_be_erased(Grid* grid, TElem* elem)
{
	assert((m_pGrid == grid) && "ERROR in GenericSelector::elem_to_be_erased(...): selectors grid and elements grid do not match!");

	if(is_selected(elem))
		m_selectedElements.erase(m_aaElemIterator[elem], elem->shared_pipe_section());

}

//	grid callbacks
template <class TElem>
void
GenericSelector<TElem>::
registered_at_grid(Grid* grid)
{
	if(m_pGrid != NULL)
		m_pGrid->unregister_observer(this);

	m_pGrid = grid;

	grid->attach_to<TElem>(m_aElemIterator, false);
	m_aaElemIterator.access(*grid, m_aElemIterator);

//	initialize all iterators attached to base objects with the invalid iterator
	for(TElemIterator iter = grid->begin<TElem>();
		iter != grid->end<TElem>(); iter++)
	{
		m_aaElemIterator[*iter] = m_invalidContainer.begin();
	}

//	clear the list
	m_selectedElements.clear();
}

template <class TElem>
void
GenericSelector<TElem>::
unregistered_from_grid(Grid* grid)
{
	assert((grid == m_pGrid) && "ERROR in GenericSelector::unregistered_from_grid(...): grids do not match!");
	grid->detach_from<TElem>(m_aElemIterator);
	m_selectedElements.clear();
	m_pGrid = NULL;
}

template <class TElem>
void
GenericSelector<TElem>::
elements_to_be_cleared(Grid* grid)
{
	assert((grid == m_pGrid) && "ERROR in GenericSelector::elements_to_be_cleared(...): grids do not match!");
	clear_selection();
}

//	vertex callbacks
template <class TElem>
void
GenericSelector<TElem>::
vertex_created(Grid* grid, VertexBase* vrt, GeometricObject* pParent)
{
//	this method is only called if TElem == VertexBase...
	elem_created(grid, reinterpret_cast<TElem*>(vrt), pParent);
}

template <class TElem>
void
GenericSelector<TElem>::
vertex_to_be_erased(Grid* grid, VertexBase* vrt)
{
//	this method is only called if TElem == VertexBase...
	elem_to_be_erased(grid, reinterpret_cast<TElem*>(vrt));
}

//	edge callbacks
template <class TElem>
void
GenericSelector<TElem>::
edge_created(Grid* grid, EdgeBase* edge, GeometricObject* pParent)
{
//	this method is only called if TElem == EdgeBase...
	elem_created(grid,reinterpret_cast<TElem*>(edge), pParent);
}

template <class TElem>
void
GenericSelector<TElem>::
edge_to_be_erased(Grid* grid, EdgeBase* edge)
{
//	this method is only called if TElem == EdgeBase...
	elem_to_be_erased(grid,reinterpret_cast<TElem*>(edge));
}

//	face callbacks
template <class TElem>
void
GenericSelector<TElem>::
face_created(Grid* grid, Face* face, GeometricObject* pParent)
{
//	this method is only called if TElem == Face...
	elem_created(grid,reinterpret_cast<TElem*>(face), pParent);
}

template <class TElem>
void
GenericSelector<TElem>::
face_to_be_erased(Grid* grid, Face* face)
{
//	this method is only called if TElem == Face...
	elem_to_be_erased(grid,reinterpret_cast<TElem*>(face));
}

//	volume callbacks
template <class TElem>
void
GenericSelector<TElem>::
volume_created(Grid* grid, Volume* vol, GeometricObject* pParent)
{
//	this method is only called if TElem == Volume...
	elem_created(grid,reinterpret_cast<TElem*>(vol), pParent);
}

template <class TElem>
void
GenericSelector<TElem>::
volume_to_be_erased(Grid* grid, Volume* vol)
{
//	this method is only called if TElem == Volume...
	elem_to_be_erased(grid,reinterpret_cast<TElem*>(vol));
}

//	explicit instantiation
template class GenericSelector<VertexBase>;
template class GenericSelector<EdgeBase>;
template class GenericSelector<Face>;
template class GenericSelector<Volume>;


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	implementation of selector
Selector::
Selector()
{
	m_bAutoselectionEnabled = false;
	m_pGrid = NULL;
}

Selector::
Selector(Grid& grid)
{
	m_bAutoselectionEnabled = false;
	m_pGrid = NULL;
	assign_grid(grid);
}

Selector::
Selector(const Selector& sel)
{
	assert(!"WARNING in Selector::Selector(const Selector& sel): Copy-Constructor not yet implemented!");
	LOG("WARNING in Selector::Selector(const Selector& sel): Copy-Constructor not yet implemented! Expect unexpected behavior." << endl);
}

Selector::
~Selector()
{
	if(m_pGrid)
	{
	//	unregister observer
		m_pGrid->unregister_observer(this);
	}
}

void
Selector::
assign_grid(Grid& grid)
{
	grid.register_observer(this, OT_GRID_OBSERVER);
}

Grid*
Selector::
get_assigned_grid()
{
	return m_pGrid;
}

void
Selector::
enable_autoselection(bool bEnable)
{
	m_bAutoselectionEnabled = bEnable;
	m_vertexSelector.enable_autoselection(bEnable);
	m_edgeSelector.enable_autoselection(bEnable);
	m_faceSelector.enable_autoselection(bEnable);
	m_volumeSelector.enable_autoselection(bEnable);
}

void
Selector::
enable_selection_inheritance(bool bEnable)
{
	m_bSelectionInheritanceEnabled = bEnable;
	m_vertexSelector.enable_selection_inheritance(bEnable);
	m_edgeSelector.enable_selection_inheritance(bEnable);
	m_faceSelector.enable_selection_inheritance(bEnable);
	m_volumeSelector.enable_selection_inheritance(bEnable);
}

//	geometric-object-collection
GeometricObjectCollection
Selector::
get_geometric_object_collection()
{
//TODO: ugly casts! GenericSelector should store its selected elements
//		in a GeometricObjectSectionContainer!
	return GeometricObjectCollection(
			(GeometricObjectSectionContainer*)&m_vertexSelector.m_selectedElements,
			(GeometricObjectSectionContainer*)&m_edgeSelector.m_selectedElements,
			(GeometricObjectSectionContainer*)&m_faceSelector.m_selectedElements,
			(GeometricObjectSectionContainer*)&m_volumeSelector.m_selectedElements);
}

//	grid callbacks
void
Selector::
registered_at_grid(Grid* grid)
{
	if(m_pGrid != NULL)
		m_pGrid->unregister_observer(this);

	m_pGrid = grid;
	m_vertexSelector.assign_grid(*m_pGrid);
	m_edgeSelector.assign_grid(*m_pGrid);
	m_faceSelector.assign_grid(*m_pGrid);
	m_volumeSelector.assign_grid(*m_pGrid);
}

void
Selector::
unregistered_from_grid(Grid* grid)
{
	assert((grid == m_pGrid) && "ERROR in Selector::unregistered from grid(...): grids do not match!");
	m_pGrid->unregister_observer(&m_vertexSelector);
	m_pGrid->unregister_observer(&m_edgeSelector);
	m_pGrid->unregister_observer(&m_faceSelector);
	m_pGrid->unregister_observer(&m_volumeSelector);
	m_pGrid = NULL;
}

}//	end of namespace
