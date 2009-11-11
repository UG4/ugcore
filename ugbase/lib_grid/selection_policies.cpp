// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y09 m11 d09

#include "selection_policies.h"
#include "grid/grid.h"
#include "multi_grid.h"

using namespace std;

namespace ug
{
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	implementation of GridSelectionPolicy
template <class TElem>
GridSelectionPolicy<TElem>::
GridSelectionPolicy() : m_aElemIterator(false), 
						m_baseObjectType(geometry_traits<TElem>::BASE_OBJECT_TYPE_ID)
{
	m_pGrid = NULL;
	m_invalidContainer.push_back(NULL);
}

template <class TElem>
GridSelectionPolicy<TElem>::
~GridSelectionPolicy()
{
//	clean up
	new_grid(NULL);
}

template <class TElem>
void
GridSelectionPolicy<TElem>::
new_grid(Grid* gridNew)
{
	if(m_pGrid)
	{
	//	unregister attachment.
		m_pGrid->detach_from<TElem>(m_aElemIterator);
		m_selectedElements.clear();
		m_pGrid = NULL;
	}
	
	m_pGrid = gridNew;
	if(m_pGrid)
	{
		m_pGrid->attach_to<TElem>(m_aElemIterator, false);
		m_aaElemIterator.access(*m_pGrid, m_aElemIterator);

	//	initialize all iterators attached to base objects with the invalid iterator
		for(TElemIterator iter = m_pGrid->begin<TElem>();
			iter != m_pGrid->end<TElem>(); iter++)
		{
			m_aaElemIterator[*iter] = m_invalidContainer.begin();
		}

	//	clear the list
		m_selectedElements.clear();
	}
}

template <class TElem>
void
GridSelectionPolicy<TElem>::
select(TElem* elem)
{
	assert(m_pGrid && "ERROR in GridSelectionPolicy::select(...): selector not registered at any grid!");

	if(!is_selected(elem))
	{
		m_aaElemIterator[elem] = m_selectedElements.insert(elem, elem->shared_pipe_section());
	}
}

template <class TElem>
void
GridSelectionPolicy<TElem>::
deselect(TElem* elem)
{
	assert(m_pGrid && "ERROR in GridSelectionPolicy::deselect(...): selector not registered at any grid!");
	if(is_selected(elem))
	{
		m_selectedElements.erase(m_aaElemIterator[elem], elem->shared_pipe_section());
		m_aaElemIterator[elem] = m_invalidContainer.begin();
	}
}

template <class TElem>
bool
GridSelectionPolicy<TElem>::
is_selected(TElem* elem)
{
	assert(m_pGrid && "ERROR in GridSelectionPolicy::is_selected(...): selector not registered at any grid!");
	if(m_aaElemIterator[elem] == (m_invalidContainer.begin()))
		return false;
	return true;
}

template <class TElem>
bool
GridSelectionPolicy<TElem>::
is_selected(GeometricObject* elem)
{
	assert(m_pGrid && "ERROR in GridSelectionPolicy::is_selected(...): selector not registered at any grid!");
	TElem* pElem = dynamic_cast<TElem*>(elem);
	if(pElem != NULL)
	{
		if(m_aaElemIterator[elem] == (m_invalidContainer.begin()))
			return false;
		return true;
	}
	return false;
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

//	explicit instantiation
template class GridSelectionPolicy<VertexBase>;
template class GridSelectionPolicy<EdgeBase>;
template class GridSelectionPolicy<Face>;
template class GridSelectionPolicy<Volume>;


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	implementation of MultiGridSelectionPolicy
template <class TElem>
MultiGridSelectionPolicy<TElem>::
MultiGridSelectionPolicy() : m_aElemIterator(false), 
							m_baseObjectType(geometry_traits<TElem>::BASE_OBJECT_TYPE_ID)
{
	m_pMultiGrid = NULL;
	m_invalidContainer.push_back(NULL);
}

template <class TElem>
MultiGridSelectionPolicy<TElem>::
~MultiGridSelectionPolicy()
{
//	clean up
	if(m_pMultiGrid)
	{
	//	unregister attachment.
		m_pMultiGrid->detach_from<TElem>(m_aElemIterator);
	//	delete section containers
		for(size_t i = 0; i < m_vSections.size(); ++i)
			delete m_vSections[i];
	}
}

template <class TElem>
void
MultiGridSelectionPolicy<TElem>::
new_grid(MultiGrid* gridNew)
{
	if(m_pMultiGrid)
	{
	//	unregister attachment.
		m_pMultiGrid->detach_from<TElem>(m_aElemIterator);
	//	clear sections
		for(size_t i = 0; i < m_vSections.size(); ++i)
			m_vSections[i]->clear();

		m_pMultiGrid = NULL;
	}
	
	m_pMultiGrid = gridNew;
	if(m_pMultiGrid)
	{
		m_pMultiGrid->attach_to<TElem>(m_aElemIterator, false);
		m_aaElemIterator.access(*m_pMultiGrid, m_aElemIterator);

	//	initialize all iterators attached to base objects with the invalid iterator
		for(TElemIterator iter = m_pMultiGrid->begin<TElem>();
			iter != m_pMultiGrid->end<TElem>(); iter++)
		{
			m_aaElemIterator[*iter] = m_invalidContainer.begin();
		}
	}
}

template <class TElem>
void
MultiGridSelectionPolicy<TElem>::
select(TElem* elem)
{
	assert(m_pMultiGrid && "ERROR in MultiGridSelectionPolicy::select(...): selector not registered at any grid!");

	if(!is_selected(elem))
	{
	//	get the level of the element
		int level = m_pMultiGrid->get_level(elem);
		m_aaElemIterator[elem] = get_section(level).insert(elem, elem->shared_pipe_section());
	}
}

template <class TElem>
void
MultiGridSelectionPolicy<TElem>::
deselect(TElem* elem)
{
	assert(m_pMultiGrid && "ERROR in MultiGridSelectionPolicy::deselect(...): selector not registered at any grid!");
	if(is_selected(elem))
	{
		int level = m_pMultiGrid->get_level(elem);
		get_section(level).erase(m_aaElemIterator[elem], elem->shared_pipe_section());
		m_aaElemIterator[elem] = m_invalidContainer.begin();
	}
}

template <class TElem>
bool
MultiGridSelectionPolicy<TElem>::
is_selected(TElem* elem)
{
	assert(m_pMultiGrid && "ERROR in MultiGridSelectionPolicy::is_selected(...): selector not registered at any grid!");
	if(m_aaElemIterator[elem] == (m_invalidContainer.begin()))
		return false;
	return true;
}

template <class TElem>
bool
MultiGridSelectionPolicy<TElem>::
is_selected(GeometricObject* elem)
{
	assert(m_pMultiGrid && "ERROR in MultiGridSelectionPolicy::is_selected(...): selector not registered at any grid!");
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
MultiGridSelectionPolicy<TElem>::
grow_sections(int newSize)
{
//	create new SectionContainers and push them to the list,
//	until there are enough of them.
	while((int)m_vSections.size() < newSize)
		m_vSections.push_back(new ElemSectionContainer);
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

//	explicit instantiation
template class MultiGridSelectionPolicy<VertexBase>;
template class MultiGridSelectionPolicy<EdgeBase>;
template class MultiGridSelectionPolicy<Face>;
template class MultiGridSelectionPolicy<Volume>;

}//	end of namespace
