// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y09 m11 d09

#ifndef __H__LIB_GRID__SELECTION_POLICIES_IMPL__
#define __H__LIB_GRID__SELECTION_POLICIES_IMPL__

namespace ug
{

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	implementation of GridSelectionPolicy
template <class TElem>
template <class TSelElem>
void
GridSelectionPolicy<TElem>::
clear_selection()
{
	if(geometry_traits<TSelElem>::BASE_OBJECT_TYPE_ID == m_baseObjectType)
	{
		int pipeSection = geometry_traits<TSelElem>::SHARED_PIPE_SECTION;
		for(TElemIterator iter = m_selectedElements.section_begin(pipeSection);
				iter != static_cast<TElemIterator>(m_selectedElements.section_end(pipeSection)); iter++)
			m_aaElemIterator[*iter] = m_invalidContainer.begin();

		if(pipeSection == -1)
			m_selectedElements.clear();
		else
			m_selectedElements.clear_section(pipeSection);
	}
}

template <class TElem>
template <class TSelElem>
uint
GridSelectionPolicy<TElem>::
num_selected()
{
	if(geometry_traits<TSelElem>::BASE_OBJECT_TYPE_ID == m_baseObjectType)
	{
		int secInd = geometry_traits<TSelElem>::SHARED_PIPE_SECTION;
		if(secInd == -1)
			return m_selectedElements.num_elements();
		else
			return m_selectedElements.num_elements(secInd);
	}
	return 0;
}

template <class TElem>
template <class TSelElem>
typename geometry_traits<TSelElem>::iterator
GridSelectionPolicy<TElem>::
begin()
{
	if(geometry_traits<TSelElem>::BASE_OBJECT_TYPE_ID == m_baseObjectType)
		return iterator_cast<typename geometry_traits<TSelElem>::iterator >(
				m_selectedElements.section_begin(geometry_traits<TSelElem>::SHARED_PIPE_SECTION));

	return iterator_cast<typename geometry_traits<TSelElem>::iterator>(
			m_selectedElements.end());
}

template <class TElem>
template <class TSelElem>
typename geometry_traits<TSelElem>::iterator
GridSelectionPolicy<TElem>::
end()
{
	if(geometry_traits<TSelElem>::BASE_OBJECT_TYPE_ID == m_baseObjectType)
		return iterator_cast<typename geometry_traits<TSelElem>::iterator >(
				m_selectedElements.section_end(geometry_traits<TSelElem>::SHARED_PIPE_SECTION));

	return iterator_cast<typename geometry_traits<TSelElem>::iterator>(
			m_selectedElements.end());
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	implementation of MultiGridSelectionPolicy
template <class TElem>
template <class TSelElem>
void
MultiGridSelectionPolicy<TElem>::
clear_selection()
{
//	iterate over all layers and clear the elements of the given type in each.
	if(geometry_traits<TSelElem>::BASE_OBJECT_TYPE_ID == m_baseObjectType)
	{
		int pipeSection = geometry_traits<TSelElem>::SHARED_PIPE_SECTION;

		for(size_t i = 0; i < m_vSections.size(); ++i)
		{
			ElemSectionContainer& selectedElements = *m_vSections[i];
			for(TElemIterator iter = selectedElements.section_begin(pipeSection);
					iter != static_cast<TElemIterator>(selectedElements.section_end(pipeSection)); iter++)
				m_aaElemIterator[*iter] = m_invalidContainer.begin();

			if(pipeSection == -1)
				selectedElements.clear();
			else
				selectedElements.clear_section(pipeSection);
		}
	}
}

template <class TElem>
template <class TSelElem>
void
MultiGridSelectionPolicy<TElem>::
clear_selection(int level)
{
//	iterate over all layers and clear the elements of the given type in each.
	if(geometry_traits<TSelElem>::BASE_OBJECT_TYPE_ID == m_baseObjectType)
	{
		int pipeSection = geometry_traits<TSelElem>::SHARED_PIPE_SECTION;

		ElemSectionContainer& selectedElements = get_section(level);
		for(TElemIterator iter = selectedElements.section_begin(pipeSection);
				iter != static_cast<TElemIterator>(selectedElements.section_end(pipeSection)); iter++)
			m_aaElemIterator[*iter] = m_invalidContainer.begin();

		if(pipeSection == -1)
			selectedElements.clear();
		else
			selectedElements.clear_section(pipeSection);
	}
}

template <class TElem>
template <class TSelElem>
uint
MultiGridSelectionPolicy<TElem>::
num_selected()
{
	int numSel = 0;
//	we have to sum the number of selected elements in each layer.
	if(geometry_traits<TSelElem>::BASE_OBJECT_TYPE_ID == m_baseObjectType)
	{
		int secInd = geometry_traits<TSelElem>::SHARED_PIPE_SECTION;
		if(secInd == -1)
		{
			for(size_t i = 0; i < m_vSections.size(); ++i)
				numSel += m_vSections[i]->num_elements();
		}
		else
		{
			for(size_t i = 0; i < m_vSections.size(); ++i)
				numSel += m_vSections[i]->num_elements(secInd);
		}
	}
	return numSel;
}

template <class TElem>
template <class TSelElem>
uint
MultiGridSelectionPolicy<TElem>::
num_selected(int level)
{
	if(geometry_traits<TSelElem>::BASE_OBJECT_TYPE_ID == m_baseObjectType)
	{
		int secInd = geometry_traits<TSelElem>::SHARED_PIPE_SECTION;
		if(secInd == -1)
			return get_section(level).num_elements();
		else
			return get_section(level).num_elements(secInd);
	}
	return 0;
}

template <class TElem>
template <class TSelElem>
typename geometry_traits<TSelElem>::iterator
MultiGridSelectionPolicy<TElem>::
begin(int level)
{
	if(geometry_traits<TSelElem>::BASE_OBJECT_TYPE_ID == m_baseObjectType)
		return iterator_cast<typename geometry_traits<TSelElem>::iterator >(
				get_section(level).section_begin(geometry_traits<TSelElem>::SHARED_PIPE_SECTION));

	return iterator_cast<typename geometry_traits<TSelElem>::iterator>(
			get_section(level).end());
}

template <class TElem>
template <class TSelElem>
typename geometry_traits<TSelElem>::iterator
MultiGridSelectionPolicy<TElem>::
end(int level)
{
	if(geometry_traits<TSelElem>::BASE_OBJECT_TYPE_ID == m_baseObjectType)
		return iterator_cast<typename geometry_traits<TSelElem>::iterator >(
				get_section(level).section_end(geometry_traits<TSelElem>::SHARED_PIPE_SECTION));

	return iterator_cast<typename geometry_traits<TSelElem>::iterator>(
			get_section(level).end());
}

}//	end of namespace

#endif
