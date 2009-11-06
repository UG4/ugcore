// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y08 m11 d19

#ifndef __H__LIBGRID__SELECTOR_IMPL__
#define __H__LIBGRID__SELECTOR_IMPL__

#include <cassert>
//#include "selector.h"

namespace ug
{

template <class TElem>
template <class TIterator>
void
GenericElementSelector<TElem>::
select(TIterator iterBegin, TIterator iterEnd)
{
	assert(m_pGrid && "ERROR in GenericElementSelector::select(...): selector not registered at any grid!");
	while(iterBegin != iterEnd)
	{
		select(*iterBegin);
		iterBegin++;
	}
}

template <class TElem>
template <class TIterator>
void
GenericElementSelector<TElem>::
deselect(TIterator iterBegin, TIterator iterEnd)
{
	assert(m_pGrid && "ERROR in GenericElementSelector::deselect(...): selector not registered at any grid!");
	while(iterBegin != iterEnd)
	{
		deselect(*iterBegin);
		iterBegin++;
	}
}

template <class TElem>
template <class TSelElem>
void
GenericElementSelector<TElem>::
select_all()
{
	assert(m_pGrid && "ERROR in GenericElementSelector::select(...): selector not registered at any grid!");

	if(geometry_traits<TSelElem>::BASE_OBJECT_TYPE_ID == m_baseObjectType)
	{
		for(typename geometry_traits<TSelElem>::iterator iter = m_pGrid->begin<TSelElem>();
				iter != m_pGrid->end<TSelElem>(); iter++)
			select(*iter);
	}
}

template <class TElem>
template <class TSelElem>
void
GenericElementSelector<TElem>::
clear_selection()
{
	assert(m_pGrid && "ERROR in GenericElementSelector::clear_selection(...): selector not registered at any grid!");

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
GenericElementSelector<TElem>::
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
GenericElementSelector<TElem>::
begin()
{
	assert(m_pGrid && "ERROR in GenericElementSelector::begin(...): selector not registered at any grid!");

	if(geometry_traits<TSelElem>::BASE_OBJECT_TYPE_ID == m_baseObjectType)
		return iterator_cast<typename geometry_traits<TSelElem>::iterator >(
				m_selectedElements.section_begin(geometry_traits<TSelElem>::SHARED_PIPE_SECTION));

	return iterator_cast<typename geometry_traits<TSelElem>::iterator>(
			m_selectedElements.end());
}

template <class TElem>
template <class TSelElem>
typename geometry_traits<TSelElem>::iterator
GenericElementSelector<TElem>::
end()
{
	assert(m_pGrid && "ERROR in GenericElementSelector::end(...): selector not registered at any grid!");
	if(geometry_traits<TSelElem>::BASE_OBJECT_TYPE_ID == m_baseObjectType)
		return iterator_cast<typename geometry_traits<TSelElem>::iterator >(
				m_selectedElements.section_end(geometry_traits<TSelElem>::SHARED_PIPE_SECTION));

	return iterator_cast<typename geometry_traits<TSelElem>::iterator>(
			m_selectedElements.end());
}



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	implementation of selector
template <class TSelElem>
void
Selector::
select_all()
{
	TSelElem* pTmp(NULL);
	select_all<TSelElem>(pTmp);
}

template <class TSelElem>
void
Selector::
clear_selection()
{
	TSelElem* pTmp(NULL);
	clear_selection<TSelElem>(pTmp);
}

template <class TSelElem>
typename geometry_traits<TSelElem>::iterator
Selector::
begin()
{
	TSelElem* pTmp(NULL);
	return begin<TSelElem>(pTmp);
}

template <class TSelElem>
typename geometry_traits<TSelElem>::iterator
Selector::
end()
{
	TSelElem* pTmp(NULL);
	return end<TSelElem>(pTmp);
}

}//	end of namespace

#endif
