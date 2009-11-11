// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y08 m11 d19

#ifndef __H__LIBGRID__SELECTOR_IMPL__
#define __H__LIBGRID__SELECTOR_IMPL__

#include <cassert>
//#include "selector.h"

namespace ug
{
template <class TElem, class SelectionPolicy>
template <class TIterator>
void
GenericElementSelector<TElem, SelectionPolicy>::
select(TIterator iterBegin, TIterator iterEnd)
{
	assert(m_pGrid && "ERROR in GenericElementSelector::select(...): selector not registered at any grid!");
	while(iterBegin != iterEnd)
	{
		SelectionPolicy::select(*iterBegin);
		iterBegin++;
	}
}

template <class TElem, class SelectionPolicy>
template <class TIterator>
void
GenericElementSelector<TElem, SelectionPolicy>::
deselect(TIterator iterBegin, TIterator iterEnd)
{
	assert(m_pGrid && "ERROR in GenericElementSelector::deselect(...): selector not registered at any grid!");
	while(iterBegin != iterEnd)
	{
		SelectionPolicy::deselect(*iterBegin);
		iterBegin++;
	}
}

template <class TElem, class SelectionPolicy>
template <class TSelElem>
void
GenericElementSelector<TElem, SelectionPolicy>::
select_all()
{
	assert(m_pGrid && "ERROR in GenericElementSelector::select(...): selector not registered at any grid!");

	if(geometry_traits<TSelElem>::BASE_OBJECT_TYPE_ID == geometry_traits<TElem>::BASE_OBJECT_TYPE_ID)
	{
		for(typename geometry_traits<TSelElem>::iterator iter = m_pGrid->template begin<TSelElem>();
				iter != m_pGrid->template end<TSelElem>(); iter++)
			SelectionPolicy::select(*iter);
	}
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	implementation of selector
template <class TElementSelectors>
template <class TSelElem>
void
GenericSelector<TElementSelectors>::
select_all()
{
	TSelElem* pTmp(NULL);
	select_all<TSelElem>(pTmp);
}

template <class TElementSelectors>
template <class TSelElem>
void
GenericSelector<TElementSelectors>::
clear_selection()
{
	TSelElem* pTmp(NULL);
	clear_selection<TSelElem>(pTmp);
}
/*
template <class TElementSelectors>
template <class TSelElem>
typename geometry_traits<TSelElem>::iterator
GenericSelector<TElementSelectors>::
begin()
{
	TSelElem* pTmp(NULL);
	return begin<TSelElem>(pTmp);
}

template <class TElementSelectors>
template <class TSelElem>
typename geometry_traits<TSelElem>::iterator
GenericSelector<TElementSelectors>::
end()
{
	TSelElem* pTmp(NULL);
	return end<TSelElem>(pTmp);
}
*/
}//	end of namespace

#endif
