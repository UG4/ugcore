// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y10 m02 d16

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	...
////////////////////////////////////////////////////////////////////////

#ifndef __H__LIBGRID__SELECTOR_GRID_ELEM_IMPL__
#define __H__LIBGRID__SELECTOR_GRID_ELEM_IMPL__

#include <cassert>

namespace ug
{

template <class BaseElem>
template <class TElem>
inline ISelector::SectionContainer&
TElemSelector<BaseElem>::get_section_container()
{
	assert(((int)geometry_traits<TElem>::BASE_OBJECT_TYPE_ID
			== (int)geometry_traits<BaseElem>::BASE_OBJECT_TYPE_ID)
			&& "Element type does not match BaseElem type");

	return m_elements;
}

template <class BaseElem>
template <class TElem>
inline int
TElemSelector<BaseElem>::get_section_index()
{
	return geometry_traits<TElem>::SHARED_PIPE_SECTION;
}

template <class BaseElem>
template <class TElem>
inline void
TElemSelector<BaseElem>::clear()
{
	if(m_pGrid){
	//	mark all elements as deselected
		typename geometry_traits<TElem>::iterator iter;
		for(iter = begin<TElem>(); iter != end<TElem>(); ++iter)
			mark_deselected(*iter);

	//	clear the section
		const int sInd = get_section_index<TElem>();
		if(sInd < 0)
			get_section_container<TElem>().clear();
		else
			get_section_container<TElem>().clear_section(sInd);
	}
}

template <class BaseElem>
template <class TElem>
inline uint 
TElemSelector<BaseElem>::num()
{
	const int sInd = get_section_index<TElem>();
	if(sInd < 0)
		return get_section_container<TElem>().num_elements();
	else
		return get_section_container<TElem>().num_elements(sInd);
}

template <class BaseElem>
inline uint 
TElemSelector<BaseElem>::num()
{
	return num<BaseElem>();
}

//	empty
template <class BaseElem>
inline bool 
TElemSelector<BaseElem>::empty()
{
	return num() == 0;
}

template <class BaseElem>
template <class TElem>
inline bool 
TElemSelector<BaseElem>::empty()
{
	return num<TElem>() == 0;
}

//	begin
template <class BaseElem>
template <class TElem>
inline typename geometry_traits<TElem>::iterator
TElemSelector<BaseElem>::begin()
{
	const int sInd = get_section_index<TElem>();
	if(sInd < 0)
		return iterator_cast<typename geometry_traits<TElem>::iterator>(
								get_section_container<TElem>().begin());
	else
		return iterator_cast<typename geometry_traits<TElem>::iterator>(
					get_section_container<TElem>().section_begin(sInd));
}

template <class BaseElem>
inline typename TElemSelector<BaseElem>::BaseElemIterator
TElemSelector<BaseElem>::begin()
{
	return iterator_cast<BaseElemIterator>(m_elements.begin());
}

//	end
template <class BaseElem>
template <class TElem>
inline typename geometry_traits<TElem>::iterator
TElemSelector<BaseElem>::end()
{
	const int sInd = get_section_index<TElem>();
	if(sInd < 0)
		return iterator_cast<typename geometry_traits<TElem>::iterator>(
									get_section_container<TElem>().end());
	else
		return iterator_cast<typename geometry_traits<TElem>::iterator>(
								get_section_container<TElem>().section_end(sInd));
}

template <class BaseElem>
inline typename TElemSelector<BaseElem>::BaseElemIterator
TElemSelector<BaseElem>::end()
{
	return iterator_cast<BaseElemIterator>(m_elements.end());
}

}//	end of namespace

#endif
