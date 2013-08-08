// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y10 m02 d15

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	...
////////////////////////////////////////////////////////////////////////

#ifndef __H__LIBGRID__SELECTOR_GRID_IMPL__
#define __H__LIBGRID__SELECTOR_GRID_IMPL__

#include <cassert>

namespace ug
{

template <class TElem>
inline int
Selector::get_section_index() const
{
	return geometry_traits<TElem>::CONTAINER_SECTION;
}

template <class TElem>
inline void
Selector::clear()
{
	if(m_pGrid){
	//	mark all elements as deselected
		typename geometry_traits<TElem>::iterator iter;
		for(iter = begin<TElem>(); iter != end<TElem>(); ++iter)
			mark_deselected(*iter);

	//	clear the section
		const int sInd = get_section_index<TElem>();
		if(sInd < 0)
			section_container<TElem>().clear();
		else
			section_container<TElem>().clear_section(sInd);
	}
}

template <class TElem>
inline size_t
Selector::num() const
{
	const int sInd = get_section_index<TElem>();
	if(sInd < 0)
		return section_container<TElem>().num_elements();
	else
		return section_container<TElem>().num_elements(sInd);
}

inline size_t
Selector::num() const
{
	return num<VertexBase>() + num<EdgeBase>() + num<Face>() + num<Volume>();
}

//	empty
inline bool 
Selector::empty() const
{
	return num() == 0;
}

template <class TElem>
inline bool 
Selector::empty() const
{
	return num<TElem>() == 0;
}

//	begin
template <class TElem>
inline typename geometry_traits<TElem>::iterator
Selector::begin()
{
	const int sInd = get_section_index<TElem>();
	if(sInd < 0)
		return iterator_cast<typename geometry_traits<TElem>::iterator>(
								section_container<TElem>().begin());
	else
		return iterator_cast<typename geometry_traits<TElem>::iterator>(
					section_container<TElem>().section_begin(sInd));
}

//	const begin
template <class TElem>
inline typename geometry_traits<TElem>::const_iterator
Selector::begin() const
{
	const int sInd = get_section_index<TElem>();
	if(sInd < 0)
		return iterator_cast<typename geometry_traits<TElem>::const_iterator>(
								section_container<TElem>().begin());
	else
		return iterator_cast<typename geometry_traits<TElem>::const_iterator>(
					section_container<TElem>().section_begin(sInd));
}

//	end
template <class TElem>
inline typename geometry_traits<TElem>::iterator
Selector::end()
{
	const int sInd = get_section_index<TElem>();
	if(sInd < 0)
		return iterator_cast<typename geometry_traits<TElem>::iterator>(
									section_container<TElem>().end());
	else
		return iterator_cast<typename geometry_traits<TElem>::iterator>(
								section_container<TElem>().section_end(sInd));
}

//	const end
template <class TElem>
inline typename geometry_traits<TElem>::const_iterator
Selector::end() const
{
	const int sInd = get_section_index<TElem>();
	if(sInd < 0)
		return iterator_cast<typename geometry_traits<TElem>::const_iterator>(
									section_container<TElem>().end());
	else
		return iterator_cast<typename geometry_traits<TElem>::const_iterator>(
								section_container<TElem>().section_end(sInd));
}

template <class TElem>
TElem*
Selector::front()
{
	const int sInd = get_section_index<TElem>();
	return static_cast<TElem*>(section_container<TElem>().front(sInd));
}

template <class TElem>
TElem*
Selector::back()
{
	const int sInd = get_section_index<TElem>();
	return static_cast<TElem*>(section_container<TElem>().back(sInd));
}

////////////////////////////////////////
//	for compatibility with MGSelector


inline size_t Selector::
num_levels() const
{
	return 1;
}

inline uint Selector::
num(size_t) const
{
	return num();
}

template <class TElem>
inline size_t Selector::
num(size_t) const
{
	return num<TElem>();
}

inline bool Selector::
empty(size_t) const
{
	return empty();
}

template <class TElem>
inline bool Selector::
empty(size_t) const
{
	return empty<TElem>();
}

template <class TElem>
inline typename geometry_traits<TElem>::iterator
Selector::begin(size_t)
{
	return begin<TElem>();
}

//	end
///	calls end<TElem>();
template <class TElem>
inline typename geometry_traits<TElem>::iterator
Selector::end(size_t)
{
	return end<TElem>();
}

template <class TElem>
typename Grid::traits<TElem>::SectionContainer&
Selector::
section_container()
{
	return SectionContainerSelector<typename geometry_traits<TElem>::geometric_base_object>::
			section_container(m_vertices, m_edges, m_faces, m_volumes);
}


template <class TElem>
const typename Grid::traits<TElem>::SectionContainer&
Selector::
section_container() const
{
	return SectionContainerSelector<typename geometry_traits<TElem>::geometric_base_object>::
			section_container(m_vertices, m_edges, m_faces, m_volumes);
}

}//	end of namespace

#endif
