// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y10 m02 d15

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	...
////////////////////////////////////////////////////////////////////////

#ifndef __H__LIBGRID__SELECTOR_MULTI_GRID_IMPL__
#define __H__LIBGRID__SELECTOR_MULTI_GRID_IMPL__

#include <cassert>

namespace ug
{

template <class TElem>
inline int
MGSelector::get_section_index() const
{
	return geometry_traits<TElem>::CONTAINER_SECTION;
}

inline void
MGSelector::level_required(int level)
{
//	create new SectionContainers and push them to the list,
//	until there are enough of them.
	while((int)m_levels.size() <= level){
		add_level();
	}
}

template <class TElem>
inline void
MGSelector::clear(int level)
{
	if(m_pGrid){
	//	mark all elements as deselected
		typename geometry_traits<TElem>::iterator iter;
		for(iter = begin<TElem>(level); iter != end<TElem>(level); ++iter)
			mark_deselected(*iter);

	//	clear the section
		const int sInd = get_section_index<TElem>();
		if(sInd < 0)
			section_container<TElem>(level).clear();
		else
			section_container<TElem>(level).clear_section(sInd);
	}
}

template <class TElem>
inline void
MGSelector::clear()
{
	if(m_pGrid){
		for(size_t i = 0; i < num_levels(); ++i)
			clear<TElem>(i);
	}
}

template <class TElem>
inline size_t
MGSelector::num(int level) const
{
	const int sInd = get_section_index<TElem>();
	if(level < (int)num_levels()){
		if(sInd < 0)
			return section_container<TElem>(level).num_elements();
		else
			return section_container<TElem>(level).num_elements(sInd);
	}
	return 0;
}

inline size_t
MGSelector::num(int level) const
{
	return num<Vertex>(level) + num<EdgeBase>(level)
			+ num<Face>(level) + num<Volume>(level);
}

template <class TElem>
inline size_t
MGSelector::num() const
{
	size_t n = 0;
	for(size_t i = 0; i < num_levels(); ++i)
		n += num<TElem>(i);
	return n;
}

inline size_t
MGSelector::num() const
{
	return num<Vertex>() + num<EdgeBase>()
			+ num<Face>() + num<Volume>();
}

//	empty
inline bool 
MGSelector::empty(int level) const
{
	return num(level) == 0;
}

template <class TElem>
inline bool 
MGSelector::empty(int level) const
{
	return num<TElem>(level) == 0;
}

inline bool 
MGSelector::empty() const
{
	return num() == 0;
}

template <class TElem>
inline bool 
MGSelector::empty() const
{
	return num<TElem>() == 0;
}

//	begin
template <class TElem>
inline typename MGSelector::traits<TElem>::iterator
MGSelector::begin()
{
//	get the begin iterator of the first level which is not empty.
//	if all levels are empty the returned iterator is equal to the end iterator
//	of the top-level.
	size_t lvl = 0;
	while(empty<TElem>(lvl) && (lvl < top_level()))
		++lvl;

	return typename MGSelector::traits<TElem>::iterator(this, lvl, begin<TElem>(lvl));
}

template <class TElem>
inline typename MGSelector::traits<TElem>::const_iterator
MGSelector::begin() const
{
	//	get the begin iterator of the first level which is not empty.
//	if all levels are empty the returned iterator is equal to the end iterator
//	of the top-level.
	size_t lvl = 0;
	while(empty<TElem>(lvl) && (lvl < top_level()))
		++lvl;

	return typename MGSelector::traits<TElem>::const_iterator(this, lvl, begin<TElem>(lvl));
}

template <class TElem>
inline typename MGSelector::traits<TElem>::level_iterator
MGSelector::begin(int level)
{
	const int sInd = get_section_index<TElem>();

	if(sInd < 0)
		return iterator_cast<typename geometry_traits<TElem>::iterator>(
								section_container<TElem>(level).begin());
	else
		return iterator_cast<typename geometry_traits<TElem>::iterator>(
					section_container<TElem>(level).section_begin(sInd));
}

template <class TElem>
inline typename MGSelector::traits<TElem>::const_level_iterator
MGSelector::begin(int level) const
{
	const int sInd = get_section_index<TElem>();

	if(sInd < 0)
		return iterator_cast<typename geometry_traits<TElem>::const_iterator>(
								section_container<TElem>(level).begin());
	else
		return iterator_cast<typename geometry_traits<TElem>::const_iterator>(
					section_container<TElem>(level).section_begin(sInd));
}


template <class TElem>
inline typename MGSelector::traits<TElem>::iterator
MGSelector::end()
{
	size_t l = top_level();
	return typename MGSelector::traits<TElem>::iterator(this, l, end<TElem>(l));
}

template <class TElem>
inline typename MGSelector::traits<TElem>::const_iterator
MGSelector::end() const
{
	size_t l = top_level();
	return typename MGSelector::traits<TElem>::const_iterator(this, l, end<TElem>(l));
}

template <class TElem>
inline typename MGSelector::traits<TElem>::level_iterator
MGSelector::end(int level)
{
	const int sInd = get_section_index<TElem>();

	if(sInd < 0)
		return iterator_cast<typename geometry_traits<TElem>::iterator>(
								section_container<TElem>(level).end());
	else
		return iterator_cast<typename geometry_traits<TElem>::iterator>(
					section_container<TElem>(level).section_end(sInd));
}

template <class TElem>
inline typename MGSelector::traits<TElem>::const_level_iterator
MGSelector::end(int level) const
{
	const int sInd = get_section_index<TElem>();

	if(sInd < 0)
		return iterator_cast<typename geometry_traits<TElem>::const_iterator>(
								section_container<TElem>(level).end());
	else
		return iterator_cast<typename geometry_traits<TElem>::const_iterator>(
					section_container<TElem>(level).section_end(sInd));
}

template <class TElem>
TElem*
MGSelector::front(int level)
{
	const int sInd = get_section_index<TElem>();
	return static_cast<TElem*>(section_container<TElem>(level).front(sInd));
}

template <class TElem>
TElem*
MGSelector::back(int level)
{
	const int sInd = get_section_index<TElem>();
	return static_cast<TElem*>(section_container<TElem>(level).back(sInd));
}

template <class TElem>
typename Grid::traits<TElem>::SectionContainer&
MGSelector::
section_container(int level)
{
	assert(level >= 0 && "bad level index.");
	level_required(level);
	Level* lev = m_levels[level];
	return SectionContainerSelector<typename geometry_traits<TElem>::grid_base_object>::
			section_container(lev->m_vertices, lev->m_edges, lev->m_faces, lev->m_volumes);
}


template <class TElem>
const typename Grid::traits<TElem>::SectionContainer&
MGSelector::
section_container(int level) const
{
	assert((level >= 0) && (level < (int)m_levels.size()) && "bad level index.");
	const Level* lev = m_levels[level];
	return SectionContainerSelector<typename geometry_traits<TElem>::grid_base_object>::
			section_container(lev->m_vertices, lev->m_edges, lev->m_faces, lev->m_volumes);
}

}//	end of namespace

#endif
