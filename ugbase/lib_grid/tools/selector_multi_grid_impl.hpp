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
inline MGSelector::SectionContainer&
MGSelector::get_section_container(int level)
{
	assert(level >= 0 && "bad level index.");

	const int baseObjID = geometry_traits<TElem>::BASE_OBJECT_TYPE_ID;

	assert((baseObjID >= 0) && (baseObjID < NUM_GEOMETRIC_BASE_OBJECTS) &&
			"no section-container associated with TElem.");

	level_required(level);
	return m_levels[level]->m_elements[baseObjID];
}

template <class TElem>
inline int
MGSelector::get_section_index()
{
	return geometry_traits<TElem>::SHARED_PIPE_SECTION;
}
inline void
MGSelector::level_required(int level)
{
//	create new SectionContainers and push them to the list,
//	until there are enough of them.
	while((int)m_levels.size() <= level)
		m_levels.push_back(new Level);
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
			get_section_container<TElem>(level).clear();
		else
			get_section_container<TElem>(level).clear_section(sInd);
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
inline uint 
MGSelector::num(int level)
{
	const int sInd = get_section_index<TElem>();
	if(sInd < 0)
		return get_section_container<TElem>(level).num_elements();
	else
		return get_section_container<TElem>(level).num_elements(sInd);
}

inline uint 
MGSelector::num(int level)
{
	return num<VertexBase>(level) + num<EdgeBase>(level)
			+ num<Face>(level) + num<Volume>(level);
}

template <class TElem>
inline uint 
MGSelector::num()
{
	uint n = 0;
	for(uint i = 0; i < num_levels(); ++i)
		n += num<TElem>(i);
	return n;
}

inline uint 
MGSelector::num()
{
	return num<VertexBase>() + num<EdgeBase>()
			+ num<Face>() + num<Volume>();
}

//	empty
inline bool 
MGSelector::empty(int level)
{
	return num(level) == 0;
}

template <class TElem>
inline bool 
MGSelector::empty(int level)
{
	return num<TElem>(level) == 0;
}

inline bool 
MGSelector::empty()
{
	return num() == 0;
}

template <class TElem>
inline bool 
MGSelector::empty()
{
	return num<TElem>() == 0;
}

//	begin
template <class TElem>
inline typename geometry_traits<TElem>::iterator
MGSelector::begin(int level)
{
	const int baseObjID = geometry_traits<TElem>::BASE_OBJECT_TYPE_ID;
	assert((baseObjID >= 0) && (baseObjID < NUM_GEOMETRIC_BASE_OBJECTS));
	
	const int sInd = get_section_index<TElem>();

	level_required(level);
	if(sInd < 0)
		return iterator_cast<typename geometry_traits<TElem>::iterator>(
			m_levels[level]->m_elements[baseObjID].begin());
	else
		return iterator_cast<typename geometry_traits<TElem>::iterator>(
			m_levels[level]->m_elements[baseObjID].section_begin(sInd));
}

//	end
template <class TElem>
inline typename geometry_traits<TElem>::iterator
MGSelector::end(int level)
{
	const int baseObjID = geometry_traits<TElem>::BASE_OBJECT_TYPE_ID;
	assert((baseObjID >= 0) && (baseObjID < NUM_GEOMETRIC_BASE_OBJECTS));
	
	const int sInd = get_section_index<TElem>();

	level_required(level);
	if(sInd < 0)
		return iterator_cast<typename geometry_traits<TElem>::iterator>(
			m_levels[level]->m_elements[baseObjID].end());
	else
		return iterator_cast<typename geometry_traits<TElem>::iterator>(
			m_levels[level]->m_elements[baseObjID].section_end(sInd));
}


}//	end of namespace

#endif
