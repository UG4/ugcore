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
inline Selector::SectionContainer&
Selector::get_section_container()
{
	const int baseObjID = geometry_traits<TElem>::BASE_OBJECT_TYPE_ID;

	assert((baseObjID >= 0) && (baseObjID < NUM_GEOMETRIC_BASE_OBJECTS) &&
			"no section-container associated with TElem.");
	return m_elements[baseObjID];
}

template <class TElem>
inline int
Selector::get_section_index()
{
	return geometry_traits<TElem>::SHARED_PIPE_SECTION;
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
			get_section_container<TElem>().clear();
		else
			get_section_container<TElem>().clear_section(sInd);
	}
}

template <class TElem>
inline uint 
Selector::num()
{
	const int sInd = get_section_index<TElem>();
	if(sInd < 0)
		return get_section_container<TElem>().num_elements();
	else
		return get_section_container<TElem>().num_elements(sInd);
}

inline uint 
Selector::num()
{
	return num<VertexBase>() + num<EdgeBase>() + num<Face>() + num<Volume>();
}

//	empty
inline bool 
Selector::empty()
{
	return num() == 0;
}

template <class TElem>
inline bool 
Selector::empty()
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
								get_section_container<TElem>().begin());
	else
		return iterator_cast<typename geometry_traits<TElem>::iterator>(
					get_section_container<TElem>().section_begin(sInd));
}

//	end
template <class TElem>
inline typename geometry_traits<TElem>::iterator
Selector::end()
{
	const int sInd = get_section_index<TElem>();
	if(sInd < 0)
		return iterator_cast<typename geometry_traits<TElem>::iterator>(
									get_section_container<TElem>().end());
	else
		return iterator_cast<typename geometry_traits<TElem>::iterator>(
								get_section_container<TElem>().section_end(sInd));
}


}//	end of namespace

#endif
