// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y10 m02 d15

#ifndef __H__LIBGRID__SELECTOR_INTERFACE_IMPL__
#define __H__LIBGRID__SELECTOR_INTERFACE_IMPL__

namespace ug
{
inline bool
ISelector::elements_are_supported(uint shElements) const
{
	return (m_supportedElements & shElements) == shElements;
}

template <class TElem>
inline void ISelector::select(TElem* elem){
	if(!is_selected(elem))
		mark_selected(elem, add_to_list(elem));
}

inline void ISelector::select(GeometricObject* elem){
	int elemID = elem->base_object_type_id();
	switch(elemID){
		case VERTEX:
			select(static_cast<VertexBase*>(elem));
			break;
		case EDGE:
			select(static_cast<EdgeBase*>(elem));
			break;
		case FACE:
			select(static_cast<Face*>(elem));
			break;
		case VOLUME:
			select(static_cast<Volume*>(elem));
			break;
		default:
			LOG("  ERROR: Bad Element Type in ISelector::select. Aborting.\n");
			assert(0);
			break;
	}
}

template <class TIterator>
inline void ISelector::select(TIterator iterBegin, TIterator iterEnd)
{
	while(iterBegin != iterEnd){
		select(*iterBegin);
		iterBegin++;
	}
}


template <class TElem>
inline void ISelector::deselect(TElem* elem){
	if(is_selected(elem)){
		erase_from_list(elem);
		mark_deselected(elem);
	}
}

inline void ISelector::deselect(GeometricObject* elem){
	int elemID = elem->base_object_type_id();
	switch(elemID){
		case VERTEX:
			deselect(static_cast<VertexBase*>(elem));
			break;
		case EDGE:
			deselect(static_cast<EdgeBase*>(elem));
			break;
		case FACE:
			deselect(static_cast<Face*>(elem));
			break;
		case VOLUME:
			deselect(static_cast<Volume*>(elem));
			break;
	}
}

template <class TIterator>
inline void ISelector::deselect(TIterator iterBegin, TIterator iterEnd)
{
	while(iterBegin != iterEnd){
		typename TIterator::value_type& val = *iterBegin;
		++iterBegin;
		deselect(val);
	}
}


bool ISelector::is_selected(GeometricObject* elem) const{
	int elemID = elem->base_object_type_id();
	switch(elemID){
		case VERTEX:
			return is_selected(static_cast<VertexBase*>(elem));
		case EDGE:
			return is_selected(static_cast<EdgeBase*>(elem));
		case FACE:
			return is_selected(static_cast<Face*>(elem));
		case VOLUME:
			return is_selected(static_cast<Volume*>(elem));
	}
	return false;
}

}//	end of namespace

#endif
