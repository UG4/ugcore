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
inline void ISelector::select(TElem* elem, byte status){
	if(status != 0){
		if(!is_selected(elem)){
			add_to_list(elem);
		}
		mark_selected(elem, status);
	}
	else
		deselect(elem);
}

inline void ISelector::select(GridObject* elem, byte status){
	int elemID = elem->base_object_id();
	switch(elemID){
		case VERTEX:
			select(static_cast<Vertex*>(elem), status);
			break;
		case EDGE:
			select(static_cast<Edge*>(elem), status);
			break;
		case FACE:
			select(static_cast<Face*>(elem), status);
			break;
		case VOLUME:
			select(static_cast<Volume*>(elem), status);
			break;
		default:
			LOG("  ERROR: Bad Element Type in ISelector::select. Aborting.\n");
			assert(0);
			break;
	}
}

template <class TIterator>
inline void ISelector::select(TIterator iterBegin, TIterator iterEnd, byte status)
{
	while(iterBegin != iterEnd){
		select(*iterBegin, status);
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

inline void ISelector::deselect(GridObject* elem){
	int elemID = elem->base_object_id();
	switch(elemID){
		case VERTEX:
			deselect(static_cast<Vertex*>(elem));
			break;
		case EDGE:
			deselect(static_cast<Edge*>(elem));
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
		typename TIterator::value_type val = *iterBegin;
		++iterBegin;
		deselect(val);
	}
}


byte ISelector::get_selection_status(GridObject* elem) const{
	int elemID = elem->base_object_id();
	switch(elemID){
		case VERTEX:
			return get_selection_status(static_cast<Vertex*>(elem));
		case EDGE:
			return get_selection_status(static_cast<Edge*>(elem));
		case FACE:
			return get_selection_status(static_cast<Face*>(elem));
		case VOLUME:
			return get_selection_status(static_cast<Volume*>(elem));
	}
	return 0;
}

}//	end of namespace

#endif
