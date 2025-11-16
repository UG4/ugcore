/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

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
inline void ISelector::select(TElem* elem, byte_t status){
	if(status != 0){
		if(!is_selected(elem)){
			add_to_list(elem);
		}
		mark_selected(elem, status);
	}
	else
		deselect(elem);
}

inline void ISelector::select(GridObject* elem, byte_t status){
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
inline void ISelector::select(TIterator iterBegin, TIterator iterEnd, byte_t status)
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


byte_t ISelector::get_selection_status(GridObject* elem) const{
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
