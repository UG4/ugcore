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

#ifndef __H__LIB_GRID__DISTRIBUTED_GRID_IMPL__
#define __H__LIB_GRID__DISTRIBUTED_GRID_IMPL__

#include <vector>
#include <utility>

namespace ug
{

template <class TElem>
bool DistributedGridManager::
is_interface_element(TElem* elem)
{
	return elem_info(elem).get_status() & ES_IN_INTERFACE;
}

template<class TElem>
inline bool DistributedGridManager::
is_in_horizontal_interface(TElem* elem) const
{
	byte status = get_status(elem);
	return 	(status & (ES_H_MASTER | ES_H_SLAVE)) != 0;
}

template<class TElem>
inline bool DistributedGridManager::
is_in_vertical_interface(TElem* elem) const
{
	byte status = get_status(elem);
	return 	(status & (ES_V_MASTER | ES_V_SLAVE)) != 0;
}

template<class TElem>
inline bool DistributedGridManager::
is_ghost(TElem* elem) const
{
	byte status = get_status(elem);
	return 	(status & (ES_V_MASTER | ES_H_MASTER | ES_H_SLAVE))
			== ES_V_MASTER;

	//would require update_ghost_states
	//return contains_status(elem, ES_GHOST);
}

template <class TElem>
void DistributedGridManager::
collect_interface_entries(
				std::vector<std::pair<int, size_t> >& vEntriesOut,
				TElem* elem, byte statusType, bool clearContainer)
{
//TODO: make sure that the localIDs match the position at which
//		an element is stored in the interface
	typedef ElementInfo<TElem> ElemInfo;
	ElemInfo& info = elem_info(elem);

	if(clearContainer)
		vEntriesOut.clear();

	for(typename ElemInfo::EntryIterator iter = info.entries_begin();
		iter != info.entries_end(); ++iter)
	{
		if((info.get_interface_type(iter) & statusType) == statusType){
			vEntriesOut.push_back(std::make_pair(info.get_target_proc(iter),
											info.get_local_id(iter)));
		}
	}
}

}// end of namespace

#endif
