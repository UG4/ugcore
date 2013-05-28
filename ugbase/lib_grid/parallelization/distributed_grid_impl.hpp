// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y10 m03 d31

#ifndef __H__LIB_GRID__DISTRIBUTED_GRID_IMPL__
#define __H__LIB_GRID__DISTRIBUTED_GRID_IMPL__

#include <vector>

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
			vEntriesOut.push_back(make_pair(info.get_target_proc(iter),
											info.get_local_id(iter)));
		}
	}
}

}// end of namespace

#endif
