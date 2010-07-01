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
is_ghost(TElem* elem)
{
	byte status = get_status(elem);
	if((status & ES_VERTICAL_MASTER))
		UG_LOG("vm-");
	if((status & ES_MASTER))
		UG_LOG("hm-");
	if((status & ES_SLAVE))
		UG_LOG("hs-");
	return 	(status & (ES_VERTICAL_MASTER | ES_MASTER | ES_SLAVE))
			== ES_VERTICAL_MASTER;
}

template <class TElem>
void DistributedGridManager::
collect_interface_entries(
				std::vector<std::pair<int, size_t> >& vEntriesOut,
				TElem* elem)
{
//TODO: make sure that the localIDs match the position at which
//		an element is stored in the interface
	typedef ElementInfo<TElem> ElemInfo;
	ElemInfo& info = elem_info(elem);

	vEntriesOut.clear();
	
	for(typename ElemInfo::EntryIterator iter = info.entries_begin();
		iter != info.entries_end(); ++iter)
	{
		vEntriesOut.push_back(make_pair(info.get_target_proc(iter),
										info.get_local_id(iter)));
	}
}

}// end of namespace

#endif
