// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 25.05.2011 (m,d,y)
 
#include "partition_map.h"

using namespace std;

namespace ug{

////////////////////////////////////////////////////////////////////////////////
//	IMPLEMENTATION OF PartitionMap
PartitionMap::PartitionMap()
{
	m_shPartitions = make_sp(new SubsetHandler());
}

void PartitionMap::clear()
{
	m_targetProcs.clear();
	m_shPartitions->clear();
}

void PartitionMap::assign_grid(Grid& grid)
{
	if(&grid != m_shPartitions->grid())
		m_shPartitions->assign_grid(grid);
}

SmartPtr<SubsetHandler> PartitionMap::get_partition_handler()
{return m_shPartitions;}

void PartitionMap::add_target_proc(int tarProcRank)
{m_targetProcs.push_back(tarProcRank);}

void PartitionMap::add_target_procs(int first, int num)
{
	for(int i = 0; i < num; ++i)
		add_target_proc(first + i);
}

size_t PartitionMap::num_target_procs()
{return m_targetProcs.size();}

int PartitionMap::get_target_proc(size_t index)
{
	if(index < m_targetProcs.size())
		return m_targetProcs[index];
	UG_LOG("BAD INDEX in PartitionMap::get_target_proc: " << index);
	if(num_target_procs() > 0){
		UG_LOG("    Max valid index: " << num_target_procs() - 1 << endl);
	}
	else{
		UG_LOG("    No target processes available.\n");
	}
	return -1;
}

int* PartitionMap::get_target_procs()
{return &m_targetProcs.front();}

std::vector<int>& PartitionMap::get_target_proc_vec()
{return m_targetProcs;}

bool PartitionMap::change_target_proc(size_t index, int newRank)
{
//	make sure that the given index is valid
	if(index >= num_target_procs()){
		UG_LOG("WARNING in PartitionMap::change_target_proc: Bad index given.\n");
		return false;
	}

	m_targetProcs[index] = newRank;
	return true;
}

int PartitionMap::find_target_proc(int procRank)
{
	for(size_t i = 0; i < m_targetProcs.size(); ++i){
		if(m_targetProcs[i] == procRank)
			return i;
	}
	return -1;
}

void PartitionMap::shift_target_procs(int offset)
{
	for(size_t i = 0; i < m_targetProcs.size(); ++i){
		m_targetProcs[i] += offset;
	}
}

}// end of namespace
