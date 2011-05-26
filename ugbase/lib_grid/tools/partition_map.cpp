// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 25.05.2011 (m,d,y)
 
#include "partition_map.h"

using namespace std;

namespace ug{

////////////////////////////////////////////////////////////////////////////////
//	IMPLEMENTATION OF PartitionMap
void PartitionMap::clear()
{
	m_targetProcs.clear();
	m_shPartitions.clear();
}

void PartitionMap::assign_grid(Grid& grid)
{
	if(&grid != m_shPartitions.get_assigned_grid())
		m_shPartitions.assign_grid(grid);
}

SubsetHandler& PartitionMap::get_partition_handler()
{return m_shPartitions;}

void PartitionMap::add_target_proc(int targetProcRank)
{m_targetProcs.push_back(targetProcRank);}

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

///	returns the index at which the given process lies. -1 if it doesn't exist.
int PartitionMap::find_target_proc(int procRank)
{
	for(size_t i = 0; i < m_targetProcs.size(); ++i){
		if(m_targetProcs[i] == procRank)
			return i;
	}
	return -1;
}

}// end of namespace
