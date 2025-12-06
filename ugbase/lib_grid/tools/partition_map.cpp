/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
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

size_t PartitionMap::num_target_procs() const {return m_targetProcs.size();}

int PartitionMap::get_target_proc(size_t index) const {
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

int PartitionMap::find_target_proc(int procRank) const {
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
