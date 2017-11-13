/*
 * Copyright (c) 2017:  G-CSC, Goethe University Frankfurt
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

#include "common/log.h"
#include "pcl_multi_group_communicator.h"

using namespace std;

namespace pcl {

MultiGroupCommunicator::
MultiGroupCommunicator (ProcessCommunicator com) :
	m_com (com)
{

}

MultiGroupCommunicator::
MultiGroupCommunicator (const std::vector<bool>& participates,
                        ProcessCommunicator com) :
	m_com (com)
{
	reinit (participates);
}

void MultiGroupCommunicator::
reinit (const std::vector<bool>& participates)
{
	m_participates = participates;

//	build group list
	vector<int>	memIndFromGrpInd (participates.size(), -1);

	m_memberships.clear();
	for(size_t i = 0; i < participates.size(); ++i){
		if(participates[i]){
			memIndFromGrpInd [i] = (int) m_memberships.size();
			m_memberships.push_back (i);
		}
	}


//	distribute membership infos between all processes
	vector<int> memListSizes;
	vector<int> memListOffsets;
	vector<size_t> memLists;
	m_com.allgatherv (memLists, m_memberships, &memListSizes, &memListOffsets);


//	collect sizes of all groups in which this process participates
	vector<size_t>	grpSizes (m_memberships.size(), 0);

	for(size_t iproc = 0; iproc < m_com.size(); ++iproc){
		const int numMems = memListSizes [iproc];
		const int o = memListOffsets [iproc];
		for(int i = 0; i < numMems; ++i){
			const size_t igrp = memLists [o + i];
			if (participates [igrp]){
				++ grpSizes[memIndFromGrpInd[igrp]];
			}
		}
	}


//	compute local offsets to from sizes
	m_groupOffsets.resize (m_memberships.size() + 1);
	m_groupOffsets[0] = 0;

	for(size_t imem = 0; imem < m_memberships.size(); ++imem){
		m_groupOffsets [imem + 1] = m_groupOffsets [imem] + grpSizes [imem];
	}


//	fill m_groupMembers by collecting all processes which participate in a group
	m_groupMembers.resize (m_groupOffsets.back());
	vector<size_t> grpMemFillCount (m_memberships.size(), 0);

	for(size_t iproc = 0; iproc < m_com.size(); ++iproc){
		const size_t s = memListSizes [iproc];
		const size_t o = memListOffsets [iproc];

		for(size_t i = 0; i < s; ++i){
			const size_t igrp = memLists [o + i];
			if(participates [igrp]) {
				const int imem = memIndFromGrpInd [igrp];
				m_groupMembers [m_groupOffsets [imem] + grpMemFillCount [imem]] = iproc;
				++ grpMemFillCount [imem];
			}
		}
	}
}

bool MultiGroupCommunicator::
participates (size_t igrp) const
{
	return m_participates [igrp];
}


size_t MultiGroupCommunicator::
num_groups () const
{
	return m_participates.size();
}


size_t MultiGroupCommunicator::
num_memberships () const
{
	return m_memberships.size();
}

size_t MultiGroupCommunicator::
membership_group_index (size_t imem) const
{
	return m_memberships [imem];
}


}//	end of namespace
