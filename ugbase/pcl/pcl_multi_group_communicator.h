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

#ifndef __H__pcl_multi_group_communicator
#define __H__pcl_multi_group_communicator

#include <vector>

#include "pcl_process_communicator.h"


namespace pcl {

///	communicator for simultaneous data exchange between many small groups
/**	The MultiGroupCommunicator allows for a process to be member in various small groups
 * and to perform collective operations on those groups.
 *
 * \warning The communicator is currently in development and not ready to use!
 *			Segfaults may arise during or after application. Debugging required!
 */
class MultiGroupCommunicator {
public:
	MultiGroupCommunicator (ProcessCommunicator com = ProcessCommunicator());
	MultiGroupCommunicator (const std::vector<bool>& participates,
	                        ProcessCommunicator com = ProcessCommunicator());

	void reinit (const std::vector<bool>& participates);

	bool participates (size_t igrp) const;
	size_t num_groups () const;
	
	size_t num_memberships () const;
	size_t membership_group_index (size_t imem) const;

///	performs an allreduce between all groups
/** sendBuf and recvBuf have to be of length
 * 'countPerGroup' * num_memberships()
 */
	template<typename T>
	void allreduce (const T *sendBuf, T *recvBuf,
				    size_t countPerGroup, ReduceOperation op) const;

private:
	ProcessCommunicator	m_com;
	std::vector<bool> m_participates; ///< size: #groups
	std::vector<size_t>	m_memberships; ///< size: #memberships. Holds indices to groups in which the process participates
	std::vector<int> m_groupMembers; ///< size: m_groupOffsets.back(). Consecutively holds proc-indices of each group in which the process participates.
	std::vector<size_t>	m_groupOffsets;	///< size: #memberships+1. Offset of each group in m_groupMembers. The last entry always holds m_groupMembers.size().
	std::map<int, ug::BinaryBuffer>	m_binBufs; ///< used for sending/receiving data from/to processes

};

}//	end of namespace


////////////////////////////////
//	include implementation
#include "pcl_multi_group_communicator_impl.hpp"


#endif