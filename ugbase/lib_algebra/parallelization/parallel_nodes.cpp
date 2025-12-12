/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
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

#include "parallel_nodes.h"

using namespace std;

namespace ug {

ParallelNodes::ParallelNodes(ConstSmartPtr<AlgebraLayouts> layout, size_t s)
{
	m_layout = layout;
	GenerateGlobalAlgebraIDs(layout->comm(), m_localToGlobal, s, layout->master(), layout->slave());

	for(size_t i=0; i<m_localToGlobal.size(); ++i)
	{
		if(m_localToGlobal[i].first != pcl::ProcRank())
			m_globalToLocal.insert(pair (m_localToGlobal[i], i));
	}

	m_OLtype.resize(s);
	for(auto iter = layout->master().begin(); iter != layout->master().end(); ++iter)
	{
		const IndexLayout::Interface &interface = layout->master().interface(iter);
		for(auto iter = interface.begin(); iter != interface.end(); ++iter)
		{
			size_t i = interface.get_element(iter);
			m_OLtype[i].set_master();
		}
	}
	for(auto iter = layout->slave().begin(); iter != layout->slave().end(); ++iter)
	{
		const IndexLayout::Interface &interface = layout->slave().interface(iter);
		for(auto iter = interface.begin(); iter != interface.end(); ++iter)
		{
			size_t i = interface.get_element(iter);
			UG_ASSERT(m_OLtype[i].is_master() == false, i);
			m_OLtype[i].set_slave();
		}
	}

	AddLayout(totalMasterLayout, layout->master());
	AddLayout(totalSlaveLayout, layout->slave());

	for(auto it = layout->master().begin(); it != layout->master().end(); ++it)
		masterPIDs.insert(layout->master().proc_id(it));
	for(auto it = layout->slave().begin(); it != layout->slave().end(); ++it)
		slavePIDs.insert(layout->slave().proc_id(it));
	create_mark_map(layout->master());

	m_originalSize = s;
}

size_t ParallelNodes::get_local_index_if_available(const AlgebraID &globalIndex, bool &bHasIndex) const
{
	// UG_DLOG(LIB_ALG_MATRIX, 4, "get_index_if_available for global index " << globalIndex << ": ");
	if(globalIndex.first == pcl::ProcRank())
	{
		// UG_DLOG(LIB_ALG_MATRIX, 4, "is on this processor\n");
		bHasIndex = true;
		return globalIndex.second;
	}
	auto it = m_globalToLocal.find(globalIndex);
	if(it == m_globalToLocal.end())
	{
		// UG_DLOG(LIB_ALG_MATRIX, 4, "index not found\n");
		bHasIndex = false;
		return (size_t) -1;
	}
	else
	{
		// UG_DLOG(LIB_ALG_MATRIX, 4, "index is " << it->second << "\n");
		bHasIndex = true;
		return it->second;
	}
}

size_t ParallelNodes::get_local_index_or_create_new(const AlgebraID &globalIndex, int distanceToMasterOrInner, bool &bCreated)
{
	// UG_DLOG(LIB_ALG_MATRIX, 4, "get_index_or_create_new for global index " << globalIndex.first << " | " << globalIndex.second << ": ");
	if(globalIndex.first == pcl::ProcRank())
	{
		// UG_DLOG(LIB_ALG_MATRIX, 4, "is on this processor\n");
		bCreated = false;
		return globalIndex.second;
	}
	else
	{
		pair<iterator, bool> ret = m_globalToLocal.insert(pair (globalIndex, m_localToGlobal.size()));
		if(ret.second)
		{
			UG_DLOG(LIB_ALG_MATRIX, 4, "created new index " << m_localToGlobal.size() << " (global Index = " << globalIndex << "\n");
			m_localToGlobal.push_back(globalIndex);
			m_OLtype.emplace_back(distanceToMasterOrInner);
			bCreated = true;
		}
		else
		{
			// UG_DLOG(LIB_ALG_MATRIX, 4, "already has index " << ret.first->second << "\n");
			bCreated = false;
		}
		return ret.first->second;
	}
}

size_t ParallelNodes::get_local_index_or_create_new(const AlgebraID &globalIndex, int distanceToMasterOrInner)
{
	bool b;
	return get_local_index_or_create_new(globalIndex, distanceToMasterOrInner, b);
}

size_t ParallelNodes::global_to_local(const AlgebraID &globalIndex) const
{
	bool hasIndex;
	size_t index = get_local_index_if_available(globalIndex, hasIndex);
	UG_ASSERT(hasIndex, "global id " << globalIndex << " has no associated local id");
	return index;
}


void ParallelNodes::insert_into_interface_sorted(vector<size_t> &v, IndexLayout::Interface &interface)
{
	sort_by_global_id(v);
	for(size_t i=0; i<v.size(); i++)
		interface.push_back(v[i]);
}

void ParallelNodes::insert_into_layout_sorted(map<int, set<size_t> > &m, IndexLayout &layout)
{
	for(auto it = m.begin(); it != m.end(); ++it)
	{
		int pid = it->first;
		vector<size_t> v;
		v.reserve(it->second.size());
		v.insert(v.end(), it->second.begin(), it->second.end());
		if(v.size() != 0)
			insert_into_interface_sorted(v, layout.interface(pid));
	}
}

void ParallelNodes::insert_into_layout_sorted(map<int, vector<size_t> > &m, IndexLayout &layout)
{
	for(auto it = m.begin(); it != m.end(); ++it)
	{
		int pid = it->first;
		vector<size_t> &v = it->second;
		if(v.size() != 0)
			insert_into_interface_sorted(v, layout.interface(pid));
	}
}

void ParallelNodes::sort_interface(IndexLayout::Interface &interface)
{
	vector<size_t> v;
	v.reserve(interface.size());
	for(auto iter = interface.begin(); iter != interface.end(); ++iter)
		v.push_back(interface.get_element(iter));

	IndexLayout::Interface interface2;
	insert_into_interface_sorted(v, interface2);
	interface2.swap(interface);
}

void ParallelNodes::sort_layout(IndexLayout &layout)
{
	for(auto iter = layout.begin(); iter != layout.end(); ++iter)
		sort_interface(layout.interface(iter));
}


void ParallelNodes::create_mark_map(const IndexLayout &masterLayout)
{
	PROFILE_FUNC();
	UG_DLOG(LIB_ALG_MATRIX, 4, "\n\nGenerateOverlapClass::CreateMarks\n");
	notified.clear();
	for(auto iter = masterLayout.begin(); iter != masterLayout.end(); ++iter)
	{
		const IndexLayout::Interface &interface = masterLayout.interface(iter);
		int pid = masterLayout.proc_id(iter);

		set<size_t> &mark = notified[pid];
		for(auto iter2 = interface.begin(); iter2 != interface.end(); ++iter2)
		{
			size_t localIndex = interface.get_element(iter2);
			mark.insert(localIndex);
		}
	}
}

size_t ParallelNodes::create_slave_node(const AlgebraID &globalID, int distanceToMasterOrInner)
{
	bool bCreated;
	size_t localIndex = get_local_index_or_create_new(globalID, distanceToMasterOrInner, bCreated);

	if(bCreated)
	{
		newSlaves[globalID.master_proc()].insert(localIndex);
		UG_DLOG(LIB_ALG_MATRIX, 4, "created new slave with GID " << globalID << " with local index " << localIndex << " and distanceToMasterOrInner=" << distanceToMasterOrInner << "\n");
	}
	return localIndex;
}

void ParallelNodes::create_node(const AlgebraID &globalID, size_t localIndex, int pid)
{
	if(globalID.master_proc() == pid) // obviously the owner knows that his node is on his processor
		return;

	if(notified[pid].find(localIndex) != notified[pid].end())
	{
		UG_DLOG(LIB_ALG_MATRIX, 4, " not sending notification because already send one\n");
		return;
	}

	/*IF_DEBUG(LIB_ALG_MATRIX, 4)
	{
		newNodesOn[pid].insert(globalID);
	}*/

	notified[pid].insert(localIndex);

	if(globalID.master_proc() == pcl::ProcRank())
	{
		// add i to mNewLayout[pid] if not already in some interface to pid.

		// only add nodes once to interface
		UG_DLOG(LIB_ALG_MATRIX, 4, " adding node " << localIndex << " to master interface to processor " << pid << "\n");
		newMasters[pid].insert(localIndex);
	}
	else
	{
		//	prepare to send a notification to owner of (globalColIndex.master_proc()) that globalColIndex has a copy on pid if not already done so

		NewSlaveNotification notification(globalID, pid);
		UG_DLOG(LIB_ALG_MATRIX, 4, "sending " << globalID.master_proc() << " " << notification << "\n");
		newSlaveNotifications[globalID.master_proc()].push_back(notification);
	}
}

	/*void check_new_nodes(pcl::InterfaceCommunicator<IndexLayout> &communicator,
				IndexLayout &sendingLayout, IndexLayout &receivingLayout)
	{
		set<int> sendingPIDs;
		for(IndexLayout::iterator it = sendingLayout.begin(); it != sendingLayout.end(); ++it)
			sendingPIDs.insert(sendingLayout.proc_id(it));

		set<int> receivingPIDs;
		for(IndexLayout::iterator it = receivingLayout.begin(); it != receivingLayout.end(); ++it)
			receivingPIDs.insert(receivingLayout.proc_id(it));
		check_new_nodes(commnunicator, sendingPIDs, receivingPIDs);
	}

	void check_new_nodes(vector<int> sendingPIDs, vector<int> receivingPIDs)
	{
		for(size_t i=0; i<sendingPIDs.size(); i++)
		{
			BinaryBuffer buf;
			int pid = sendingPIDs[i];
			Serialize(buf, newNodesOn[pid]);
			communicator.send_raw(pid, buf.buffer(), buf.write_pos(), false);
		}

		for(size_t i=0; i<receivingPIDs.size(); i++)
			communicator.receive_raw(pid, newNodesBufferMap[receivingPIDs[i]]);

		communicator.communicate();
	}*/

void ParallelNodes::issue(pcl::InterfaceCommunicator<IndexLayout> &communicator)
{
	UG_DLOG(LIB_ALG_MATRIX, 4, "NewLayoutCreator::issue\n");
	// notifications for new Master Nodes
	for(auto it = slavePIDs.begin(); it != slavePIDs.end(); ++it)
	{
		BinaryBuffer buf;
		int pid = *it;
		for(size_t i=0; i<newSlaveNotifications[pid].size(); i++)
			Serialize(buf, newSlaveNotifications[pid][i]);
		UG_DLOG(LIB_ALG_MATRIX, 4, " sending " << newSlaveNotifications[pid].size() << " notifications to processor " << pid << " (" << buf.write_pos() << " bytes)\n");
		communicator.send_raw(pid, buf.buffer(), buf.write_pos(), false);
	}
	newSlaveNotifications.clear();

	for(auto it = masterPIDs.begin(); it != masterPIDs.end(); ++it)
	{
		int pid = *it;
		UG_DLOG(LIB_ALG_MATRIX, 4, " issue receive from processor " << pid << "\n");
		communicator.receive_raw(pid, notificationBufferMap[pid]);
	}
}

void ParallelNodes::process()
{
	UG_DLOG(LIB_ALG_MATRIX, 4, "NewLayoutCreator::process\n");
	for(auto it = masterPIDs.begin(); it != masterPIDs.end(); ++it)
	{
		int pid = *it;
		BinaryBuffer &buf = notificationBufferMap[pid];
		UG_DLOG(LIB_ALG_MATRIX, 4, "received " << buf.write_pos() << " bytes of notification from pid " << pid << ":\n");
		while(!buf.eof())
		{
			NewSlaveNotification notification;
			Deserialize(buf, notification);
			UG_ASSERT(notification.id.master_proc() == pcl::ProcRank(), notification.id << ", pid = " << pcl::ProcRank());
			UG_DLOG(LIB_ALG_MATRIX, 4, notification << "\n");

			set<size_t> &mark = notified[notification.newSlaveOnPID];
			if(mark.find(notification.id.index_on_master()) == mark.end())
			{
				newMasters[notification.newSlaveOnPID].insert(notification.id.index_on_master());
				mark.insert(notification.id.index_on_master());
			}

		}
	}
	notificationBufferMap.clear();

	/*for(BufferMap::iterator it = newNodesBufferMap.begin(); it != newNodesBufferMap.end(); ++it)
	{
		int pid = it->first;
		BinaryBuffer &buf = it->second;
		set<AlgebraID> newNodes;
		Deserialize(buf, newNodes);
		bool bCreated;
		for(set<AlgebraID>::iterator it = newNodes.begin(); it != newNodes.end(); ++it)
		{
			AlgebraID &globalID = (*it);
			int localIndex = get_local_index_or_create_new(globalID, distanceToMasterOrInner, bCreated);
			if(bCreated)
				newSlaves[pid].insert(localIndex);
		}
	}*/

	insert_into_layout_sorted(newMasters, totalMasterLayout);
	insert_into_layout_sorted(newSlaves, totalSlaveLayout);
	for(auto it = newMasters.begin(); it != newMasters.end(); ++it)
	{
		/*UG_LOG("new masters on " << it->first << ": ");
		for(set<size_t>::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2)
			UG_LOG(*it2 << " ");
		UG_LOG("\n");*/
		masterPIDs.insert(it->first);
	}
	for(auto it = newSlaves.begin(); it != newSlaves.end(); ++it)
	{
		/*UG_LOG("new slaves on " << it->first << ": ");
		for(set<size_t>::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2)
			UG_LOG(*it2 << " ");
		UG_LOG("\n");*/
		slavePIDs.insert(it->first);
	}
}

void ParallelNodes::add_new_layouts_to(IndexLayout &newMasterLayout, IndexLayout &newSlaveLayout)
{
	insert_into_layout_sorted(newMasters, newMasterLayout);
	insert_into_layout_sorted(newSlaves, newSlaveLayout);

	newMasters.clear();
	newSlaves.clear();
}

#ifdef NAE_APPEND
/// Appending masters (without comm skills)
void ParallelNodes::append_nodes_without_comm(size_t nnodes)
{

	const size_t nold = m_localToGlobal.size();
	const size_t nnew = nold + nnodes;

	const int myRank = pcl::ProcRank();
	UG_ASSERT(myRank==0, "TODO: Implement for real parallel case")

	m_localToGlobal.resize(nnew);
	for(size_t i=nold; i<nnew; ++i)
	{
		m_localToGlobal[i] = AlgebraID(myRank, i);
		if(m_localToGlobal[i].first != myRank)
					m_globalToLocal.insert(pair<AlgebraID, size_t> (m_localToGlobal[i], i));
	}

}
#endif
} // namespace ug
