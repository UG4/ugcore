/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__LIB_ALGEBRA__PARALLELIZATION__NEW_LAYOUT_CREATOR_H_
#define __H__LIB_ALGEBRA__PARALLELIZATION__NEW_LAYOUT_CREATOR_H_

namespace ug {

class NewLayoutCreator
{
private:
	void create_mark_map(IndexLayout &masterLayout)
	{
		PROFILE_FUNC();
		UG_DLOG(LIB_ALG_MATRIX, 4, "\n\nGenerateOverlapClass::CreateMarks\n");
		notified.clear();
		for(IndexLayout::iterator iter = masterLayout.begin(); iter != masterLayout.end(); ++iter)
		{
			IndexLayout::Interface &interface = masterLayout.interface(iter);
			int pid = masterLayout.proc_id(iter);

			std::set<size_t> &mark = notified[pid];
			for(IndexLayout::Interface::iterator iter2 = interface.begin(); iter2 != interface.end(); ++iter2)
			{
				size_t localIndex = interface.get_element(iter2);
				mark.insert(localIndex);
			}
		}
	}
	struct NewSlaveNotification
	{
		NewSlaveNotification() = default;
		NewSlaveNotification(const AlgebraID &_id, int _newSlaveOnPID)
		{
			id = _id;
			newSlaveOnPID = _newSlaveOnPID;
		}
		friend std::ostream &operator << (std::ostream &out, const NewSlaveNotification &n)
		{
			out << "notification that node " << n.id << " is slave on " << n.newSlaveOnPID << " ";
			return out;
		}
		AlgebraID id;
		int newSlaveOnPID;
	};

private:
	using BufferMap = std::map<int, BinaryBuffer>;

	ParallelNodes &PN;
	std::set<int> masterPIDs;
	std::set<int> slavePIDs;

	BufferMap notificationBufferMap;

	/// map for marking nodes
	/**
	 * for each processor we need to have a list which of our master nodes exist on their processor
	 * this is important because we sometimes will need to add them to interfaces
	 *
	 * the map serves two functions:
	 * - knowing which processor has copies of our own master nodes, and so constructing correct interfaces
	 * - knowing which notifications are sent
	 *
	 * - if i is a master node (that is PN.local_to_global(i).master_proc() == pcl::ProcRank())
	 *   then i in notified[pid] means: process pid already knows that i is a slave node on his processor,
	 *   and this processor knows that i the associated master
	 *   that means: i is in a master interface on this processor and in a slave interface on pid.
	 * - if is is not a master node
	 *   then i in notified[pid] means: we already sent a notification to the owner of i (processor PN.local_to_global(i).master_proc())
	 *   that processor pid has a copy of i.
	 *

	 */
	std::map<int, std::set<size_t> > notified; // notified[pid] is the set of indices which are slave on pid.
	std::map<int, std::vector<NewSlaveNotification> > newSlaveNotifications; // newSlaveNotifications is a vector of pids
	std::map<int, std::set<size_t> > newMasters;
	std::map<int, std::set<size_t> > newSlaves;

public:

	NewLayoutCreator(ParallelNodes &_PN, IndexLayout &masterLayout, IndexLayout &slaveLayout)
		: PN(_PN)
	{
		for(auto it = slaveLayout.begin(); it != slaveLayout.end(); ++it)
			slavePIDs.insert(slaveLayout.proc_id(it));
		for(auto it = masterLayout.begin(); it != masterLayout.end(); ++it)
			masterPIDs.insert(masterLayout.proc_id(it));
		create_mark_map(masterLayout);
	}
	size_t create_slave_node(const AlgebraID &globalID, int distanceToMasterOrInner)
	{
		bool bCreated;
		size_t localIndex = PN.get_local_index_or_create_new(globalID, distanceToMasterOrInner, bCreated);

		if(bCreated)
		{
			newSlaves[globalID.master_proc()].insert(localIndex);
			UG_DLOG(LIB_ALG_MATRIX, 4, "created new slave with GID " << globalID << " with local index " << localIndex << " and distanceToMasterOrInner=" << distanceToMasterOrInner << "\n");
		}
		return localIndex;
	}

	//////////////////////////////////////////////////////////////////////////////////////////

	void create_node(const AlgebraID &globalID, int pid)
	{
		create_node(globalID, PN.global_to_local(globalID), pid);
	}

	void create_node(size_t localIndex, int pid)
	{
		create_node(PN.local_to_global(localIndex), localIndex, pid);
	}


	void create_node(const AlgebraID &globalID, size_t localIndex, int pid)
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
		std::set<int> sendingPIDs;
		for(IndexLayout::iterator it = sendingLayout.begin(); it != sendingLayout.end(); ++it)
			sendingPIDs.insert(sendingLayout.proc_id(it));

		std::set<int> receivingPIDs;
		for(IndexLayout::iterator it = receivingLayout.begin(); it != receivingLayout.end(); ++it)
			receivingPIDs.insert(receivingLayout.proc_id(it));
		check_new_nodes(commnunicator, sendingPIDs, receivingPIDs);
	}

	void check_new_nodes(std::vector<int> sendingPIDs, std::vector<int> receivingPIDs)
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

	void issue(pcl::InterfaceCommunicator<IndexLayout> &communicator)
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

	void process()
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

				std::set<size_t> &mark = notified[notification.newSlaveOnPID];
				if(mark.find(notification.id.index_on_master()) == mark.end())
				{
					newMasters[notification.newSlaveOnPID].insert(notification.id.index_on_master());
					mark.insert(notification.id.index_on_master());
				}

			}
		}

		/*for(BufferMap::iterator it = newNodesBufferMap.begin(); it != newNodesBufferMap.end(); ++it)
		{
			int pid = it->first;
			BinaryBuffer &buf = it->second;
			std::set<AlgebraID> newNodes;
			Deserialize(buf, newNodes);
			bool bCreated;
			for(std::set<AlgebraID>::iterator it = newNodes.begin(); it != newNodes.end(); ++it)
			{
				AlgebraID &globalID = (*it);
				int localIndex = PN.get_local_index_or_create_new(globalID, distanceToMasterOrInner, bCreated);
				if(bCreated)
					newSlaves[pid].insert(localIndex);
			}
		}*/
	}

	void add_new_layouts_to(IndexLayout &newMasterLayout, IndexLayout &newSlaveLayout)
	{
		PN.insert_into_layout_sorted(newMasters, newMasterLayout);
		PN.insert_into_layout_sorted(newSlaves, newSlaveLayout);
		for(auto it = newMasters.begin(); it != newMasters.end(); ++it)
			masterPIDs.insert(it->first);
		for(auto it = newSlaves.begin(); it != newSlaves.end(); ++it)
			slavePIDs.insert(it->first);
		newMasters.clear();
		newSlaves.clear();
	}
};

}

#endif