/*
 * parallel_nodes.h
 *
 *  Created on: 27.10.2011
 *      Author: mrupp
 */

#ifndef UG_PARALLEL_NODES_H_
#define UG_PARALLEL_NODES_H_

#include "pcl/pcl.h"
#include <vector>
#include <set>
#include <map>

#include "lib_algebra/parallelization/parallelization_util.h"

#include "common/util/sort_util.h"
#include "common/util/binary_buffer.h"
#include "common/serialization.h"


///	serializes data from a vector<bool> into a binary stream
/*template <class TOStream>
void Serialize(TOStream& buf, const std::vector<bool>& vec)
{
	size_t size = vec.size();
	Serialize(buf, size);
	size_t sm8 = size%8;
	size_t s8 = size - sm8;
	for(size_t i=0; i<s8; i+=8)
	{
		char a=0;
		for(size_t j=0; j<8; j++)
			if(vec[i+j])
				a &= (1 << j);
		Serialize(buf, a);
	}
	if(sm8)
	{
		char a=0;
		for(size_t j=0; j<sm8; j++)
			if(vec[s8+j])
				a &= 1 << j;
		Serialize(buf, a);
	}
}

///	deserializes data from a binary stream into a vector<bool>
template <class TIStream>
void Deserialize(TIStream& buf, std::vector<bool> &vec)
{
	size_t size = 0;
	Deserialize(buf, size);
	vec.resize(size);
	size_t sm8 = size%8;
	size_t s8 = size - sm8;
	char a;
	for(size_t i=0; i<s8; i+=8)
	{
		Deserialize(buf, a);
		char a=0;
		for(size_t j=0; j<8; j++)
			vec[i+j] = (a & (1 << j));
	}
	if(sm8)
	{
		Deserialize(buf, a);
		for(size_t j=0; j<sm8; j++)
			vec[s8+j] = (a & (1 << j));
	}
}*/

namespace ug
{

class ParallelNodes
{
private:
	typedef std::map<AlgebraID,size_t>::iterator iterator;
	typedef std::map<AlgebraID,size_t>::const_iterator const_iterator;
public:
	ParallelNodes();
	ParallelNodes(const ParallelNodes &);
	ParallelNodes(pcl::InterfaceCommunicator<IndexLayout> &communicator, IndexLayout &masterLayout, IndexLayout &slaveLayout, size_t s)
	: m_communicator(communicator), m_masterLayout(masterLayout), m_slaveLayout(slaveLayout)
	{
		GenerateGlobalAlgebraIDs(communicator, m_localToGlobal, s, masterLayout, slaveLayout);

		for(size_t i=0; i<m_localToGlobal.size(); ++i)
		{
			if(m_localToGlobal[i].first != pcl::GetProcRank())
				m_globalToLocal.insert(std::pair<AlgebraID, size_t> (m_localToGlobal[i], i));
		}

		m_OLtype.resize(s);
		for(IndexLayout::iterator iter = masterLayout.begin(); iter != masterLayout.end(); ++iter)
		{
			IndexLayout::Interface &interface = masterLayout.interface(iter);
			for(IndexLayout::Interface::iterator iter = interface.begin(); iter != interface.end(); ++iter)
			{
				size_t i = interface.get_element(iter);
				m_OLtype[i].set_master();
			}
		}
		for(IndexLayout::iterator iter = slaveLayout.begin(); iter != slaveLayout.end(); ++iter)
		{
			IndexLayout::Interface &interface = slaveLayout.interface(iter);
			for(IndexLayout::Interface::iterator iter = interface.begin(); iter != interface.end(); ++iter)
			{
				size_t i = interface.get_element(iter);
				UG_ASSERT(m_OLtype[i].is_master() == false, i);
				m_OLtype[i].set_slave();
			}
		}

		AddLayout(totalMasterLayout, masterLayout);
		AddLayout(totalSlaveLayout, slaveLayout);

		for(IndexLayout::iterator it = slaveLayout.begin(); it != slaveLayout.end(); ++it)
			slavePIDs.insert(slaveLayout.proc_id(it));
		for(IndexLayout::iterator it = masterLayout.begin(); it != masterLayout.end(); ++it)
			masterPIDs.insert(masterLayout.proc_id(it));
		create_mark_map(masterLayout);

		m_originalSize = s;
	}

	size_t m_originalSize;
	size_t get_original_size()
	{
		return m_originalSize;
	}

	size_t local_size() const
	{
		return m_localToGlobal.size();
	}

	const AlgebraID &local_to_global(size_t i) const
	{
		UG_ASSERT(i < local_size(), i << " >= " << local_size());
		return m_localToGlobal[i];
	}

	const AlgebraID &operator [] (size_t i) const
	{
		return local_to_global(i);
	}

	/// returns a local index by returning a old local one or a saved created one
	size_t get_local_index_if_available(const AlgebraID &globalIndex, bool &bHasIndex) const
	{
		// UG_DLOG(LIB_ALG_MATRIX, 4, "get_index_if_available for global index " << globalIndex << ": ");
		if(globalIndex.first == pcl::GetProcRank())
		{
			// UG_DLOG(LIB_ALG_MATRIX, 4, "is on this processor\n");
			bHasIndex = true;
			return globalIndex.second;
		}
		const_iterator it = m_globalToLocal.find(globalIndex);
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

	/// get_index_or_create_new: returns a local index by creating and saving a new one or returning an old
	size_t get_local_index_or_create_new(const AlgebraID &globalIndex, int distanceToMasterOrInner, bool &bCreated)
	{
		// UG_DLOG(LIB_ALG_MATRIX, 4, "get_index_or_create_new for global index " << globalIndex.first << " | " << globalIndex.second << ": ");
		if(globalIndex.first == pcl::GetProcRank())
		{
			// UG_DLOG(LIB_ALG_MATRIX, 4, "is on this processor\n");
			bCreated = false;
			return globalIndex.second;
		}
		else
		{
			std::pair<iterator, bool> ret = m_globalToLocal.insert(std::pair<AlgebraID, size_t> (globalIndex, m_localToGlobal.size()));
			if(ret.second)
			{
				UG_DLOG(LIB_ALG_MATRIX, 4, "created new index " << m_localToGlobal.size() << "\n");
				m_localToGlobal.push_back(globalIndex);
				m_OLtype.push_back(OverlapType(distanceToMasterOrInner));
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

	size_t get_local_index_or_create_new(const AlgebraID &globalIndex, int distanceToMasterOrInner)
	{
		bool b;
		return get_local_index_or_create_new(globalIndex, distanceToMasterOrInner, b);
	}

	size_t global_to_local(const AlgebraID &globalIndex) const
	{
		bool hasIndex;
		size_t index = get_local_index_if_available(globalIndex, hasIndex);
		UG_ASSERT(hasIndex, "global id " << globalIndex << " has no associated local id");
		return index;
	}


	size_t operator [] (const AlgebraID &globalIndex) const
	{
		return global_to_local(globalIndex);
	}


	void sort_by_global_id(std::vector<size_t> &v)
	{
		sort(v.begin(), v.end(), CompareIndicesBy(m_localToGlobal) );
	}

	void insert_into_interface_sorted(std::vector<size_t> &v, IndexLayout::Interface &interface)
	{
		sort_by_global_id(v);
		for(size_t i=0; i<v.size(); i++)
			interface.push_back(v[i]);
	}

	void insert_into_layout_sorted(std::map<int, std::set<size_t> > &m, IndexLayout &layout)
	{
		for(std::map<int, std::set<size_t> >::iterator it = m.begin(); it != m.end(); ++it)
		{
			int pid = it->first;
			std::vector<size_t> v;
			v.reserve(it->second.size());
			v.insert(v.end(), it->second.begin(), it->second.end());
			if(v.size() != 0)
				insert_into_interface_sorted(v, layout.interface(pid));
		}
	}

	void insert_into_layout_sorted(std::map<int, std::vector<size_t> > &m, IndexLayout &layout)
	{
		for(std::map<int, std::vector<size_t> >::iterator it = m.begin(); it != m.end(); ++it)
		{
			int pid = it->first;
			std::vector<size_t> &v = it->second;
			if(v.size() != 0)
				insert_into_interface_sorted(v, layout.interface(pid));
		}
	}

	void sort_interface(IndexLayout::Interface &interface)
	{
		std::vector<size_t> v;
		v.reserve(interface.size());
		for(IndexLayout::Interface::iterator iter = interface.begin(); iter != interface.end(); ++iter)
			v.push_back(interface.get_element(iter));

		IndexLayout::Interface interface2;
		insert_into_interface_sorted(v, interface2);
		interface2.swap(interface);
	}

	void sort_layout(IndexLayout &layout)
	{
		for(IndexLayout::iterator iter = layout.begin(); iter != layout.end(); ++iter)
			sort_interface(layout.interface(iter));
	}

	pcl::InterfaceCommunicator<IndexLayout> &m_communicator;

	std::map<AlgebraID, size_t> m_globalToLocal;
	std::vector<AlgebraID> m_localToGlobal;
	IndexLayout &m_masterLayout;
	IndexLayout &m_slaveLayout;

	IndexLayout totalMasterLayout;
	IndexLayout totalSlaveLayout;

	IndexLayout &get_total_master_layout()
	{
		return totalMasterLayout;
	}
	IndexLayout &get_total_slave_layout()
	{
		return totalSlaveLayout;
	}

	IndexLayout &master_layout()
	{
		return m_masterLayout;
	}

	IndexLayout &slave_layout()
	{
		return m_slaveLayout;
	}

	pcl::InterfaceCommunicator<IndexLayout> &communicator()
	{
		return m_communicator;
	}


	struct OverlapType
	{
#define OT_SLAVE_FLAG 4
#define OT_MASTER_FLAG 2
#define OT_INNER_FLAG 1

		enum eOverlapType
		{
			OT_MASTER, OT_SLAVE, OT_OUTER
		};
		OverlapType()
		{
			// create as inner node
			set_inner();
		}

		OverlapType(int distanceToMasterOrInner)
		{
			set_distance_to_master_or_inner(distanceToMasterOrInner);
		}

		bool is_master_or_inner() const
		{
			return is_master() || is_inner();
		}

		bool is_master() const
		{
			return type == OT_MASTER_FLAG;
		}
		bool is_slave() const
		{
			return type == (OT_SLAVE_FLAG & (1 << 3));
		}

		bool is_inner() const
		{
			return type & OT_INNER_FLAG;
		}

		void set_inner()
		{
			type = OT_INNER_FLAG;
		}

		void set_master()
		{
			type = OT_MASTER_FLAG;
		}

		void set_slave()
		{
			type = OT_SLAVE_FLAG & (1 << 3);
		}

		void set_distance_to_master_or_inner(size_t i)
		{
			type = (i << 3);
		}

		size_t distance_to_master_or_inner() const
		{
			return type >> 3;
		}

		friend std::ostream &operator << (std::ostream &out, const OverlapType &o)
		{
			if(o.is_master())
				out << "master ";
			if(o.is_slave())
				out << "slave ";
			if(o.is_inner())
				out << "inner ";
			out << "dTMI=" << o.distance_to_master_or_inner();
			return out;
		}


		int type;
	};

	bool is_master_or_inner(size_t i) const
	{
		return m_OLtype[i].is_master_or_inner();
	}

	bool is_inner(size_t i) const
	{
		return m_OLtype[i].is_inner();
	}
	bool is_slave(size_t i) const
	{
		return m_OLtype[i].is_slave();
	}
	bool is_master(size_t i) const
	{
		return m_OLtype[i].is_master();
	}
	size_t distance_to_master_or_inner(size_t i) const
	{
		return m_OLtype[i].distance_to_master_or_inner();
	}
	void print() const
	{
		for(size_t i=0; i<local_size(); i++)
			UG_LOG(i << ": " << m_OLtype[i] << " globalID = " << local_to_global(i) << "\n");

	}

	const OverlapType &overlap_type(size_t i)
	{
		return m_OLtype[i];
	}

	std::vector<OverlapType> m_OLtype;






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
		NewSlaveNotification() : id() {}
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
	typedef std::map<int, BinaryBuffer>	BufferMap;

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
	 * - if i is a master node (that is PN.local_to_global(i).master_proc() == pcl::GetProcRank())
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


	size_t create_slave_node(const AlgebraID &globalID, int distanceToMasterOrInner)
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

	//////////////////////////////////////////////////////////////////////////////////////////

	void create_node(const AlgebraID &globalID, int pid)
	{
		create_node(globalID, global_to_local(globalID), pid);
	}

	void create_node(size_t localIndex, int pid)
	{
		create_node(local_to_global(localIndex), localIndex, pid);
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

		if(globalID.master_proc() == pcl::GetProcRank())
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
		for(std::set<int>::iterator it = slavePIDs.begin(); it != slavePIDs.end(); ++it)
		{
			BinaryBuffer buf;
			int pid = *it;
			for(size_t i=0; i<newSlaveNotifications[pid].size(); i++)
				Serialize(buf, newSlaveNotifications[pid][i]);
			UG_DLOG(LIB_ALG_MATRIX, 4, " sending " << newSlaveNotifications[pid].size() << " notifications to processor " << pid << " (" << buf.write_pos() << " bytes)\n");
			communicator.send_raw(pid, buf.buffer(), buf.write_pos(), false);
		}
		newSlaveNotifications.clear();

		for(std::set<int>::iterator it = masterPIDs.begin(); it != masterPIDs.end(); ++it)
		{
			int pid = *it;
			UG_DLOG(LIB_ALG_MATRIX, 4, " issue receive from processor " << pid << "\n");
			communicator.receive_raw(pid, notificationBufferMap[pid]);
		}
	}

	void process()
	{
		UG_DLOG(LIB_ALG_MATRIX, 4, "NewLayoutCreator::process\n");
		for(std::set<int>::iterator it = masterPIDs.begin(); it != masterPIDs.end(); ++it)
		{
			int pid = *it;
			BinaryBuffer &buf = notificationBufferMap[pid];
			UG_DLOG(LIB_ALG_MATRIX, 4, "received " << buf.write_pos() << " bytes of notification from pid " << pid << ":\n");
			while(!buf.eof())
			{
				NewSlaveNotification notification;
				Deserialize(buf, notification);
				UG_ASSERT(notification.id.master_proc() == pcl::GetProcRank(), notification.id << ", pid = " << pcl::GetProcRank());
				UG_DLOG(LIB_ALG_MATRIX, 4, notification << "\n");

				std::set<size_t> &mark = notified[notification.newSlaveOnPID];
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
			std::set<AlgebraID> newNodes;
			Deserialize(buf, newNodes);
			bool bCreated;
			for(std::set<AlgebraID>::iterator it = newNodes.begin(); it != newNodes.end(); ++it)
			{
				AlgebraID &globalID = (*it);
				int localIndex = get_local_index_or_create_new(globalID, distanceToMasterOrInner, bCreated);
				if(bCreated)
					newSlaves[pid].insert(localIndex);
			}
		}*/

		insert_into_layout_sorted(newMasters, totalMasterLayout);
		insert_into_layout_sorted(newSlaves, totalSlaveLayout);
		for(std::map<int, std::set<size_t> >::iterator it = newMasters.begin(); it != newMasters.end(); ++it)
		{
			/*UG_LOG("new masters on " << it->first << ": ");
			for(std::set<size_t>::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2)
				UG_LOG(*it2 << " ");
			UG_LOG("\n");*/
			masterPIDs.insert(it->first);
		}
		for(std::map<int, std::set<size_t> >::iterator it = newSlaves.begin(); it != newSlaves.end(); ++it)
		{
			/*UG_LOG("new slaves on " << it->first << ": ");
			for(std::set<size_t>::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2)
				UG_LOG(*it2 << " ");
			UG_LOG("\n");*/
			slavePIDs.insert(it->first);
		}
	}

	void add_new_layouts_to(IndexLayout &newMasterLayout, IndexLayout &newSlaveLayout)
	{
		insert_into_layout_sorted(newMasters, newMasterLayout);
		insert_into_layout_sorted(newSlaves, newSlaveLayout);

		newMasters.clear();
		newSlaves.clear();
	}
};


}
#endif

