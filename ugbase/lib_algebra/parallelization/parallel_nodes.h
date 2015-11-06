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

#ifndef UG_PARALLEL_NODES_H_
#define UG_PARALLEL_NODES_H_

#include "pcl/pcl.h"
#include <vector>
#include <set>
#include <map>

#include "lib_algebra/parallelization/parallelization_util.h"
#include "lib_algebra/parallelization/algebra_layouts.h"

#include "common/util/sort_util.h"
#include "common/util/binary_buffer.h"
#include "common/serialization.h"


#define NAE_APPEND
//	serializes data from a vector<bool> into a binary stream
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

/**
 * ParallelNodes is a class to provide functions when adding nodes on other processes.
 * Especially it is used in the construction of matrix overlaps,
 * where it can be that process A sends a matrix row to process B containing connections
 * to process C. We need to make sure everyone has the right parallel connections afterwards
 * without the need to have an all-to-all communication
 */
class ParallelNodes
{
private:
	typedef std::map<AlgebraID,size_t>::iterator iterator;
	typedef std::map<AlgebraID,size_t>::const_iterator const_iterator;

public:
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
private:
	/**
	 * these are "i have a new slave to your process" notifications.
	 * the sender sends "i have a new slave to your process", and
	 * the receiver now adds this node to a master layout to this process.
	 */
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
	ConstSmartPtr<AlgebraLayouts> m_layout;

	std::map<AlgebraID, size_t> m_globalToLocal;
	std::vector<AlgebraID> m_localToGlobal;

	IndexLayout totalMasterLayout;
	IndexLayout totalSlaveLayout;

	size_t m_originalSize;
	std::vector<OverlapType> m_OLtype;

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
	ParallelNodes();
	ParallelNodes(const ParallelNodes &);
	ParallelNodes(ConstSmartPtr<AlgebraLayouts> layout, size_t s);

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
		UG_COND_THROW(i >= local_size(), i << " >= " << local_size());
		return m_localToGlobal[i];
	}

	const AlgebraID &operator [] (size_t i) const
	{
		return local_to_global(i);
	}

	/// returns a local index by returning a old local one or a saved created one
	size_t get_local_index_if_available(const AlgebraID &globalIndex, bool &bHasIndex) const;

	/// get_index_or_create_new: returns a local index by creating and saving a new one or returning an old
	size_t get_local_index_or_create_new(const AlgebraID &globalIndex, int distanceToMasterOrInner, bool &bCreated);


	size_t get_local_index_or_create_new(const AlgebraID &globalIndex, int distanceToMasterOrInner);

	size_t global_to_local(const AlgebraID &globalIndex) const;

	size_t operator [] (const AlgebraID &globalIndex) const
	{
		return global_to_local(globalIndex);
	}


	void sort_by_global_id(std::vector<size_t> &v)
	{
		sort(v.begin(), v.end(), CompareIndicesBy(m_localToGlobal) );
	}

	void insert_into_interface_sorted(std::vector<size_t> &v, IndexLayout::Interface &interface);

	void insert_into_layout_sorted(std::map<int, std::set<size_t> > &m, IndexLayout &layout);

	void insert_into_layout_sorted(std::map<int, std::vector<size_t> > &m, IndexLayout &layout);

	void sort_interface(IndexLayout::Interface &interface);

	void sort_layout(IndexLayout &layout);


	IndexLayout &get_total_master_layout()
	{
		return totalMasterLayout;
	}
	IndexLayout &get_total_slave_layout()
	{
		return totalSlaveLayout;
	}

	const IndexLayout &master_layout() const
	{
		return m_layout->master();
	}

	const IndexLayout &slave_layout() const
	{
		return m_layout->slave();
	}

	pcl::InterfaceCommunicator<IndexLayout> &comm() const
	{
		return m_layout->comm();
	}

	const pcl::ProcessCommunicator &proc_comm() const
	{
		return m_layout->proc_comm();
	}

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

	/// returns the overlap type (inner, master, slave or distanceToMasterOrInner=X)
	const OverlapType &overlap_type(size_t i)
	{
		return m_OLtype[i];
	}

	/**
	 * when receiving nodes which get slaves, we call this function to make sure
	 * there is a slave interface to master process.
	 * @param globalID					global ID of possibly new slave
	 * @param distanceToMasterOrInner	distance to nearest master node
	 * @return new local index
	 */
	size_t create_slave_node(const AlgebraID &globalID, int distanceToMasterOrInner);
	//////////////////////////////////////////////////////////////////////////////////////////

	/**
	 * create a new node in this ParallelNode structure on processor pid.
	 * note that the node referred to by globalID has a processor where it is master on.
	 * that is globalID.master(). Now if we create a node on a processor which is no
	 * globalID.master(), we need to inform the process globalID.master() that he's getting
	 * a new node. The process pid will notice by himself that he needs to add globalID into
	 * the slave layout to globalID.master().
	 * This is especially necessary when we are sending matrix rows since a processor A can send
	 * a row to processor B with connections to master nodes on processor C.
	 * @param globalID
	 * @param localIndex
	 * @param pid			where the node is created
	 */
	void create_node(const AlgebraID &globalID, size_t localIndex, int pid);

	/**
	 * @param globalID
	 * @param pid
	 * \sa create_node
	 */
	void create_node(const AlgebraID &globalID, int pid)
	{
		create_node(globalID, global_to_local(globalID), pid);
	}

	/**
	 * @param localIndex
	 * @param pid
	 * \sa create_node
	 */
	void create_node(size_t localIndex, int pid)
	{
		create_node(local_to_global(localIndex), localIndex, pid);
	}


	/**
	 * write all 'i have a new slave to your process' notification into send buffers
	 * issue
	 * - sending of send buffers
	 * - receive of notification
	 * @param communicator
	 */
	void issue(pcl::InterfaceCommunicator<IndexLayout> &communicator);

	/**
	 * call this when communication has been made.
	 * this function processes the received data: here it is
	 * the notification about that another processor has a slave to our process,
	 * so we need to add a node to masterLayout
	 */
	void process();

	/**
	 * PN created some new master and slave nodes. now we add them to another layout
	 * (we communicate the changes to the outside world)
	 * @param newMasterLayout
	 * @param newSlaveLayout
	 */
	void add_new_layouts_to(IndexLayout &newMasterLayout, IndexLayout &newSlaveLayout);
#ifdef NAE_APPEND
	void append_nodes_without_comm(size_t n);
#endif
private:
	void create_mark_map(const IndexLayout &masterLayout);

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


};


}
#endif

