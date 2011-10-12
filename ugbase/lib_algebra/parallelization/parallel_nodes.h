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


///	serializes data from a vector<bool> into a binary stream
template <class TOStream>
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
}

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
	ParallelNodes(pcl::ParallelCommunicator<IndexLayout> communicator, IndexLayout &masterLayout, IndexLayout &slaveLayout, size_t s)
	: m_masterLayout(masterLayout), m_slaveLayout(slaveLayout)
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

	void sort_layouts(IndexLayout &layout)
	{
		for(IndexLayout::iterator iter = layout.begin(); iter != layout.end(); ++iter)
			sort_interface(layout.interface(iter));
	}


	std::map<AlgebraID, size_t> m_globalToLocal;
	std::vector<AlgebraID> m_localToGlobal;
	IndexLayout &m_masterLayout;
	IndexLayout &m_slaveLayout;

	std::vector<IndexLayout> OLMasterLayouts;
	std::vector<IndexLayout> OLSlaveLayouts;


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

	void print()
	{
		for(size_t i=0; i<local_size(); i++)
			UG_LOG(i << ": " << overlap_type(i) << " globalID = " << local_to_global(i) << "\n");

	}

	const OverlapType &overlap_type(size_t i)
	{
		return m_OLtype[i];
	}

	std::vector<OverlapType> m_OLtype;
};


}
#endif

