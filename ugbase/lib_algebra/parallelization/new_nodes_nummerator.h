/*
 * new_node_nummerator.h
 *
 *  Created on: 21.01.2010
 *      Author: Martin Rupp
 */

#ifndef __H__LIB_ALGEBRA__PARALLELIZATION__NEW_NODE_NUMMERATOR__
#define __H__LIB_ALGEBRA__PARALLELIZATION__NEW_NODE_NUMMERATOR__

#include <map>
#include "parallelization_util.h"

namespace ug
{

///\ingroup lib_algebra_parallelization

///\brief Helper class to assign
/**
 * This class provides methods which will get you indices for
 * global indices. These can be new indices or old local indices
 *
 */
class NewNodesNummerator
{
private:
	std::map<AlgebraID, size_t> m_indicesMap;
	typedef std::map<AlgebraID,size_t>::iterator iterator;
	typedef std::map<AlgebraID,size_t>::const_iterator const_iterator;
	iterator m_it;

	std::vector<AlgebraID> &m_globalIDs;

public:
	/// constructor. new_indices_start is the index of the first newly create local index
	NewNodesNummerator(std::vector<AlgebraID> &globalIDs) : m_globalIDs(globalIDs)
	{
		for(size_t i=0; i<globalIDs.size(); ++i)
		{
			if(globalIDs[i].first != pcl::GetProcRank() || globalIDs[i].second != i)
				m_indicesMap.insert(std::pair<AlgebraID, size_t> (globalIDs[i], i));
		}
	}


	/// get_index_or_create_new: returns a local index by creating and saving a new one or returning an old
	size_t get_index_or_create_new(const AlgebraID &globalIndex, bool &bCreated)
	{
		// UG_LOG("get_index_or_create_new for global index " << globalIndex.first << " | " << globalIndex.second << ": ");
		if(globalIndex.first == pcl::GetProcRank())
		{
			// UG_LOG("is on this processor\n");
			bCreated = false;
			return globalIndex.second;
		}
		else
		{
			std::pair<iterator, bool> ret = m_indicesMap.insert(std::pair<AlgebraID, size_t> (globalIndex, m_globalIDs.size()));
			if(ret.second)
			{
				// UG_LOG("created new index " << m_globalIDs.size() << "\n");
				m_globalIDs.push_back(globalIndex);
				bCreated = true;
			}
			else
			{
				// UG_LOG("already has index " << ret.first->second << "\n");
				bCreated = false;
			}
			return ret.first->second;
		}
	}

	size_t get_index_or_create_new(const AlgebraID &globalIndex)
	{
		bool b;
		return get_index_or_create_new(globalIndex, b);
	}

	size_t operator [] (const AlgebraID &globalIndex) const
	{
		bool hasIndex;
		size_t index = get_index_if_available(globalIndex, hasIndex);
		UG_ASSERT(hasIndex, "global id " << globalIndex << " has no associated local id");
		return index;
	}

	/// returns a local index by returning a old local one or a saved created one
	size_t get_index_if_available(const AlgebraID &globalIndex, bool &bHasIndex) const
	{
		// UG_LOG("get_index_if_available for global index " << globalIndex << ": ");
		if(globalIndex.first == pcl::GetProcRank())
		{
			// UG_LOG("is on this processor\n");
			bHasIndex = true;
			return globalIndex.second;
		}
		const_iterator it = m_indicesMap.find(globalIndex);
		if(it == m_indicesMap.end())
		{
			// UG_LOG("index not found\n");
			bHasIndex = false;
			return (size_t) -1;
		}
		else
		{
			// UG_LOG("index is " << it->second << "\n");
			bHasIndex = true;
			return it->second;
		}
	}

	/// returns the index of the last new created local index+1
	size_t get_new_indices_size() const
	{
		return m_globalIDs.size();
	}

};


}
#endif
