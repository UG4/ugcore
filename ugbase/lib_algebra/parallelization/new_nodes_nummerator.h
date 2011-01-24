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
	iterator m_it;
	size_t m_newIndex;

public:
	/// constructor. new_indices_start is the index of the first newly create local index
	NewNodesNummerator(size_t new_indices_start, std::vector<AlgebraID> &global_ids)
	{
		m_newIndex = new_indices_start;
		for(size_t i=0; i<global_ids.size(); ++i)
		{
			if(global_ids[i].first != pcl::GetProcRank() || global_ids[i].second != i)
				m_indicesMap.insert(pair<AlgebraID, size_t> (global_ids[i], i));
		}
	}


	/// get_index_or_create_new: returns a local index by creating and saving a new one or returning an old
	size_t get_index_or_create_new(AlgebraID &global_index)
	{
		// UG_LOG("get_index_or_create_new for global index " << global_index.first << " | " << global_index.second << ": ");
		if(global_index.first == pcl::GetProcRank())
		{
			// UG_LOG("is on this processor\n");
			return global_index.second;
		}
		else
		{
			std::pair<iterator, bool> ret = m_indicesMap.insert(pair<AlgebraID, size_t> (global_index, m_newIndex));
			if(ret.second)
			{
				// UG_LOG("created new index " << m_newIndex << "\n");
				m_newIndex++;
			}
			else
			{
				// UG_LOG("already has index " << m_newIndex << "\n");
			}
			return ret.first->second;
		}
	}

	/// returns a local index by returning a old local one or a saved created one
	size_t get_index_if_available(AlgebraID &global_index, bool &has_index)
	{
		// UG_LOG("get_index_if_available for global index " << global_index.first << " | " << global_index.second << ": ");
		if(global_index.first == pcl::GetProcRank())
		{
			// UG_LOG("is on this processor\n");
			has_index = true;
			return global_index.second;
		}
		iterator it = m_indicesMap.find(global_index);
		if(it == m_indicesMap.end())
		{
			// UG_LOG("index not found\n");
			has_index = false;
			return -1;
		}
		else
		{
			// UG_LOG("index is " << it->second << "\n");
			has_index = true;
			return it->second;
		}
	}

	/// returns the index of the last new created local index+1
	size_t get_new_indices_size()
	{
		return m_newIndex;
	}

};


}
#endif
