/*
 * subset_index_group.cpp
 *
 *  Created on: 13.07.2010
 *      Author: andreasvogel
 */

#include "./subset_index_group.h"
#include "common/common.h"
#include <algorithm>

namespace ug{

bool SubsetIndexGroup::add_subset(int si)
{
	if(si < 0)
		{UG_LOG("Subset indices must be non-negative.\n");return false;}

	std::vector<int>::iterator iter;
	iter = find(m_vSubset.begin(), m_vSubset.end(), si);
	if(iter != m_vSubset.end())
		 {UG_LOG("Index already contained in SubsetIndexGroup.\n");return false;}

	m_vSubset.push_back(si);
	sort(m_vSubset.begin(), m_vSubset.end());
	return true;
}

bool SubsetIndexGroup::remove_subset(int si)
{
	if(si < 0)
		{UG_LOG("Subset indices must be non-negative.\n");return false;}

	std::vector<int>::iterator iter;
	iter = find(m_vSubset.begin(), m_vSubset.end(), si);
	if(iter == m_vSubset.end())
		{UG_LOG("Index not contained in SubsetIndexGroup.\n");return false;}

	m_vSubset.erase(iter);
	return true;
}

bool SubsetIndexGroup::containes_subset(int si) const
{
	//TODO: since we subetindices are sorted be increasing number, one could optimize the search
	std::vector<int>::const_iterator iter;
	iter = find(m_vSubset.begin(), m_vSubset.end(), si);
	if(iter == m_vSubset.end()) return false;
	return true;
}

} // end namespace ug
