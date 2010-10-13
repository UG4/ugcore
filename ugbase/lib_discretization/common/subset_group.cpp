/*
 * subset_group.cpp
 *
 *  Created on: 13.07.2010
 *      Author: andreasvogel
 */

#include "./subset_group.h"
#include "common/common.h"
#include "lib_discretization/domain_util.h"
#include <algorithm>
#include <cstring>

namespace ug{

bool SubsetGroup::add_subset(int si)
{
	if(!is_init())
		{UG_LOG("No SubsetHandler set. Cannot use SubsetGroup without SubsetHandler.\n");return false;}

	if(si < 0)
		{UG_LOG("Subset indices must be non-negative.\n");return false;}

	std::vector<int>::iterator iter;
	iter = find(m_vSubset.begin(), m_vSubset.end(), si);
	if(iter != m_vSubset.end())
		 {UG_LOG("Index already contained in SubsetGroup.\n");return false;}

	m_vSubset.push_back(si);
	sort(m_vSubset.begin(), m_vSubset.end());
	return true;
}

bool SubsetGroup::add_subset(const char* name)
{
	if(!is_init())
		{UG_LOG("No SubsetHandler set. Cannot use SubsetGroup without SubsetHandler.\n");return false;}

	size_t found = 0;

//	Search for name in Subset list
	for(int si = 0; si < m_pSH->num_subsets(); ++si)
	{
		if(strcmp (name, m_pSH->get_subset_name(si)) == 0)
		{
			found++;
			add_subset(si);
		}
	}

// 	if not found, return false
	if(found == 0) return false;
	return true;
}

void SubsetGroup::add_all_subsets()
{
	for(int si = 0; si < m_pSH->num_subsets(); ++si)
	{
		add_subset(si);
	}
}

bool SubsetGroup::remove_subset(int si)
{
	if(!is_init())
		{UG_LOG("No SubsetHandler set. Cannot use SubsetGroup without SubsetHandler.\n");return false;}

	if(si < 0)
		{UG_LOG("Subset indices must be non-negative.\n");return false;}

	std::vector<int>::iterator iter;
	iter = find(m_vSubset.begin(), m_vSubset.end(), si);
	if(iter == m_vSubset.end())
		{UG_LOG("Index not contained in SubsetGroup.\n");return false;}

	m_vSubset.erase(iter);
	return true;
}

bool SubsetGroup::remove_subset(const char* name)
{
	if(!is_init())
		{UG_LOG("No SubsetHandler set. Cannot use SubsetGroup without SubsetHandler.\n");return false;}

	size_t found = 0;

//	Search for name in Subset list
	for(int si = 0; si < m_pSH->num_subsets(); ++si)
	{
		if(strcmp (name, m_pSH->get_subset_name(si)) == 0)
		{
			found++;
			remove_subset(si);
		}
	}

// 	if not found, return false
	if(found == 0) return false;
	return true;
}

///	name of subset
const char* SubsetGroup::get_subset_name(size_t i) const
{
	if(!is_init())
		{UG_LOG("No SubsetHandler set. Cannot use SubsetGroup without SubsetHandler.\n");return false;}

//	Check, that subset exist
	if(i >= num_subsets())
		throw(ERROR_BadIndexInSubsetGroup(i));

	return m_pSH->get_subset_name(m_vSubset[i]);
}

int SubsetGroup::get_subset_dimension(size_t i) const
{
	if(!is_init())
		{UG_LOG("No SubsetHandler set. Cannot use SubsetGroup without SubsetHandler.\n");return false;}

//	Check, that subset exist
	if(i >= num_subsets())
		throw(ERROR_BadIndexInSubsetGroup(i));

	return DimensionOfSubset(*m_pSH, m_vSubset[i]);
}

int SubsetGroup::get_subset_dimension() const
{
	if(!is_init())
		{UG_LOG("No SubsetHandler set. Cannot use SubsetGroup without SubsetHandler.\n");return false;}

//	without subsets no dimension
	if(num_subsets() == 0) return -1;

	int dim = get_subset_dimension(0);

	for(size_t i = 0; i < num_subsets(); ++i)
	{
		int test_dim = get_subset_dimension(i);
		if(dim != test_dim)
			return -1;
	}

	return dim;
}

bool SubsetGroup::containes_subset(int si) const
{
	if(!is_init())
		{UG_LOG("No SubsetHandler set. Cannot use SubsetGroup without SubsetHandler.\n");return false;}

	//TODO: since we subetindices are sorted be increasing number, one could optimize the search
	std::vector<int>::const_iterator iter;
	iter = find(m_vSubset.begin(), m_vSubset.end(), si);
	if(iter == m_vSubset.end()) return false;
	return true;
}

bool SubsetGroup::containes_subset(const char* name) const
{
	if(!is_init())
		{UG_LOG("No SubsetHandler set. Cannot use SubsetGroup without SubsetHandler.\n");return false;}

//	Search for name in Subset list
	for(int si = 0; si < m_pSH->num_subsets(); ++si)
	{
		if(strcmp (name, m_pSH->get_subset_name(si)) == 0)
			return true;
	}

	return false;
}

} // end namespace ug
