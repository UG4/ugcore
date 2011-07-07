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

bool SubsetGroup::add(int si)
{
	if(!is_init())
	{
		UG_LOG("No SubsetHandler set. "
				"Cannot use SubsetGroup without SubsetHandler.\n");
		return false;
	}

	if(si < 0)
		{UG_LOG("Subset indices must be non-negative.\n");return false;}

	std::vector<int>::iterator iter;
	iter = find(m_vSubset.begin(), m_vSubset.end(), si);
	if(iter != m_vSubset.end()) return true;

	m_vSubset.push_back(si);
	sort(m_vSubset.begin(), m_vSubset.end());
	return true;
}

bool SubsetGroup::add(const char* name)
{
	if(!is_init())
	{
		UG_LOG("No SubsetHandler set. "
				"Cannot use SubsetGroup without SubsetHandler.\n");
		return false;
	}

	size_t found = 0;

//	Search for name in Subset list
	for(int si = 0; si < m_pSH->num_subsets(); ++si)
	{
		if(strcmp (name, m_pSH->get_subset_name(si)) == 0)
		{
			found++;
			add(si);
		}
	}

// 	if not found, return false
	if(found == 0) return false;
	return true;
}

bool SubsetGroup::add(const SubsetGroup& ssGroup)
{
	if(!is_init())
	{
		UG_LOG("No SubsetHandler set. "
				"Cannot use SubsetGroup without SubsetHandler.\n");
		return false;
	}

//	check that underlying subset handlers are equal
	if(m_pSH != ssGroup.get_subset_handler())
	{
		UG_LOG("Underlying subset handler does not match. Cannot add"
				" subsets to subset group.\n");
		return false;
	}

//	add all subsets
	for(size_t i = 0; i < ssGroup.num_subsets(); ++i)
		add(ssGroup[i]);

//	we're done
	return true;
}


bool SubsetGroup::add_all()
{
	if(!is_init())
	{
		UG_LOG("No SubsetHandler set. "
				"Cannot use SubsetGroup without SubsetHandler.\n");
		return false;
	}

	for(int si = 0; si < m_pSH->num_subsets(); ++si)
	{
		add(si);
	}

	return true;
}

bool SubsetGroup::remove(int si)
{
	if(!is_init())
	{
		UG_LOG("No SubsetHandler set. "
				"Cannot use SubsetGroup without SubsetHandler.\n");
		return false;
	}

	if(si < 0)
		{UG_LOG("Subset indices must be non-negative.\n");return false;}

	std::vector<int>::iterator iter;
	iter = find(m_vSubset.begin(), m_vSubset.end(), si);
	if(iter == m_vSubset.end())
		{UG_LOG("Index not contained in SubsetGroup.\n");return false;}

	m_vSubset.erase(iter);
	return true;
}

bool SubsetGroup::remove(const char* name)
{
	if(!is_init())
	{
		UG_LOG("No SubsetHandler set. "
				"Cannot use SubsetGroup without SubsetHandler.\n");
		return false;
	}

	size_t found = 0;

//	Search for name in Subset list
	for(int si = 0; si < m_pSH->num_subsets(); ++si)
	{
		if(strcmp (name, m_pSH->get_subset_name(si)) == 0)
		{
			found++;
			remove(si);
		}
	}

// 	if not found, return false
	if(found == 0) return false;
	return true;
}

///	name of subset
const char* SubsetGroup::name(size_t i) const
{
	if(!is_init())
	{
		UG_LOG("No SubsetHandler set. "
				"Cannot use SubsetGroup without SubsetHandler.\n");
		return false;
	}

//	Check, that subset exist
	if(i >= num_subsets())
		throw(ERROR_BadIndexInSubsetGroup(i));

	return m_pSH->get_subset_name(m_vSubset[i]);
}

///	returns if a subset is a regular grid
bool SubsetGroup::regular_grid(size_t i) const
{
	if(!is_init())
	{
		UG_LOG("No SubsetHandler set. "
				"Cannot use SubsetGroup without SubsetHandler.\n");
		return false;
	}

//	Check, that subset exist
	if(i >= num_subsets())
		throw(ERROR_BadIndexInSubsetGroup(i));

	return SubsetIsRegularGrid(*m_pSH, m_vSubset[i]);
}


int SubsetGroup::dim(size_t i
#ifdef UG_PARALLEL
                     ,pcl::ProcessCommunicator* pProcCom
#endif
                     ) const
{
	if(!is_init())
	{
		UG_LOG("No SubsetHandler set. "
				"Cannot use SubsetGroup without SubsetHandler.\n");
		return false;
	}

//	Check, that subset exist
	if(i >= num_subsets())
		throw(ERROR_BadIndexInSubsetGroup(i));

#ifdef UG_PARALLEL
	return DimensionOfSubset(*m_pSH, m_vSubset[i], pProcCom);
#else
	return DimensionOfSubset(*m_pSH, m_vSubset[i]);
#endif
}

int SubsetGroup::get_local_highest_subset_dimension() const
{
	if(!is_init())
	{
		UG_LOG("No SubsetHandler set. "
				"Cannot use SubsetGroup without SubsetHandler.\n");
		return false;
	}

//	without subsets no dimension
	if(num_subsets() == 0) return -1;

//	get first dimension
	int d = dim(0);

//	loop other dimensions and compare
	for(size_t i = 0; i < num_subsets(); ++i)
	{
		int test_dim = dim(i);
		if(d < test_dim)
			d = test_dim;
	}

//	returnn dimension
	return d;
}

bool SubsetGroup::contains(int si) const
{
	if(!is_init())
	{
		UG_LOG("No SubsetHandler set. "
				"Cannot use SubsetGroup without SubsetHandler.\n");
		return false;
	}

	//TODO: since we subetindices are sorted be increasing number, one could optimize the search
	std::vector<int>::const_iterator iter;
	iter = find(m_vSubset.begin(), m_vSubset.end(), si);
	if(iter == m_vSubset.end()) return false;
	return true;
}

bool SubsetGroup::contains(const char* name) const
{
	if(!is_init())
	{
		UG_LOG("No SubsetHandler set. "
				"Cannot use SubsetGroup without SubsetHandler.\n");
		return false;
	}

//	Search for name in Subset list
	for(int si = 0; si < m_pSH->num_subsets(); ++si)
	{
		if(strcmp (name, m_pSH->get_subset_name(si)) == 0)
			return true;
	}

	return false;
}

} // end namespace ug
