/*
 * function_group.cpp
 *
 *  Created on: 13.07.2010
 *      Author: andreasvogel
 */

#include "./function_group.h"
#include "common/common.h"
#include <algorithm>
#include <cstring>

namespace ug{

bool FunctionGroup::add(size_t fct)
{
	if(!is_init())
	{
		UG_LOG("No FunctionPattern set. Cannot use "
				"FunctionGroup without FunctionPattern.\n");
		return false;
	}

	if(fct >= m_pFunctionPattern->num_fct())
	{
		UG_LOG("Unique function ID " <<fct << " not conatined in underlying"
		       " function pattern (with num_fct=" <<
		       m_pFunctionPattern->num_fct() << ".\n");
		return false;
	}

	std::vector<size_t>::iterator iter;
	iter = find(m_vFunction.begin(), m_vFunction.end(), fct);
	if(iter != m_vFunction.end()) return true;

	m_vFunction.push_back(fct);
	return true;
}

bool FunctionGroup::add(const char* name)
{
	if(!is_init())
	{
		UG_LOG("No FunctionPattern set. Cannot use "
				"FunctionGroup without FunctionPattern.\n");
		return false;
	}

	size_t found = 0;

//	Search for name in FunctionPattern
	for(size_t fct = 0; fct < m_pFunctionPattern->num_fct(); ++fct)
	{
		if(strcmp (name, m_pFunctionPattern->name(fct)) == 0)
		{
			found++;
			add(fct);
		}
	}

// 	if not found, return false
	if(found == 0)
	{
		UG_LOG("No function with name found in underlying Function Pattern.\n");
		return false;
	}
	return true;
}

bool FunctionGroup::add(const FunctionGroup& fctGroup)
{
	if(!is_init())
	{
		UG_LOG("No FunctionPattern set. Cannot use "
				"FunctionGroup without FunctionPattern.\n");
		return false;
	}

	if(m_pFunctionPattern != fctGroup.get_function_pattern())
	{
		UG_LOG("Underlying function pattern does not match. Cannot"
				" add function group. \n");
		return false;
	}

	for(size_t i = 0; i < fctGroup.num_fct(); ++i)
	{
		add(fctGroup[i]);
	}
	return true;
}

bool FunctionGroup::add_all()
{
	if(!is_init())
	{
		UG_LOG("No FunctionPattern set. Cannot use "
				"FunctionGroup without FunctionPattern.\n");
		return false;
	}

	for(size_t i = 0; i < m_pFunctionPattern->num_fct(); ++i)
		add(i);

	return true;
}

void FunctionGroup::sort()
{
	std::sort(m_vFunction.begin(), m_vFunction.end());
}

bool FunctionGroup::remove(size_t fct)
{
	std::vector<size_t>::iterator iter;
	iter = find(m_vFunction.begin(), m_vFunction.end(), fct);
	if(iter == m_vFunction.end())
		{UG_LOG("Function not contained in FunctionGroup.\n");return false;}

	m_vFunction.erase(iter);
	return true;
}

bool FunctionGroup::remove(const char* name)
{
	if(!is_init())
	{
		UG_LOG("No FunctionPattern set. Cannot use "
				"FunctionGroup without FunctionPattern.\n");
		return false;
	}

	size_t found = 0;

//	Search for name in FunctionPattern
	for(size_t fct = 0; fct < m_pFunctionPattern->num_fct(); ++fct)
	{
		if(strcmp (name, m_pFunctionPattern->name(fct)) == 0)
		{
			found++;
			remove(fct);
		}
	}

// 	if not found, return false
	if(found == 0) return false;
	return true;
}

///	name of function
const char* FunctionGroup::name(size_t i) const
{
	if(!is_init())
	{
		UG_LOG("No FunctionPattern set. Cannot use "
				"FunctionGroup without FunctionPattern.\n");
		return false;
	}

//	Check, that subset exist
	if(i >= num_fct())
		throw(ERROR_BadIndexInFunctionGroup(i));

	return m_pFunctionPattern->name(m_vFunction[i]);
}


/// returns the trial space of the discrete function fct
LocalShapeFunctionSetID FunctionGroup::local_shape_function_set_id(size_t i) const
{
	if(!is_init())
	{
		UG_LOG("No FunctionPattern set. Cannot use "
				"FunctionGroup without FunctionPattern.\n");
		throw(ERROR_BadIndexInFunctionGroup(i));
	}

//	Check, that subset exist
	if(i >= num_fct())
		throw(ERROR_BadIndexInFunctionGroup(i));

	return m_pFunctionPattern->local_shape_function_set_id(m_vFunction[i]);
}

///	dimension of function
int FunctionGroup::dim(size_t i) const
{
	if(!is_init())
	{
		UG_LOG("No FunctionPattern set. Cannot use "
				"FunctionGroup without FunctionPattern.\n");
		return false;
	}

//	Check, that subset exist
	if(i >= num_fct())
		throw(ERROR_BadIndexInFunctionGroup(i));

	return m_pFunctionPattern->dim(m_vFunction[i]);
}

///	dimension of all functions in this group
int FunctionGroup::dim() const
{
	if(!is_init())
	{
		UG_LOG("No FunctionPattern set. Cannot use "
				"FunctionGroup without FunctionPattern.\n");
		return false;
	}

//	without functions no dimension
	if(num_fct() == 0) return -1;

	int d = dim(0);

	for(size_t i = 0; i < num_fct(); ++i)
	{
		int test_dim = dim(i);
		if(d != test_dim)
			return -1;
	}

	return d;
}

bool FunctionGroup::contains(size_t fct) const
{
	std::vector<size_t>::const_iterator iter;
	iter = find(m_vFunction.begin(), m_vFunction.end(), fct);
	if(iter == m_vFunction.end()) return false;
	return true;
}

bool FunctionGroup::contains(const FunctionGroup& fctGroup) const
{
// 	loop all functions
	for(size_t i = 0; i < fctGroup.num_fct(); ++i)
	{
	//	if one function is not contained, return false
		if(!contains(fctGroup[i]))
			return false;
	}

// 	all functions contained
	return true;
}


size_t FunctionGroup::local_index(size_t fct) const
{
	for(size_t i = 0; i < m_vFunction.size(); ++i)
	{
		if(fct == m_vFunction[i]) return i;
	}

	throw(ERROR_BadIndexInFunctionGroup(fct));
	return (size_t)-1;
}

} // end namespace ug
