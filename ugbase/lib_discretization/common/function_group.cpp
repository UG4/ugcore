/*
 * function_group.cpp
 *
 *  Created on: 13.07.2010
 *      Author: andreasvogel
 */

#include "./function_group.h"
#include "common/common.h"
#include <algorithm>

namespace ug{

bool FunctionGroup::add_function(size_t fct)
{
	if(!is_init())
	{
		UG_LOG("No FunctionPattern set. Cannot use "
				"FunctionGroup without FunctionPattern.\n");
		return false;
	}

	std::vector<size_t>::iterator iter;
	iter = find(m_vFunction.begin(), m_vFunction.end(), fct);
	if(iter != m_vFunction.end())
		 {UG_LOG("Function already contained in FunctionGroup.\n");return false;}

	m_vFunction.push_back(fct);
	return true;
}

bool FunctionGroup::add_function(const char* name)
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
			add_function(fct);
		}
	}

// 	if not found, return false
	if(found == 0) return false;
	return true;
}

bool FunctionGroup::remove_function(size_t fct)
{
	std::vector<size_t>::iterator iter;
	iter = find(m_vFunction.begin(), m_vFunction.end(), fct);
	if(iter == m_vFunction.end())
		{UG_LOG("Function not contained in FunctionGroup.\n");return false;}

	m_vFunction.erase(iter);
	return true;
}

bool FunctionGroup::remove_function(const char* name)
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
			remove_function(fct);
		}
	}

// 	if not found, return false
	if(found == 0) return false;
	return true;
}

///	name of function
const char* FunctionGroup::get_function_name(size_t i) const
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
int FunctionGroup::get_function_dimension(size_t i) const
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
int FunctionGroup::get_function_dimension() const
{
	if(!is_init())
	{
		UG_LOG("No FunctionPattern set. Cannot use "
				"FunctionGroup without FunctionPattern.\n");
		return false;
	}

//	without functions no dimension
	if(num_fct() == 0) return -1;

	int dim = get_function_dimension(0);

	for(size_t i = 0; i < num_fct(); ++i)
	{
		int test_dim = get_function_dimension(i);
		if(dim != test_dim)
			return -1;
	}

	return dim;
}

bool FunctionGroup::containes_function(size_t fct) const
{
	std::vector<size_t>::const_iterator iter;
	iter = find(m_vFunction.begin(), m_vFunction.end(), fct);
	if(iter == m_vFunction.end()) return false;
	return true;
}

} // end namespace ug
