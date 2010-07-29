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
	std::vector<size_t>::iterator iter;
	iter = find(m_vFunction.begin(), m_vFunction.end(), fct);
	if(iter != m_vFunction.end())
		 {UG_LOG("Function already contained in FunctionGroup.\n");return false;}

	m_vFunction.push_back(fct);
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

bool FunctionGroup::containes_function(size_t fct) const
{
	std::vector<size_t>::const_iterator iter;
	iter = find(m_vFunction.begin(), m_vFunction.end(), fct);
	if(iter == m_vFunction.end()) return false;
	return true;
}

} // end namespace ug
