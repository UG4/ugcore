/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
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

#include "function_group.h"
#include "common/common.h"
#include "common/util/string_util.h"
#include <algorithm>
#include <cstring>

using namespace std;

namespace ug{

FunctionGroup::FunctionGroup() : m_spFunctionPattern(nullptr) {clear();}

FunctionGroup::FunctionGroup(ConstSmartPtr<FunctionPattern> spFuncPattern)
	: m_spFunctionPattern(spFuncPattern)
{
	clear();
}

FunctionGroup::FunctionGroup(ConstSmartPtr<FunctionPattern> spFuncPattern, const char* name)
	: m_spFunctionPattern(spFuncPattern)
{
	clear();
	add(name);
}

FunctionGroup::FunctionGroup(ConstSmartPtr<FunctionPattern> spFuncPattern, const std::string& name)
	: m_spFunctionPattern(spFuncPattern)
{
	clear();
	add(name);
}

FunctionGroup::FunctionGroup(ConstSmartPtr<FunctionPattern> spFuncPattern, const std::vector<std::string>& vName)
	: m_spFunctionPattern(spFuncPattern)
{
	clear();
	add(vName);
}

void FunctionGroup::set_function_pattern(ConstSmartPtr<FunctionPattern> spFuncPattern)
{
	m_spFunctionPattern = spFuncPattern;
	clear();
}


void FunctionGroup::add(size_t fct)
{
	if(!is_init())
		UG_THROW("Cannot use FunctionGroup without FunctionPattern.");

	if(fct >= m_spFunctionPattern->num_fct())
		UG_THROW("Unique function ID " <<fct << " not contained in "
		               "underlying function pattern (with size=" <<
		               	  m_spFunctionPattern->num_fct() << ".");

	auto iter = find(m_vFunction.begin(), m_vFunction.end(), fct);
	if(iter != m_vFunction.end()) return;

	m_vFunction.push_back(fct);
}

void FunctionGroup::add(const string& name)
{
	string tName = TrimString(name);

	if(!is_init())
		UG_THROW("Cannot use FunctionGroup without FunctionPattern.");

	size_t found = 0;

//	Search for name in FunctionPattern
	for(size_t fct = 0; fct < m_spFunctionPattern->num_fct(); ++fct)
	{
		if(strcmp (tName.c_str(), m_spFunctionPattern->name(fct)) == 0)
		{
			found++;
			add(fct);
		}
	}

// 	if not found, return false
	if(found == 0)
		UG_THROW("FunctionGroup: no function '"<<tName<<"' in Function Pattern.");
}

void FunctionGroup::add(const char* name)
{
	add(TokenizeString(name));
}

void FunctionGroup::add(const vector<string>& vName)
{
	for (size_t i = 0; i < vName.size(); ++i) {
		add(vName[i]);
	}
}

void FunctionGroup::add(const FunctionGroup& fctGroup)
{
	if(!is_init())
		UG_THROW("Cannot use FunctionGroup without FunctionPattern.");

	if(m_spFunctionPattern != fctGroup.function_pattern())
		UG_THROW("Underlying function pattern does not match. Cannot"
				" add function group.");

	for(size_t i = 0; i < fctGroup.size(); ++i)
		add(fctGroup[i]);
}

void FunctionGroup::add_all()
{
	if(!is_init())
		UG_THROW("Cannot use FunctionGroup without FunctionPattern.");

	for(size_t i = 0; i < m_spFunctionPattern->num_fct(); ++i)
		add(i);
}

void FunctionGroup::sort()
{
	std::sort(m_vFunction.begin(), m_vFunction.end());
}

void FunctionGroup::remove(size_t fct)
{
	auto iter = find(m_vFunction.begin(), m_vFunction.end(), fct);
	if(iter == m_vFunction.end())
		UG_THROW("Function "<<fct<<" not contained in FunctionGroup.");

	m_vFunction.erase(iter);
}

void FunctionGroup::remove(const string& name)
{
	string tName  = TrimString(name);

	if(!is_init())
		UG_THROW("Cannot use FunctionGroup without FunctionPattern.");

	size_t found = 0;

//	Search for name in FunctionPattern
	for(size_t fct = 0; fct < m_spFunctionPattern->num_fct(); ++fct)
	{
		if(strcmp (tName.c_str(), m_spFunctionPattern->name(fct)) == 0)
		{
			found++;
			remove(fct);
		}
	}

// 	if not found, return false
	if(found == 0)
		UG_THROW("FunctionGroup: no function '"<<tName<<"' in Function Pattern.");
}

void FunctionGroup::remove(const char* name)
{
	remove(TokenizeString(name));
}

void FunctionGroup::remove(const vector<string>& vName)
{
	for (size_t i = 0; i < vName.size(); ++i) {
		remove(vName[i]);
	}
}

const char* FunctionGroup::name(size_t i) const
{
	if(!is_init())
		UG_THROW("Cannot use FunctionGroup without FunctionPattern.");

//	Check, that subset exist
	if(i >= size())
		UG_THROW("Function index "<<i<<" not valid.");

	return m_spFunctionPattern->name(m_vFunction[i]);
}

std::string FunctionGroup::names() const
{
	if(!is_init())
		UG_THROW("Cannot use FunctionGroup without FunctionPattern.");

	std::string s;
	for(size_t i = 0; i < size(); ++i){
		if(i > 0) s.append(", ");
		s.append(name(i));
	}
	return s;
}


LFEID FunctionGroup::local_finite_element_id(size_t i) const
{
	if(!is_init())
		UG_THROW("Cannot use FunctionGroup without FunctionPattern.");

//	Check, that subset exist
	if(i >= size())
		UG_THROW("Function index "<<i<<" not valid.");

	return m_spFunctionPattern->local_finite_element_id(m_vFunction[i]);
}

LFEID FunctionGroup::lfeid(size_t i) const{
	return local_finite_element_id(i);
}

int FunctionGroup::dim(size_t i) const
{
	if(!is_init())
		UG_THROW("Cannot use FunctionGroup without FunctionPattern.");

//	Check, that subset exist
	if(i >= size())
		UG_THROW("Function index "<<i<<" not valid.");

	return m_spFunctionPattern->dim(m_vFunction[i]);
}

int FunctionGroup::dim() const
{
	if(!is_init())
		UG_THROW("Cannot use FunctionGroup without FunctionPattern.");

//	without functions no dimension
	if(size() == 0) return -1;

	int d = dim(0);

	for(size_t i = 0; i < size(); ++i)
	{
		int test_dim = dim(i);
		if(d != test_dim)
			return -1;
	}

	return d;
}

bool FunctionGroup::contains(size_t fct) const
{
	auto iter = find(m_vFunction.begin(), m_vFunction.end(), fct);
	if(iter == m_vFunction.end()) return false;
	return true;
}

bool FunctionGroup::contains(const FunctionGroup& fctGroup) const
{
// 	loop all functions
	for(size_t i = 0; i < fctGroup.size(); ++i)
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

	UG_THROW("Function index "<<fct<<" not valid.");
}

} // end namespace ug
