/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
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

#include "global_function.h"
#include <iostream>
#include <algorithm>
#include <string>
#include "class_helper.h"
#include "class_name_provider.h"
#include "common/util/string_util.h"

namespace ug
{
namespace bridge
{

ExportedFunctionBase::
ExportedFunctionBase(	const std::string& funcName, const std::string& funcOptions,
						const std::string& retValInfos, const std::string& paramInfos,
						const std::string& tooltip, const std::string& help)
: m_name(funcName), m_methodOptions(funcOptions),
  m_retValInfos(retValInfos), m_paramInfos(paramInfos),
  m_tooltip(tooltip), m_help(help)
{
//	Tokenize string for return value (separated by '|')
	tokenize(m_retValInfos, m_vRetValInfo, '|');

//	Tokenize string for parameters into infos per one parameter (separated by '#')
	std::vector<std::string> vParamInfoTmp;
	tokenize(m_paramInfos, vParamInfoTmp, '#');
	m_vvParamInfo.resize(vParamInfoTmp.size());

//	Tokenize each info-string of one parameter into single infos (separated by '|')
	for(size_t i = 0; i < vParamInfoTmp.size(); ++i)
		tokenize(vParamInfoTmp[i], m_vvParamInfo[i], '|');

//	check name   //
///////////////////

//	check for non-allowed character
	size_t found = m_name.find_first_not_of("abcdefghijklmnopqrstuvwxyz"
											 "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
											 "_"
											 "0123456789");
	if (found!=std::string::npos)
	{
		UG_ERR_LOG("Non-allowed character '"<<m_name[found]<<"' "<<
			   "contained at position "<<int(found)<<" in registered function/method "
			   "'"<<m_name<<"'.\nFunction names must match the regular expression: "
			   "[a-zA-Z_][a-zA-Z_0-9]*, \ni.e. only alphabetic characters, numbers "
			   " and '_' are allowed; no numbers at the beginning.\n");
		UG_THROW_REGISTRY_ERROR(m_name,"Function Name must only contain [a-zA-Z_][a-zA-Z_0-9]*.");
	}

//	check that no number at the beginning
	found = m_name.find_first_of("0123456789");
	if (found!=std::string::npos && found == 0)
	{
		UG_ERR_LOG("Function Name "<<m_name<<" starts with a number.\nThis is "
				" not allowed. Please change naming.\n");
		UG_THROW_REGISTRY_ERROR(m_name, "Function Name must not start with number.");
	}
};

bool ExportedFunctionBase::check_consistency(const std::string &classname) const
{
//	flag to indicate, that unnamed parameter is found
	bool bUndeclared = false, bUndeclaredParameter = false, bUndeclaredReturn = false;

//	loop method parameters
	for(int j=0; j<params_in().size(); j++)
		if(!params_in().parameter_named(j))
		{
			bUndeclared = true;
			bUndeclaredParameter = true;
		}

//	loop return values
	for(int j=0; j<params_out().size(); j++)
		if(!params_out().parameter_named(j))
		{
			bUndeclared = true;
			bUndeclaredReturn = true;
		}

//	print error message
	if(bUndeclared)
	{
		UG_ERR_LOG("#### Registry ERROR: Unregistered Class used in ");
		if(!classname.empty()){ UG_ERR_LOG("Method: '");}
		else UG_ERR_LOG("global Function: '")

		if(params_out().size() > 0){
			UG_ERR_LOG(ParameterToString(params_out(), 0) << " ");
		}
		else { UG_ERR_LOG("void ");}
		if(!classname.empty()) UG_ERR_LOG(classname << ":");
		UG_ERR_LOG(name() << "(");
		const size_t ssize = static_cast<size_t>(params_in().size());
		for(size_t i=0; i < ssize; ++i)
		{
			if(i>0) UG_ERR_LOG(", ");
			UG_ERR_LOG(ParameterToString(params_in(), i));
			if(i < num_parameter()) UG_ERR_LOG(" " << parameter_name(i));
		}
		UG_ERR_LOG(")':");
	}

	bool bNext = false;
	if(bUndeclaredParameter)
	{
		UG_ERR_LOG(" for Parameter ");
		for(int j=0; j<params_in().size(); j++)
			if(!params_in().parameter_named(j))
			{
				if(bNext) UG_ERR_LOG(", ")
				UG_ERR_LOG(j+1);
				bNext = true;
			}
	}

	if(bNext) UG_ERR_LOG(", ")
	if(bUndeclaredReturn)
		for(int j=0; j<params_out().size(); j++)
			if(!params_out().parameter_named(j))
				UG_ERR_LOG("for Return value ");

//	check if undeclared parameter has been found
	if(bUndeclared) {UG_ERR_LOG("\n"); return false;}

//	everything ok
	return true;
}

void ExportedFunctionBase::tokenize(const std::string& str,
                                    std::vector<std::string>& tokens,
                                    const char delimiter)
{
	tokens.clear();
	std::stringstream tokenstream;
	tokenstream << str;
	std::string token;

	while ( std::getline (tokenstream, token, delimiter ) )
	{
		tokens.push_back(TrimString(token));
	}
}

} // end namespace bridge
} // end namespace ug
