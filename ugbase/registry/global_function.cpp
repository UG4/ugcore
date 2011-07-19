/*
 * global_function.cpp
 *
 *  Created on: 11.07.2011
 *      Author: andreasvogel
 */

#include "global_function.h"
#include <iostream>
#include <algorithm>
#include <string>
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
		UG_LOG("Non-allowed character '"<<m_name[found]<<"' "<<
			   "contained at position "<<int(found)<<" in registered function/method "
			   "'"<<m_name<<"'.\nFunction names must match the regular expression: "
			   "[a-zA-Z_][a-zA-Z_0-9]*, \ni.e. only alphabetic characters, numbers "
			   " and '_' are allowed; no numbers at the beginning.\n");
		throw(REGISTRY_ERROR_Message("Function Name must only contain [a-zA-Z_][a-zA-Z_0-9]*."));
	}

//	check that no number at the beginning
	found = m_name.find_first_of("0123456789");
	if (found!=std::string::npos && found == 0)
	{
		UG_LOG("Function Name "<<m_name<<" starts with a number.\nThis is "
				" not allowed. Please change naming.\n");
		throw(REGISTRY_ERROR_Message("Function Name must not start with number."));
	}
};

bool ExportedFunctionBase::check_consistency(std::string classname) const
{
//	flag to indicate, that unnamed parameter is found
	bool bUndeclaredParameterFound = false;

//	loop method parameters
	for(int j=0; j<params_in().size(); j++)
	{
		if(!params_in().parameter_named(j))
		{
		//	print error output, indicate parameter by 1, ..., NumParams
			if(!bUndeclaredParameterFound)
			{
				bUndeclaredParameterFound = true;
				UG_LOG("#### Registry ERROR: Unregistered Class used in ");
				if(!classname.empty()){ UG_LOG("Method: '");}
				else UG_LOG("global Function: '")
				PrintFunctionInfo(*this, false, classname.c_str());
				UG_LOG("': Parameter " << j+1);
			}
			else
			{	UG_LOG(", " << j+1);	}
		}
	}

//	loop return values
	for(int j=0; j<params_out().size(); j++)
	{
	//	print error output
		if(!params_out().parameter_named(j))
		{
			UG_ASSERT(j == 0, "Not more than one return value can appear.");
			if(!bUndeclaredParameterFound)
			{
				bUndeclaredParameterFound = true;
				UG_LOG("#### Registry ERROR: Unregistered Class used in ");
				if(!classname.empty()){UG_LOG("Method: '");}
				else UG_LOG("global Function: '");
				PrintFunctionInfo(*this, false, classname.c_str());
				UG_LOG("': Return value ");
			}
			else
			{	UG_LOG(", Return value ");	}
		}
	}

//	check if undeclared parameter has been found
	if(bUndeclaredParameterFound) {UG_LOG("\n"); return false;}

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
