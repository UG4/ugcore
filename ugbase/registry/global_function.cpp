/*
 * global_function.cpp
 *
 *  Created on: 11.07.2011
 *      Author: andreasvogel
 */

#include "global_function.h"
#include <iostream>

namespace ug
{
namespace bridge
{

ExportedFunctionBase::
ExportedFunctionBase(	const char* funcName, const char* funcOptions,
						const char* retValInfos, const char* paramInfos,
						const char* tooltip, const char* help)
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
	{
		tokenize(vParamInfoTmp[i], m_vvParamInfo[i], '|');
	}
};

bool ExportedFunctionBase::check_consistency(const char *classname) const
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
				std::cout<<"#### Registry ERROR: Unregistered Class used in ";
				if(classname) std::cout<<"Method: '";
				else std::cout<<"global Function: '";
				PrintFunctionInfo(*this, false, classname);
				std::cout<<"': Parameter " << j+1;
			}
			else
			{	std::cout<<", " << j+1;	}
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
				std::cout<<"#### Registry ERROR: Unregistered Class used in ";
				if(classname) std::cout<<"Method: '";
				else std::cout<<"global Function: '";
				PrintFunctionInfo(*this, false, classname);
				std::cout<<"': Return value ";
			}
			else
			{	std::cout<<", Return value ";	}
		}
	}

//	check if undeclared parameter has been found
	if(bUndeclaredParameterFound) {std::cout<<"\n"; return false;}

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
		tokens.push_back(trim(token));
	}
}

std::string ExportedFunctionBase::trim(const std::string& str)
{
	const size_t start = str.find_first_not_of(" \t");
	const size_t end = str.find_last_not_of(" \t");
	if(start == std::string::npos || end == std::string::npos) return "";
	return str.substr(start, end - start + 1);
}

} // end namespace bridge
} // end namespace ug
