/*
 * class.cpp
 *
 *  Created on: 26.09.2011
 *      Author: andreasvogel
 */

#include "class.h"

namespace ug{
namespace bridge{

////////////////////////////////////////////////////////////////////////////////
// ExportedConstructor
////////////////////////////////////////////////////////////////////////////////

ExportedConstructor::
ExportedConstructor(ProxyFunc pf,
                    const std::string& className, const std::string& options,
                    const std::string& paramInfos,
                    const std::string& tooltip, const std::string& help)
: m_proxy_func(pf), m_className(className),
  m_options(options), m_paramInfos(paramInfos), m_tooltip(tooltip), m_help(help)
{
//	Tokenize string for parameters into infos per one parameter (separated by '#')
	std::vector<std::string> vParamInfoTmp;
	tokenize(m_paramInfos, vParamInfoTmp, '#');
	m_vvParamInfo.resize(vParamInfoTmp.size());

//	Tokenize each info-string of one parameter into single infos (separated by '|')
	for(size_t i = 0; i < vParamInfoTmp.size(); ++i)
		tokenize(vParamInfoTmp[i], m_vvParamInfo[i], '|');
};

bool ExportedConstructor::check_consistency(std::string classname) const
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
				UG_LOG("Constructor of class "<<classname.c_str());
				UG_LOG("': Parameter " << j+1);
			}
			else
			{	UG_LOG(", " << j+1);	}
		}
	}

//	check if undeclared parameter has been found
	if(bUndeclaredParameterFound) {UG_LOG("\n"); return false;}

//	everything ok
	return true;
}

void ExportedConstructor::tokenize(const std::string& str,
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

////////////////////////////////////////////////////////////////////////////////
// Interface Exported Class
////////////////////////////////////////////////////////////////////////////////

bool IExportedClass::check_consistency() const
{
//	get class name vector of all parents
	const std::vector<const char*>* vClassNames = class_names();

//	check if class name vector correct
	if(vClassNames==NULL)
	{
		UG_LOG("ERROR in 'IExportedClass::check_consistency':"
				" Class name vector of parent classes missing for "
				"class '"<<this->name()<<"'.\n");
		return false;
	}

//	loop all base classes
	for(size_t i = 0; i < (*vClassNames).size(); ++i)
	{
	//	get name of base class
		const char* baseName = (*vClassNames)[i];

	//	check the name
		if(baseName == NULL || strlen(baseName) == 0 || baseName[0] == '[')
		{
			if(i>0){
			UG_LOG("ERROR in 'IExportedClass::check_consistency':"
					" base class "<<i<<" of class '"<<this->name()<<
					"' has not been named.\n");
				return false;
			}
			else{
			UG_LOG("ERROR in 'IExportedClass::check_consistency':"
					" Class '"<<this->name()<<"' has not been named.\n");
				return false;
			}
		}
	}

//	everything ok
	return true;
}


} // end namespace ug
} // end namespace bridge
