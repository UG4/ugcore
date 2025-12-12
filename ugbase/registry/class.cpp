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

#include "class.h"

namespace ug {
namespace bridge {

////////////////////////////////////////////////////////////////////////////////
// ExportedConstructor
////////////////////////////////////////////////////////////////////////////////

ExportedConstructor::
ExportedConstructor(ProxyFunc pf,
                    const std::string& className, const std::string& options,
                    const std::string& paramInfos,
                    const std::string& tooltip, const std::string& help)
: m_proxy_func(pf), m_className(className),
  m_options(options), m_paramInfos(paramInfos), m_vvParamInfo(0),
  m_tooltip(tooltip), m_help(help)
#ifdef PROFILE_BRIDGE
  ,m_dpi((m_className + "(...)").c_str(), true, "registry", false)
#endif
{

//	Tokenize string for parameters into infos per one parameter (separated by '#')
	std::vector<std::string> vParamInfoTmp;
	tokenize(m_paramInfos, vParamInfoTmp, '#');
	m_vvParamInfo.resize(vParamInfoTmp.size());

//	Tokenize each info-string of one parameter into single infos (separated by '|')
	for(size_t i = 0; i < vParamInfoTmp.size(); ++i)
		tokenize(vParamInfoTmp[i], m_vvParamInfo[i], '|');
};

bool ExportedConstructor::check_consistency(const std::string &classname) const
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
				UG_ERR_LOG("#### Registry ERROR: Unregistered Class used in ");
				UG_ERR_LOG("Constructor of class " << classname.c_str());
				UG_ERR_LOG("': Parameter " << j+1);
			}
			else
			{	UG_ERR_LOG(", " << j+1);	}
		}
	}

//	check if undeclared parameter has been found
	if(bUndeclaredParameterFound) {UG_ERR_LOG("\n"); return false;}

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
	if(vClassNames==nullptr)
	{
		UG_ERR_LOG("#### Registry ERROR:"
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
		if(baseName == nullptr || *baseName == '\0' || baseName[0] == '[')
		{
			if(i>0){
			UG_ERR_LOG("#### Registry ERROR:"
					" base class "<<i<<" of class '"<<this->name()<<
					"' has not been named.\n");
				return false;
			}
			else{
			UG_ERR_LOG("#### Registry ERROR:"
					" Class '"<<this->name()<<"' has not been named.\n");
				return false;
			}
		}
	}

//	everything ok
	return true;
}

////////////////////////////////////////////////////////////////////////////////
// Implementation of ExportedClassBaseImpl
////////////////////////////////////////////////////////////////////////////////


ExportedClassBaseImpl::
ExportedClassBaseImpl(const std::string& tooltip)
	: m_destructor(nullptr), m_tooltip(tooltip), m_constructAsSmartPtr(false)
{
}

ExportedClassBaseImpl::
~ExportedClassBaseImpl()
{
//	delete constructors
	for(size_t i = 0; i < m_vConstructor.size(); ++i)
		delete (m_vConstructor[i].m_constructor);

//  delete methods
	for(size_t i = 0; i < m_vMethod.size(); ++i)
		delete m_vMethod[i];

	for(size_t i = 0; i < m_vConstMethod.size(); ++i)
		delete m_vConstMethod[i];
}

const std::string& ExportedClassBaseImpl::
tooltip() const
{
	return m_tooltip;
}

size_t ExportedClassBaseImpl::
num_methods() const
{
	return m_vMethod.size();
}

size_t ExportedClassBaseImpl::
num_const_methods() const
{
	return m_vConstMethod.size();
}

const ExportedMethod& ExportedClassBaseImpl::
get_method(size_t i) const
{
	return *m_vMethod[i]->get_overload(0);
}

const ExportedMethod& ExportedClassBaseImpl::
get_const_method(size_t i) const
{
	return *m_vConstMethod[i]->get_overload(0);
}

size_t ExportedClassBaseImpl::
num_overloads(size_t funcInd) const
{
	return m_vMethod[funcInd]->num_overloads();
}

size_t ExportedClassBaseImpl::
num_const_overloads(size_t funcInd) const
{
	return m_vConstMethod[funcInd]->num_overloads();
}

const ExportedMethod& ExportedClassBaseImpl::
get_overload(size_t funcInd, size_t oInd) const
{
	return *m_vMethod[funcInd]->get_overload(oInd);
}

const ExportedMethod& ExportedClassBaseImpl::
get_const_overload(size_t funcInd, size_t oInd) const
{
	return *m_vConstMethod[funcInd]->get_overload(oInd);
}

const ExportedMethodGroup& ExportedClassBaseImpl::
get_method_group(size_t ind) const
{
	return *m_vMethod[ind];
}

const ExportedMethodGroup& ExportedClassBaseImpl::
get_const_method_group(size_t ind) const
{
	return *m_vConstMethod[ind];
}

size_t ExportedClassBaseImpl::
num_constructors() const
{
	return m_vConstructor.size();
}

const ExportedConstructor& ExportedClassBaseImpl::
get_constructor(size_t i) const
{
	return *(m_vConstructor[i].m_constructor);
}

const boost::optional<ExportedConstructor&> ExportedClassBaseImpl::
get_json_constructor() const
{
	if(!is_json_constructible()) 
		return boost::none;

	for(size_t i = 0; i < m_vConstructor.size(); ++i)
		if(m_vConstructor[i].m_typeID == GetUniqueTypeID<void (*)(const char*)>())
			return *(m_vConstructor[i].m_constructor);

	return boost::none;
}

bool ExportedClassBaseImpl::
construct_as_smart_pointer() const
{
	return m_constructAsSmartPtr;
}

void ExportedClassBaseImpl::
set_construct_as_smart_pointer(bool enable)
{
	m_constructAsSmartPtr = enable;
}

bool ExportedClassBaseImpl::
is_instantiable() const
{
	return !m_vConstructor.empty();
}

void ExportedClassBaseImpl::
destroy(void* obj) const
{
	if(m_destructor != nullptr)
		(*m_destructor)(obj);
}

bool ExportedClassBaseImpl::
constructor_type_id_registered(size_t typeID) const
{
	for(size_t i = 0; i < m_vConstructor.size(); ++i)
		if(typeID == m_vConstructor[i].m_typeID)
			return true;

	return false;
}

bool ExportedClassBaseImpl::
constmethodname_registered(const std::string& name) const
{
	for(size_t i = 0; i < m_vConstMethod.size(); ++i)
		if(name == m_vConstMethod[i]->name())
			return true;

	return false;
}

bool ExportedClassBaseImpl::
methodname_registered(const std::string& name) const
{
	for(size_t i = 0; i < m_vMethod.size(); ++i)
		if(name == m_vMethod[i]->name())
			return true;

	return false;
}

ExportedMethodGroup* ExportedClassBaseImpl::
get_exported_method_group(const std::string& name)
{
	for(size_t i = 0; i < m_vMethod.size(); ++i)
		if(name == m_vMethod[i]->name())
			return m_vMethod[i];

	return nullptr;
}

const ExportedMethodGroup* ExportedClassBaseImpl::
get_exported_method_group(const std::string& name) const
{
	for(size_t i = 0; i < m_vMethod.size(); ++i)
		if(name == m_vMethod[i]->name())
			return m_vMethod[i];

	return nullptr;
}

ExportedMethodGroup* ExportedClassBaseImpl::
get_const_exported_method_group(const std::string& name)
{
	for(size_t i = 0; i < m_vConstMethod.size(); ++i)
		if(name == m_vConstMethod[i]->name())
			return m_vConstMethod[i];

	return nullptr;
}

const ExportedMethodGroup* ExportedClassBaseImpl::
get_const_exported_method_group(const std::string& name) const
{
	for(size_t i = 0; i < m_vConstMethod.size(); ++i)
		if(name == m_vConstMethod[i]->name())
			return m_vConstMethod[i];

	return nullptr;
}

} // end namespace ug
} // end namespace bridge
