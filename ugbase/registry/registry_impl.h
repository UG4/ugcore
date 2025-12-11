/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
 * Authors: Andreas Vogel, Sebastian Reiter
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

#ifndef __H__UG_BRIDGE__REGISTRY_IMPL__
#define __H__UG_BRIDGE__REGISTRY_IMPL__

#include "registry_util.h"
#include "common/util/typename.h"

namespace ug{
namespace bridge
{

//////////////////////
// global functions
//////////////////////
// template<typename TClass>
// struct AddTypeName
// {
// 	static void add(std::string &s)
// 	{
// 	#ifdef UG_POSIX
// 		if(s.length() != 0)
// 			s += std::string(", ");
// 		s += std::string("C++ Name: ") + TypeName<TClass>();
// 	#endif
// 	}
// };

template<class TFunc>
Registry& Registry::
add_function(std::string funcName, TFunc func, std::string group,
					 std::string retValInfos, std::string paramInfos,
					 std::string tooltip, std::string help)
{
	add_and_get_function(funcName, func, group, retValInfos, paramInfos,
						 tooltip, help);
	return *this;
}

template<class TFunc>
ExportedFunction* Registry::
add_and_get_function(std::string funcName, TFunc func, std::string group,
					 std::string retValInfos, std::string paramInfos,
					 std::string tooltip, std::string help)
{
//	At this point the method name contains parameters (name|param1=...).
//todo: they should be removed and specified with an extra parameter.

	std::string strippedMethodName = funcName;
	std::string methodOptions;
	std::string::size_type pos = strippedMethodName.find("|");
	if(pos != std::string::npos){
		methodOptions = strippedMethodName.substr(pos + 1, strippedMethodName.length() - pos);
		strippedMethodName = strippedMethodName.substr(0, pos);
	}

//	trim whitespaces
	strippedMethodName = TrimString(strippedMethodName);
	methodOptions = TrimString(methodOptions);

// 	check that name is not empty
	if(strippedMethodName.empty())
	{
		UG_THROW_REGISTRY_ERROR(strippedMethodName,
		"Trying to register empty function name.");
	}
	
	// check that name does not contain illegal characters
	if (!IsValidRegistryIdentifier(strippedMethodName)) {
		UG_THROW_REGISTRY_ERROR(strippedMethodName,
		"Trying to register function '" << strippedMethodName << "' that"
		<< " contains illegal characters. " << GetRegistryIdentifierMessage());
	}

//	if the function is already in use, we have to add an overload
	ExportedFunctionGroup* funcGrp = get_exported_function_group(strippedMethodName);
	if(!funcGrp)
	{
	//	we have to create a new function group
		funcGrp = new ExportedFunctionGroup(strippedMethodName);
		m_vFunction.push_back(funcGrp);
	}

//  add an overload to the function group
	ExportedFunction* ef = funcGrp->add_overload(func, &FunctionProxy<TFunc>::apply,
	                                     		 methodOptions, group,
	                                     		 retValInfos, paramInfos,
	                                     		 tooltip, help);

	if(!ef){
		UG_THROW_REGISTRY_ERROR(strippedMethodName,
		"Trying to register function name '"<<funcName
		<< "', that is already used by another function in this registry.");
	}

	return ef;
}


///////////////////
// classes
///////////////////

template <typename TClass, typename TBaseClass>
void Registry::
check_base_class(const std::string& className)
{
//	check that className is not already used
	if(classname_registered(className))
	{
		UG_THROW_REGISTRY_ERROR(className,
		"Trying to register class name '"<<className
		<< "', that is already used by another class in this registry.");
	}
// 	check that name is not empty
	if(className.empty())
	{
		UG_THROW_REGISTRY_ERROR(className,
		"Trying to register empty class name.");
	}

// 	check that base class is not same type as class
	if(typeid(TClass) == typeid(TBaseClass))
	{
		UG_THROW_REGISTRY_ERROR(className,
		"Trying to register class '"<< className
				<< "' that derives from itself.");
	}

// 	check that class derives from base class
	if(boost::is_base_of<TBaseClass, TClass>::value == false)
	{
		UG_THROW_REGISTRY_ERROR(className,
		"Trying to register class "<<className
		<< "with base class that is no base class.");
	}
}

template <typename TClass>
ExportedClass<TClass>& Registry::
add_class_(std::string className, std::string group, std::string tooltip)
{
//	check that className is not already used
	if(classname_registered(className))
	{
		UG_THROW_REGISTRY_ERROR(className,
		"Trying to register class name '"<<className
		<< "', that is already used by another class in this registry.");
	}
//	check that className is not already used as a group name
	if(groupname_registered(className))
	{
		UG_THROW_REGISTRY_ERROR(className,
		"Trying to register class name '"<<className
		<< "', that is already used by another group in this registry.");
	}
// 	check that name is not empty
	if(className.empty())
	{
		UG_THROW_REGISTRY_ERROR(className,
		"Trying to register empty class name.");
	}
	
	// check that name does not contain illegal characters
	if (!IsValidRegistryIdentifier(className)) {
		UG_THROW_REGISTRY_ERROR(className,
		"Trying to register class '" << className << "' that"
		<< " contains illegal characters. "<< GetRegistryIdentifierMessage());
	}

//	new class pointer
	ExportedClass<TClass>* newClass = NULL;

	// AddTypeName<TClass>::add(tooltip);
	newClass = new ExportedClass<TClass>(className, group, tooltip);

//	add new class to list of classes
	m_vClass.push_back(newClass);
#if defined(FEATURE_REGISTRY_CLASS_NAME_MAP) && FEATURE_REGISTRY_CLASS_NAME_MAP == 1
	m_classMap[className] = newClass;
#endif
	return *newClass;
}

template <typename TClass, typename TBaseClass>
ExportedClass<TClass>& Registry::
add_class_(std::string className, std::string group, std::string tooltip)
{
//	check that className is not already used
	if(classname_registered(className))
	{
		UG_THROW_REGISTRY_ERROR(className,
		"Trying to register class name '"<<className
		<< "', that is already used by another class in this registry.");
	}
//	check that className is not already used as a group name
	if(groupname_registered(className))
	{
		UG_THROW_REGISTRY_ERROR(className,
		"Trying to register class name '"<<className
		<< "', that is already used by another group in this registry.");
	}
// 	check that name is not empty
	if(className.empty())
	{
		UG_THROW_REGISTRY_ERROR(className,
		"Trying to register empty class name.");
	}
	
	// check that name does not contain illegal characters
	if (!IsValidRegistryIdentifier(className)) {
		UG_THROW_REGISTRY_ERROR(className,
		"Trying to register class '" << className << "' that"
		<< " contains illegal characters. "<< GetRegistryIdentifierMessage());
	}

//	check
	check_base_class<TClass, TBaseClass>(className);

//	new class pointer
	ExportedClass<TClass>* newClass = NULL;

	// AddTypeName<TClass>::add(tooltip);

//	try creation of new class
	newClass = new ExportedClass<TClass>(className, group, tooltip);

// 	set base class names
	ClassNameProvider<TClass>::template set_name<TBaseClass>(className, group);

//	add cast function
	ClassCastProvider::add_cast_func<TBaseClass, TClass>();

//	add new class to list of classes
	m_vClass.push_back(newClass);
#if defined(FEATURE_REGISTRY_CLASS_NAME_MAP) && FEATURE_REGISTRY_CLASS_NAME_MAP == 1
	m_classMap[className] = newClass;
#endif
	return *newClass;
}

template <typename TClass, typename TBaseClass1, typename TBaseClass2>
ExportedClass<TClass>& Registry::
add_class_(std::string className, std::string group, std::string tooltip)
{
//	check that className is not already used
	if(classname_registered(className))
	{
		UG_THROW_REGISTRY_ERROR(className,
		"Trying to register class name '"<<className
		<< "', that is already used by another class in this registry.");
	}
//	check that className is not already used as a group name
	if(groupname_registered(className))
	{
		UG_THROW_REGISTRY_ERROR(className,
		"Trying to register class name '"<<className
		<< "', that is already used by another group in this registry.");
	}
// 	check that name is not empty
	if(className.empty())
	{
		UG_THROW_REGISTRY_ERROR(className,
		"Trying to register empty class name.");
	}
	
	// check that name does not contain illegal characters
	if (!IsValidRegistryIdentifier(className)) {
		UG_THROW_REGISTRY_ERROR(className,
		"Trying to register class '" << className << "' that"
		<< " contains illegal characters. " << GetRegistryIdentifierMessage());
	}

//	check
	check_base_class<TClass, TBaseClass1>(className);
	check_base_class<TClass, TBaseClass2>(className);

//	new class pointer
	ExportedClass<TClass>* newClass = NULL;

	// AddTypeName<TClass>::add(tooltip);

//	try creation of new class
	newClass = new ExportedClass<TClass>(className, group, tooltip);

// 	set base class names
	ClassNameProvider<TClass>::template set_name<TBaseClass1, TBaseClass2>(className, group);

//	add cast function
	ClassCastProvider::add_cast_func<TBaseClass1, TClass>();
	ClassCastProvider::add_cast_func<TBaseClass2, TClass>();

//	add new class to list of classes
	m_vClass.push_back(newClass);
	return *newClass;
}

template <typename TClass>
ExportedClass<TClass>& Registry::
get_class_()
{
// 	get class names
	const std::string& name = ClassNameProvider<TClass>::name();

//	look for class in this registry
	for(size_t i = 0; i < m_vClass.size(); ++i)
		if(name == m_vClass[i]->name())
			return *dynamic_cast<ExportedClass<TClass>* >(m_vClass[i]);

//	not found
	UG_THROW_REGISTRY_ERROR(name,
	"Trying to get class with name '" << name
	<< "', that has not yet been registered to this Registry.");
}

}//	end of namespace
}//	end of namespace

#endif
