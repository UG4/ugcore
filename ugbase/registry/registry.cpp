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

#include <string>

#include "registry.h"
#include "registry_util.h"
#ifdef UG_FOR_LUA
#include "bindings/lua/lua_function_handle.h"
#endif

using namespace ug::bridge;

namespace ug{
namespace bridge
{

Registry::Registry()
	:m_bForceConstructionWithSmartPtr(false)
{
//	register native types as provided in ParameterStack
//	we use the c_ prefix to avoid clashes with java native types in java bindings.
	add_class_<bool>("c_bool");
	add_class_<int>("c_int");
	add_class_<size_t>("c_size_t");
	add_class_<float>("c_float");
	add_class_<double>("c_double");
	add_class_<char>("c_char");
	add_class_<const char*>("c_const_char_ptr");
	add_class_<std::string>("c_string");
	add_class_<void>("c_void");
#ifdef UG_FOR_LUA
	add_class_<LuaFunctionHandle>("c_LuaFunctionHandle");
#endif
}

Registry::Registry(const Registry& reg)
	: m_bForceConstructionWithSmartPtr(false)
{
}

Registry::~Registry()
{
//  delete registered functions
	for(size_t i = 0; i < m_vFunction.size(); ++i)
	{
		if(m_vFunction[i] != nullptr)
			delete m_vFunction[i];
	}
//  delete registered classes
	for(size_t i = 0; i < m_vClass.size(); ++i)
	{
		if(m_vClass[i] != nullptr)
			delete m_vClass[i];
	}
//	delete registered class groups
	for(size_t i = 0; i < m_vClassGroups.size(); ++i){
		delete m_vClassGroups[i];
	}
}


////////////////////////
//	callbacks
////////////////////////

void Registry::add_callback(FuncRegistryChanged callback)
{
	m_callbacksRegChanged.push_back(callback);
}

bool Registry::registry_changed()
{
	if(!m_callbacksRegChanged.empty())
		UG_THROW("Sorry, currently no listeners can be registered at the "
				" the Registry due to restrictions of the VRL binding.");

//	iterate through all callbacks and call them
	for(size_t i = 0; i < m_callbacksRegChanged.size(); ++i){
		m_callbacksRegChanged[i](this);
	}

//	ok
	return true;
}

//////////////////////
// global functions
//////////////////////

size_t Registry::num_functions() const
{
	return m_vFunction.size();
}

ExportedFunction& Registry::get_function(size_t ind)
{
	return *m_vFunction.at(ind)->get_overload(0);
}

const ExportedFunction& Registry::get_function(size_t ind) const
{
	return *m_vFunction.at(ind)->get_overload(0);
}

size_t Registry::num_overloads(size_t ind) const
{
	return m_vFunction.at(ind)->num_overloads();
}

ExportedFunction& Registry::get_overload(size_t funcInd, size_t oInd) const {
	return *m_vFunction.at(funcInd)->get_overload(oInd);
}

ExportedFunctionGroup& Registry::get_function_group(size_t ind) const {
	return *m_vFunction.at(ind);
}


///////////////////
// classes
///////////////////

size_t Registry::num_classes() const
{
	return m_vClass.size();
}

/// returns an exported function
const IExportedClass& Registry::get_class(size_t ind) const
{
	return *m_vClass.at(ind);
}

IExportedClass* Registry::get_class(const std::string& name)
{
//todo:	use a map to access classes by name.
	for(size_t i = 0; i < m_vClass.size(); ++i)
		if(name == m_vClass[i]->name())
			return m_vClass[i];

	return nullptr;
}

const IExportedClass* Registry::get_class(const std::string& name) const
{
//todo:	use a map to access classes by name.
	for(size_t i = 0; i < m_vClass.size(); ++i)
		if(name == m_vClass[i]->name())
			return m_vClass[i];

	return nullptr;
}


void Registry::set_force_construct_via_smart_pointer(bool bForceConstructionWithSmartPtr)
{
	m_bForceConstructionWithSmartPtr = bForceConstructionWithSmartPtr;
}

bool Registry::check_consistency()
{

	// list to check for duplicates in global functions.
	// comparison is not case-sensitive.
	std::vector<std::string> globalFunctionNames;
    
	//	check global functions
	size_t globFctUndef = 0;
	for(size_t i=0; i<num_functions(); i++)
	{
		//	get function
		ExportedFunctionGroup& funcGrp = get_function_group(i);
		
		// add lower case name entry
		globalFunctionNames.push_back(ToLower(funcGrp.name()));

		//	check all overloads
		for(size_t j = 0; j < funcGrp.num_overloads(); ++j){
			if(!funcGrp.get_overload(j)->check_consistency())
				globFctUndef++;
		}
	}
	
	// check for duplicate function names. comparison is not case-sensitive.
	std::vector<std::string> globFuncDuplicates = 
			FindDuplicates(globalFunctionNames);
	bool globFuncDuplicatesExist = !globFuncDuplicates.empty();
	std::string duplicateFuncMsg = "#### Registry ERROR: duplicate function names:\n";
	if (globFuncDuplicatesExist) {
		for(size_t i = 0; i < globFuncDuplicates.size();i++) {
			duplicateFuncMsg+=  + "\t" + globFuncDuplicates[i] + "\n";
		}
		duplicateFuncMsg += "#### NOTE: it is not allowed to register multiple"
				" functions with equal names. Comparison is NOT case sensitive,"
				" e.g., 'Func' and 'func' are equal.\n\n";
	}

	set_force_construct_via_smart_pointer(false);

// 	check classes and their methods
	size_t baseClassUndef = 0;
	size_t methodUndef = 0;
	size_t constructorUndef = 0;
	size_t notConstructedAsSmartPtr = 0;
	// list to check for duplicate classes.
	// comparison is not case-sensitive.
	std::vector<std::string> classNames;
	for(size_t i=0; i<num_classes(); i++)
	{
	//	get class
		const IExportedClass &c = get_class(i);
		
	//	check if class is constructed via smart pointer
		if(m_bForceConstructionWithSmartPtr)
		{
			if(c.is_instantiable() && !c.construct_as_smart_pointer())
			{
				UG_ERR_LOG("#### Registry ERROR: Class " << c.name()<<" is "
				       "instantiable, but not constructed as Smart Pointer.\n");
				notConstructedAsSmartPtr++;
			}
		}

	// 	add lower case name entry
		classNames.push_back(ToLower(c.name()));

	//	check class (e.g. that base classes have been named)
		if(!c.check_consistency()){
			UG_ERR_LOG("#### Registry ERROR: Base Class Error for class " << c.name()<<"\n");
			baseClassUndef++;
		}

	//	check methods
		for(size_t j=0; j<c.num_methods(); j++)
			if(!c.get_method(j).check_consistency(c.name()))
				methodUndef++;

		for(size_t j=0; j<c.num_const_methods(); j++)
			if(!c.get_const_method(j).check_consistency(c.name()))
				methodUndef++;

	//	check constructors
		for(size_t j=0; j<c.num_constructors(); j++)
			if(!c.get_constructor(j).check_consistency(c.name()))
				constructorUndef++;
	}
	
	// check for duplicate class names. comparison is not case-sensitive.
	std::vector<std::string> classDuplicates = 
			FindDuplicates(classNames);
	bool classDuplicatesExist = !classDuplicates.empty();
	//todo: use stringstream
	std::string duplicateClassMsg = 
				"#### Registry ERROR: duplicate class names:\n";
	if (classDuplicatesExist) {
		for(size_t i = 0; i < classDuplicates.size();i++) {
			duplicateClassMsg+= "\t" + classDuplicates[i] + "\n";
		}
		duplicateClassMsg += "#### NOTE: it is not allowed to register multiple"
				"  classes with equal names. Comparison is NOT case sensitive,"
				" e.g., 'Class' and 'class' are equal.\n\n";
	}

//	log error messages
	if(globFctUndef > 0)
		UG_ERR_LOG("#### Registry ERROR: "<<globFctUndef<<
		       " global Functions are using undeclared Classes.\n");
	if(methodUndef > 0)
		UG_ERR_LOG("#### Registry ERROR: "<<methodUndef<<
		       " Methods are using undeclared Classes.\n");
	if(baseClassUndef > 0)
		UG_ERR_LOG("#### Registry ERROR: "<<baseClassUndef<<
		       " Base-Classes are using undeclared Classes.\n");
	if(constructorUndef > 0)
		UG_ERR_LOG("#### Registry ERROR: "<<constructorUndef<<
		       " Constructors are using undeclared Classes.\n");
	if(notConstructedAsSmartPtr > 0)
		UG_ERR_LOG("#### Registry ERROR: "<<notConstructedAsSmartPtr<<
		       " Classes are not constructed as SmartPtr.\n");
	if(globFuncDuplicatesExist) {
		UG_ERR_LOG(duplicateFuncMsg);
	}
	if(classDuplicatesExist) {
		UG_ERR_LOG(duplicateClassMsg);
	}

//	return false if undefined classes or functions and/or duplicates have been found
	if(globFctUndef > 0) return false;
	if(methodUndef > 0) return false;
	if(baseClassUndef > 0) return false;
	if(constructorUndef > 0) return false;
	if(notConstructedAsSmartPtr > 0) return false;
	if(globFuncDuplicatesExist) return false;
	if(classDuplicatesExist) return false;

//	everything fine
	return true;
}


size_t Registry::num_class_groups() const
{
	return m_vClassGroups.size();
}

const ClassGroupDesc* Registry::get_class_group(size_t i) const
{
	return m_vClassGroups[i];
}

ClassGroupDesc* Registry::get_class_group(size_t i)
{
	return m_vClassGroups[i];
}

ClassGroupDesc* Registry::get_class_group(const std::string& name)
{
//todo:	use a map to quickly access classGroups by name
	for(size_t i = 0; i < m_vClassGroups.size(); ++i)
		if(name == m_vClassGroups[i]->name())
			return m_vClassGroups[i];

//	since we reached this point, no class-group with the given name exists.
	
	// check that name does not contain illegal characters
	if (!IsValidRegistryIdentifier(name)) {
		UG_THROW_REGISTRY_ERROR(name,
		"Trying add group '" << name << "' that"
		<< " contains illegal characters. "<< GetRegistryIdentifierMessage());
	}

	auto* classGroup = new ClassGroupDesc();
	classGroup->set_name(name);
	m_vClassGroups.push_back(classGroup);

	return classGroup;
}

const ClassGroupDesc* Registry::get_class_group(const std::string& name) const
{
//todo:	use a map to quickly access classGroups by name
	for(size_t i = 0; i < m_vClassGroups.size(); ++i)
		if(name == m_vClassGroups[i]->name())
			return m_vClassGroups[i];

//	since we reached this point, no class-group with the given name exists.
	return nullptr;
}

void Registry::add_class_to_group(const std::string& className,
								  const std::string& groupName,
								  const std::string& classTag)
{
//	make sure that no class with groupName exists.
	if(classname_registered(groupName)){
		UG_THROW_REGISTRY_ERROR(groupName,
		"A class with the given group name '"<<groupName<<"' already exists.");
	}

	ClassGroupDesc* groupDesc = get_class_group(groupName);
//todo:	make sure that groupDesc does not already contain className.
	IExportedClass* expClass = get_class(className);
	if(!expClass){
		UG_THROW_REGISTRY_ERROR(groupName,
		"The given class has to be registered before "
						"adding it to a group: " << className);
	}

	if(expClass)
		groupDesc->add_class(expClass, classTag);
}


bool Registry::classname_registered(const std::string& name)
{
	return get_class(name) != nullptr;
}

bool Registry::groupname_registered(const std::string& name) const {
//todo:	use a map to quickly access classGroups by name
	for(size_t i = 0; i < m_vClassGroups.size(); ++i){
		if(name == m_vClassGroups[i]->name())
			return true;
	}
	return false;
}

// returns true if functionname is already used by a function in this registry
bool Registry::functionname_registered(const std::string& name) const {
	for(size_t i = 0; i < m_vFunction.size(); ++i)
		if(name == m_vFunction[i]->name())
			return true;

	return false;
}

ExportedFunctionGroup* Registry::get_exported_function_group(const std::string& name) const {
	for(size_t i = 0; i < m_vFunction.size(); ++i)
		if(name == m_vFunction[i]->name())
			return m_vFunction[i];

	return nullptr;
}

}// end of namespace
}// end of namespace
