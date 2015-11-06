/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
 * Authors: Sebastian Reiter, Andreas Vogel
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

#ifndef __H__UG_BRIDGE__CLASS_NAME_PROVIDER_IMPL__
#define __H__UG_BRIDGE__CLASS_NAME_PROVIDER_IMPL__

#include <algorithm>
#include <string>

#include "class_name_provider.h"
#include "common/util/typename.h"

namespace ug{
namespace bridge{

////////////////////////////////////////////////////////////////////////////////
//  ClassNameProvider
////////////////////////////////////////////////////////////////////////////////

template <typename TClass>
void ClassNameProvider<TClass>::
set_name(const std::string& nameIn, const std::string& group, bool newName)
{
//	if class already named throw error
	if(newName == true && m_bForwardDeclared==false && !m_ClassNameNode.empty())
	{
		if(nameIn != name())
			UG_THROW_REGISTRY_ERROR(nameIn,
			"Trying to register class with name '"<<nameIn
			<< "', that has already been named. This is not allowed.");
	}

//	check name
	if(nameIn.empty())
		UG_THROW_REGISTRY_MSG("Registered class name has length zero");

	if(nameIn.c_str()[0] == '[')
		UG_THROW_REGISTRY_ERROR(nameIn, "Registered class name must not begin with '['");

//	check for non-allowed character
	size_t found = nameIn.find_first_not_of("abcdefghijklmnopqrstuvwxyz"
											 "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
											 "_"
											 "0123456789");
	if (found!=std::string::npos)
	{
		UG_ERR_LOG("Non-allowed character '"<<nameIn[found]<<"' "<<
		       "contained at position "<<int(found)<<" in registered Class Name "
		       "'"<<nameIn<<"'.\nClass names must match the regular expression: "
		       "[a-zA-Z_][a-zA-Z_0-9]*, \ni.e. only alphabetic characters, numbers "
		       " and '_' are allowed; no numbers at the beginning.\n");
		UG_THROW_REGISTRY_ERROR(nameIn, "Class Name must only contain [a-zA-Z_][a-zA-Z_0-9]*.");
	}

//	check that no number at the beginning
	 found = nameIn.find_first_of("0123456789");
	if (found!=std::string::npos && found == 0)
	{
		UG_ERR_LOG("Class Name "<<nameIn<<" starts with a number.\nThis is "
				" not allowed. Please change naming.\n");
		UG_THROW_REGISTRY_ERROR(nameIn,"Class Name must not start with number.");
	}


//	copy name into static string
	m_ClassNameNode.set_name(nameIn);

//	set forward declared to false, since now name set
	m_bForwardDeclared = false;

//	remember groups
	m_group = std::string(group);
};

template <typename TClass>
template <typename TParent1>
void ClassNameProvider<TClass>::
set_name(const std::string& name, const std::string& group, bool newName)
{
//	set own name
	set_name(name, group, newName);

//	add parent nodes
	m_ClassNameNode.add_base_class(ClassNameProvider<TParent1>::class_name_node());
}

template <typename TClass>
template <typename TParent1, typename TParent2>
void ClassNameProvider<TClass>::
set_name(const std::string& name, const std::string& group,bool newName)
{
//	set own name
	set_name(name, group, newName);

//	add parent nodes
	m_ClassNameNode.add_base_class(ClassNameProvider<TParent1>::class_name_node());
	m_ClassNameNode.add_base_class(ClassNameProvider<TParent2>::class_name_node());
}

template <typename TClass>
bool ClassNameProvider<TClass>::is_a(const std::string& parent, bool strict)
{
//	check if class is forward declared
	if(m_bForwardDeclared)
		UG_THROW_REGISTRY_ERROR(parent,
		"Class '"<<parent<<"' must not be foreward declared to use is_a");

//  strict comparison: must match this class name, parents are not allowed
	if(strict)
	{
	//  compare strings
		if(parent == name()) return true;

	//  names does not match
		return false;
	}

//	return if parent name is contained in tree of base classes
	return ClassNameTreeContains(m_ClassNameNode, parent);
}

template <typename TClass>
const std::string& ClassNameProvider<TClass>::name()
{
//	if name has not been set, set temporary forward declaration
	if(m_ClassNameNode.empty()) set_foreward_declared();

//	return the name of this class as stored in ClassNameNode
	return m_ClassNameNode.name();
}

template <typename TClass>
const std::vector<const char*>& ClassNameProvider<TClass>::names()
{
//	if name has not been set, set temporary forward declaration
	if(m_ClassNameNode.empty()) set_foreward_declared();

//	create names list, including this class and all base classes
//	\todo: remove this, avoid names-vector
	ExtractClassNameVec(m_names, m_ClassNameNode, true);

//	return the vector of base class names
	return m_names;
}

template <typename TClass>
void ClassNameProvider<TClass>::set_foreward_declared()
{
//	build default name using typeinfo
	std::string name("[[");
	name.append( TypeName<TClass>() );
	name.append(" (undeclared) ]]");

//	set this as current name, but remember pre-declaration
	m_ClassNameNode.set_name(name);
	m_bForwardDeclared = true;
}

template <typename TClass>
std::vector<const char*> ClassNameProvider<TClass>::m_names = std::vector<const char*>();

template <typename TClass>
std::string ClassNameProvider<TClass>::m_group = std::string("");

template <typename TClass>
bool ClassNameProvider<TClass>::m_bForwardDeclared = false;

template <typename TClass>
ClassNameNode ClassNameProvider<TClass>::m_ClassNameNode = ClassNameNode();

////////////////////////////////////////////////////////////////////////////////
//  ClassCastProvider
////////////////////////////////////////////////////////////////////////////////

template <typename TBase, typename TDerived>
void* StaticVoidCast(void* DerivVoidPtr)
{
//	cast to derived class; this assumes, that the void pointer points to the
//	beginning of the data field of the Derived object
	TDerived* pDeriv = reinterpret_cast<TDerived*>(DerivVoidPtr);

//	static case to the Derid class
	TBase* pBase = static_cast<TBase*>(pDeriv);

//	return cast to void
	return reinterpret_cast<void*>(pBase);
}

template <typename TBase, typename TDerived>
void ClassCastProvider::add_cast_func()
{
	const ClassNameNode* pBaseNode =
						&ClassNameProvider<TBase>::class_name_node();
	const ClassNameNode* pDerivNode =
						&ClassNameProvider<TDerived>::class_name_node();

	std::pair<const ClassNameNode*, const ClassNameNode*>
										namePair(pBaseNode, pDerivNode);

	m_mmCast[namePair] = &StaticVoidCast<TBase, TDerived>;
}

template <typename T>
T* ClassCastProvider::
cast_to(void* ptr, const ClassNameNode*& node)
{
//	get base name
	const std::string& baseName = ClassNameProvider<T>::name();

//	cast the plain pointer
	ptr = ClassCastProvider::cast_to_base_class(ptr, node, baseName);

//	return it
	return reinterpret_cast<T*>(ptr);
}

template <typename T>
const T* ClassCastProvider::
cast_to(const void* ptr, const ClassNameNode*& node)
{
//	get base name
	const std::string& baseName = ClassNameProvider<T>::name();

//	cast the plain pointer
	ptr = ClassCastProvider::cast_to_base_class(ptr, node, baseName);

//	return it
	return reinterpret_cast<const T*>(ptr);
}

template <typename T>
SmartPtr<T> ClassCastProvider::
cast_to(SmartPtr<void> spDerivVoid, const ClassNameNode*& node)
{
//	get base name
	const std::string& baseName = ClassNameProvider<T>::name();

//	extract plain pointer
	void* rawPtr = spDerivVoid.get();

//	cast the plain pointer
	rawPtr = ClassCastProvider::cast_to_base_class(rawPtr, node, baseName);

//	sets as pointer to the smart ptr
	spDerivVoid.set_impl<T, FreeDelete>(rawPtr);

//	return it
	return spDerivVoid.cast_reinterpret<T, FreeDelete>();
}

template <typename T>
ConstSmartPtr<T> ClassCastProvider::
cast_to(ConstSmartPtr<void> spDerivVoid, const ClassNameNode*& node)
{
//	get base name
	const std::string& baseName = ClassNameProvider<T>::name();

//	extract plain pointer
	const void* rawPtr = spDerivVoid.get();

//	cast the plain pointer
	rawPtr = ClassCastProvider::cast_to_base_class(rawPtr, node, baseName);

//	sets as pointer to the smart ptr
	spDerivVoid.set_impl<T, FreeDelete>(rawPtr);

//	return it
	return spDerivVoid.cast_reinterpret<T, FreeDelete>();
}
} // end namespace ug
} // end namespace bridge

#endif /* __H__UG_BRIDGE__CLASS_NAME_PROVIDER_IMPL__ */
