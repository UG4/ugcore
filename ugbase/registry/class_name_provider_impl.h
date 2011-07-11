/*
 * class_name_provider_impl.h
 *
 *  Created on: 08.07.2011
 *      Author: sreiter, andreasvogel
 */

#ifndef __H__UG_BRIDGE__CLASS_NAME_PROVIDER_IMPL__
#define __H__UG_BRIDGE__CLASS_NAME_PROVIDER_IMPL__

#include <algorithm>
#include <string>

#include "class_name_provider.h"

namespace ug{
namespace bridge{

////////////////////////////////////////////////////////////////////////////////
//  ClassNameProvider
////////////////////////////////////////////////////////////////////////////////

template <typename TClass>
void ClassNameProvider<TClass>::
set_name(const char* nameIn, const char* group, bool newName)
{
//	if class already named throw error
	if(newName == true && m_bForwardDeclared==false && !m_ClassNameNode.empty())
	{
		if(strcmp(nameIn, name()) != 0)
			throw(REGISTRY_ERROR_ClassAlreadyNamed(nameIn));
	}

//	check name
	if(nameIn == NULL)
		throw(REGISTRY_ERROR_Message("Registered class name is NULL"));

	if(strlen(nameIn) == 0)
		throw(REGISTRY_ERROR_Message("Registered class name has length zero"));

	if(nameIn[0] == '[')
		throw(REGISTRY_ERROR_Message("Registered class name must not begin with '['"));

//	create help string
	std::string nameStr = std::string(nameIn);

//	check for non-allowed character
	size_t found = nameStr.find_first_not_of("abcdefghijklmnopqrstuvwxyz"
											 "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
											 "_"
											 "0123456789");
	if (found!=std::string::npos)
	{
	//	info: do not use UG_LOG here, since in parallel MPI_Init has not been
	//		  called, but UG_LOG checks proc rank.
		std::cout<<"Non-allowed character '"<<nameStr[found]<<"' "<<
		       "contained at position "<<int(found)<<" in registered Class Name "
		       "'"<<nameStr<<"'.\nClass names must match the regular expression: "
		       "[a-zA-Z_][a-zA-Z_0-9]*, \ni.e. only alphabetic characters, numbers "
		       " and '_' are allowed; no numbers at the beginning.\n";
		throw(REGISTRY_ERROR_Message("Class Name must only contain [a-zA-Z0-9]."));
	}

//	check that no number at the beginning
	 found = nameStr.find_first_of("0123456789");
	if (found!=std::string::npos && found == 0)
	{
	//	info: do not use UG_LOG here, since in parallel MPI_Init has not been
	//		  called, but UG_LOG checks proc rank.
		std::cout<<"Class Name "<<nameStr<<" starts with a number.\nThis is "
				" not allowed. Please change naming.\n";
		throw(REGISTRY_ERROR_Message("Class Name must not start with number."));
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
set_name(const char* name, const char* group, bool newName)
{
//	set own name
	set_name(name, group, newName);

//	add parent nodes
	m_ClassNameNode.add_base_class(ClassNameProvider<TParent1>::class_name_node());
}

template <typename TClass>
template <typename TParent1, typename TParent2>
void ClassNameProvider<TClass>::
set_name(const char* name, const char* group,bool newName)
{
//	set own name
	set_name(name, group, newName);

//	add parent nodes
	m_ClassNameNode.add_base_class(ClassNameProvider<TParent1>::class_name_node());
	m_ClassNameNode.add_base_class(ClassNameProvider<TParent2>::class_name_node());
}

template <typename TClass>
bool ClassNameProvider<TClass>::is_a(const char* parent, bool strict)
{
//	check if class is forward declared
	if(m_bForwardDeclared)
		throw(UGFatalError("ERROR in 'ClassNameProvider::is_a':"
					"Class must not be foreward declared to use is_a"));

//  strict comparison: must match this class name, parents are not allowed
	if(strict)
	{
	//  check pointer
		if(parent == name()) return true;

	//  compare strings
		if(strcmp(parent, name()) == 0) return true;

	//  names does not match
		return false;
	}

//	return if parent name is contained in tree of base classes
	return ClassNameTreeContains(m_ClassNameNode, parent);
}

template <typename TClass>
const char* ClassNameProvider<TClass>::name()
{
//	if name has not been set, set temporary forward declaration
	if(m_ClassNameNode.empty()) set_foreward_declared();

//	return the name of this class as stored in ClassNameNode
	return m_ClassNameNode.name().c_str();
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
	name.append(typeid(TClass).name());
	name.append(" (undeclared) ]]");

//	set this as current name, but remember pre-declaration
	m_ClassNameNode.set_name(name.c_str());
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

} // end namespace ug
} // end namespace bridge

#endif /* __H__UG_BRIDGE__CLASS_NAME_PROVIDER_IMPL__ */
