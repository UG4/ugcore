
#ifndef __H__UG_BRIDGE__CLASS_NAME_PROVIDER__
#define __H__UG_BRIDGE__CLASS_NAME_PROVIDER__

#include "common/assert.h"
#include <vector>
#include <cstring>
#include <string>
#include <iostream>
#include <typeinfo>

namespace ug
{
namespace bridge
{

struct UG_REGISTRY_ERROR_ClassUnknownToRegistry {};
struct UG_REGISTRY_ERROR_ClassAlreadyNamed
{
		UG_REGISTRY_ERROR_ClassAlreadyNamed(const char* name_)
			: name(name_)
		{}
		std::string name;
};


bool ClassNameVecContains(const std::vector<const char*>& names, const char* name);

template <typename TClass>
struct ClassNameProvider
{
	/// set name of class and copy parent names
		static void set_name(const char* nameIn, const char* group = "", bool newName = false)
		{
		//	if class already named throw error
			if(newName == true && bForwardDeclared==false && !m_ownName.empty())
			{
				if(strcmp(nameIn, name()) != 0)
					throw(UG_REGISTRY_ERROR_ClassAlreadyNamed(nameIn));
			}

		//	copy name into static string
		//	This is necessary, since char* could be to temporary memory
			m_ownName = std::string(nameIn);
			UG_ASSERT(m_ownName.size() > 0 && std::isalpha(m_ownName[0]), "name must be longer than 0");

		//	remember const char* to own name in first position of names-list
			m_names.clear();
			m_names.push_back(m_ownName.c_str());

		//	remember groups
			m_group = std::string(group);
		}

	/// set name of class and copy parent names
		template <typename TParent>
		static void set_name(const char* name, const char* group = "", bool newName = false)
		{
		//	set own name
			set_name(name, group, newName);

		//	copy parent names
			const std::vector<const char*>& pnames = ClassNameProvider<TParent>::names();

			for(size_t i = 0; i < pnames.size(); ++i)
				m_names.push_back(pnames[i]);
		}

		/** check if a given class name is equal
		 * This function checks if a given class name is equal to this class name.
		 * If the flag 'strict' is set to true, the class name must match exactly.
		 * If it is set to false, also parent names (of the class hierarchy) are checked
		 * and if this class inherits from the parent class true is returned.
		 */
		static bool is_a(const char* parent, bool strict = false)
		{
			UG_ASSERT(!bForwardDeclared, "Class must not be foreward declared to use is_a");
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

			return ClassNameVecContains(m_names, parent);
		}

		/// name of this class
		static const char* name()
		{
			if(m_ownName.empty())
			{
				//throw(UG_REGISTRY_ERROR_ClassUnknownToRegistry());
				set_foreward_declared();
			}
			return m_ownName.c_str();
		}

		static void set_foreward_declared()
		{
			m_ownName = "[[";
			m_ownName.append(typeid(TClass).name());
			m_ownName.append(" (undeclared) ]]");
			m_names.push_back(m_ownName.c_str());
			bForwardDeclared = true;
		}


		/// groups
		static const std::string& group()
		{
			return m_group;
		}

		/// returns vector of all names including parent class names
		static const std::vector<const char*>& names()	
		{
			if(m_names.empty())
			{
				//throw(UG_REGISTRY_ERROR_ClassUnknownToRegistry());
				set_foreward_declared();
			}

			return m_names;
		}

	private:
		static std::string m_ownName;
		static std::vector<const char*> m_names;
		static std::string m_group;
		static bool bForwardDeclared;
};

template <typename TClass>
std::vector<const char*> ClassNameProvider<TClass>::m_names = std::vector<const char*>();

template <typename TClass>
std::string ClassNameProvider<TClass>::m_ownName = std::string("");

template <typename TClass>
std::string ClassNameProvider<TClass>::m_group = std::string("");

template <typename TClass>
bool ClassNameProvider<TClass>::bForwardDeclared = false;

} // end namespace
} // end namespace

#endif /* __H__UG_BRIDGE__CLASS_NAME_PROVIDER__ */
