
#ifndef __H__UG_BRIDGE__CLASS_NAME_PROVIDER__
#define __H__UG_BRIDGE__CLASS_NAME_PROVIDER__

#include "common/assert.h"
#include <vector>
#include <cstring>
#include <string>
#include <algorithm>
#include <iostream>
#include <typeinfo>
#include <map>
#include "common/common.h"

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


/// node for class names
class ClassNameNode
{
	public:
	///	set name
		void set_name(const char* name)
		{
			m_name = std::string(name);
			UG_ASSERT(m_name.size() > 0, "name must be longer than 0");
		}

	///	add a base class
		void add_base_class(const ClassNameNode& node)
		{
			std::vector<const ClassNameNode*>::iterator it
				= std::find(m_vBaseClass.begin(), m_vBaseClass.end(), &node);
			if(it == m_vBaseClass.end())
				m_vBaseClass.push_back(&node);
		}

	///	returns own name
		const std::string& name() const {return m_name;}

	///	returns if a name has been set
		bool empty() const {return m_name.empty();}

	///	returns number of parents
		size_t num_base_classes() const {return m_vBaseClass.size();}

	///	return a base class
		const ClassNameNode& base_class(size_t i) const {return *m_vBaseClass.at(i);}

	protected:
	///	own name
		std::string m_name;

	///	base classes
		std::vector<const ClassNameNode*> m_vBaseClass;
};

/// returns the vector containing all names in the name tree for node and its base classes
void ExtractClassNameVec(std::vector<const char*>& names,
                         const ClassNameNode& node,
                         bool clearVec = true);

///	returns if a name is contained in the name vector
bool ClassNameVecContains(const std::vector<const char*>& names, const char* name);

///	returns if a name is contained in the name tree at node or in base classes
bool ClassNameTreeContains(const ClassNameNode& node, const char* name);

/// returns an std::vector that contains in reverse order the base class that
///	must be used in the Class Hierarchy to get to the base class
bool ClassNameTreeWay(std::vector<size_t>& vWay, const ClassNameNode& node, const char* name);

template <typename TClass>
struct ClassNameProvider
{
	/// set name of class and copy parent names
		static void set_name(const char* nameIn, const char* group = "", bool newName = false)
		{
		//	if class already named throw error
			if(newName == true && bForwardDeclared==false && !m_ClassNameNode.empty())
			{
				if(strcmp(nameIn, name()) != 0)
					throw(UG_REGISTRY_ERROR_ClassAlreadyNamed(nameIn));
			}

		//	copy name into static string
			m_ClassNameNode.set_name(nameIn);

		//	remember groups
			m_group = std::string(group);

		//	create names list
		//	\todo: remove this, avoid names-vector
			ExtractClassNameVec(m_names, m_ClassNameNode, true);
		}

	/// set name of class and copy parent names
		template <typename TParent1>
		static void set_name(const char* name, const char* group = "", bool newName = false)
		{
		//	set own name
			set_name(name, group, newName);

		//	add parent nodes
			m_ClassNameNode.add_base_class(ClassNameProvider<TParent1>::class_name_node());

		//	create names list
		//	\todo: remove this, avoid names-vector
			ExtractClassNameVec(m_names, m_ClassNameNode, true);
		}

	/// set name of class and copy parent names
		template <typename TParent1, typename TParent2>
		static void set_name(const char* name, const char* group = "", bool newName = false)
		{
		//	set own name
			set_name(name, group, newName);

		//	add parent nodes
			m_ClassNameNode.add_base_class(ClassNameProvider<TParent1>::class_name_node());
			m_ClassNameNode.add_base_class(ClassNameProvider<TParent2>::class_name_node());

		//	create names list
		//	\todo: remove this, avoid names-vector
			ExtractClassNameVec(m_names, m_ClassNameNode, true);
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

			return ClassNameTreeContains(m_ClassNameNode, parent);
		}

	/// name of this class
		static const char* name()
		{
			if(m_ClassNameNode.empty())
				set_foreward_declared();

			return m_ClassNameNode.name().c_str();
		}

		static void set_foreward_declared()
		{
			std::string name("[[");
			name.append(typeid(TClass).name());
			name.append(" (undeclared) ]]");
			m_ClassNameNode.set_name(name.c_str());
			bForwardDeclared = true;

		//	create names list
		//	\todo: remove this, avoid names-vector
			ExtractClassNameVec(m_names, m_ClassNameNode, true);
		}


	/// groups
		static const std::string& group(){return m_group;}

	/// returns vector of all names including parent class names
		static const std::vector<const char*>& names()	
		{
			if(m_names.empty())
				set_foreward_declared();

			return m_names;
		}

	///	return the class name node in the class hierarchy
		static const ClassNameNode& class_name_node() {return m_ClassNameNode;}

	private:
		static std::string m_ownName;
		static std::vector<const char*> m_names;
		static std::string m_group;
		static bool bForwardDeclared;

		static ClassNameNode m_ClassNameNode;
};

template <typename TClass>
std::vector<const char*> ClassNameProvider<TClass>::m_names = std::vector<const char*>();

template <typename TClass>
std::string ClassNameProvider<TClass>::m_ownName = std::string("");

template <typename TClass>
std::string ClassNameProvider<TClass>::m_group = std::string("");

template <typename TClass>
bool ClassNameProvider<TClass>::bForwardDeclared = false;

template <typename TClass>
ClassNameNode ClassNameProvider<TClass>::m_ClassNameNode = ClassNameNode();

// 	static cast function for two classes
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

struct ClassCastProvider
{
	public:
	//	add a cast function to the registry: Casts: Derived -> Base
		template <typename TBase, typename TDerived>
		static void add_cast_func()
		{
			const ClassNameNode* pBaseNode = &ClassNameProvider<TBase>::class_name_node();
			const ClassNameNode* pDerivNode = &ClassNameProvider<TDerived>::class_name_node();

			std::pair<const ClassNameNode*, const ClassNameNode*> namePair(pBaseNode, pDerivNode);

			m_mmCast[namePair] = &StaticVoidCast<TBase, TDerived>;
		}

	///	cast a pointer to the desired base class
	/**
	 * This method casts a void pointer to a given derived class to the void
	 * pointer of a base class.
	 *
	 * \param[in]	pDerivVoid		void pointer to Derived object
	 * \param[in,out]	node		on entry: class name node corresponding to pDerivVoid
	 * 								on exit:  class name node corresponding to baseName
	 * \param[in]		baseName	name of base class the pointer should be casted to
	 * \returns		void* to base class, NULL if cast not possible
	 */
		static void* cast_to_base_class(void* pDerivVoid, const ClassNameNode*& node, const char* baseName)
		{
		//	find way to base class
			std::vector<size_t> vWay;
			if(!ClassNameTreeWay(vWay, *node, baseName))
			{
				UG_LOG("ERROR in ClassCastProvider::cast_to_base_class: Request"
						" to cast from derived class '"<< node->name()<<"' to "
						" base class '"<<baseName<<"', but no such base class in"
						" registered class hierarchy.\n");
				return NULL;
			}

			void* currPtr = pDerivVoid;
			const ClassNameNode* pCurrNode = node;

		//	cast all the way down
			while(!vWay.empty())
			{
			//	get base class to cast to
				const ClassNameNode* pBaseClassNode = &pCurrNode->base_class(vWay.back());

			//	get name pair
				std::pair<const ClassNameNode*, const ClassNameNode*> namePair(pBaseClassNode, pCurrNode);

			//	find in map
				std::map<std::pair<const ClassNameNode*, const ClassNameNode*>, CastFunc>::iterator it;
				it = m_mmCast.find(namePair);

				if(it == m_mmCast.end())
				{
					UG_LOG("ERROR in ClassCastProvider::cast_to_base_class:"
							" Request intermediate cast from derived class '" <<
							pCurrNode->name() <<"' to direct base class '"
							<<pBaseClassNode->name()<<"', but no such cast "
							" function registered.\n");
					return NULL;
				}

			//	get cast function
				CastFunc pCastFunc = it->second;

			//	cast
				currPtr = (*pCastFunc)(currPtr);

			//	set node to base class
				pCurrNode = pBaseClassNode;

			//	pop way
				vWay.pop_back();
			}

		//	write current node on exit
			node = pCurrNode;

		//	return current pointer
			return currPtr;
		}

	protected:
	//	type of cast pointer
		typedef void* (*CastFunc)(void*);

	//	cast map
		static std::map<std::pair<const ClassNameNode*, const ClassNameNode*>, CastFunc> m_mmCast;
};



} // end namespace
} // end namespace

#endif /* __H__UG_BRIDGE__CLASS_NAME_PROVIDER__ */
