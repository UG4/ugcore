/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
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


#ifndef __H__UG_BRIDGE__CLASS_NAME_PROVIDER__
#define __H__UG_BRIDGE__CLASS_NAME_PROVIDER__

#include "common/assert.h"
#include <vector>
#include <cstring>
#include <string>
#include <algorithm>
//#include <iostream>
//#include <typeinfo>
#include <map>
#include "common/common.h"
#include "common/util/smart_pointer.h"
#include "common/ug_config.h"
#include "error.h"

namespace ug {
namespace bridge {

/// \addtogroup registry
/// \{

/// node for class names
/**
 * A ClassNameNode stores the name of a registered class and pointers to
 * the ClassNameNodes of the direct base classes of this class. By traversing
 * the tree of ClassNameNodes all parent classes of a class can be found.
 */
class UG_API ClassNameNode
{
	public:
	///	constructor
		ClassNameNode();

	///	set name
		void set_name(const std::string& name);

	///	add a base class
		void add_base_class(const ClassNameNode& node);

	///	returns own name
		[[nodiscard]] const std::string& name() const {return m_name;}

	///	returns if a name has been set
		[[nodiscard]] bool empty() const {return m_name.empty();}

	///	returns if a name has been set by the user
		[[nodiscard]] bool named() const;

	///	returns number of parents
		[[nodiscard]] size_t num_base_classes() const {return m_vBaseClass.size();}

	///	return a base class
		[[nodiscard]] const ClassNameNode& base_class(size_t i) const {return *m_vBaseClass[i];}

	protected:
	///	own name
		std::string m_name;

	///	base classes
		std::vector<const ClassNameNode*> m_vBaseClass;
};

///	provides the name for a class
template <typename TClass>
class UG_API ClassNameProvider
{
	public:
	/// set name of class and copy parent names
		static void set_name(const std::string& nameIn, const std::string& group,
		                     bool newName = false);

	/// set name of class and copy parent names
		template <typename TParent1>
		static void set_name(const std::string& nameIn, const std::string& group,
		                     bool newName = false);

	/// set name of class and copy parent names
		template <typename TParent1, typename TParent2>
		static void set_name(const std::string& nameIn, const std::string& group,
		                     bool newName = false);

	/** check if a given class name is equal
	 * This function checks if a given class name is equal to this class name.
	 * If the flag 'strict' is set to true, the class name must match exactly.
	 * If it is set to false, also parent names (of the class hierarchy) are checked
	 * and if this class inherits from the parent class true is returned.
	 */
		static bool is_a(const std::string& parent, bool strict = false);

	/// name of this class
		static const std::string& name();

	/// groups
		static const std::string& group(){return m_group;}

	/// returns vector of all names including parent class names
		static const std::vector<const char*>& names();

	///	return the class name node in the class hierarchy
		static const ClassNameNode& class_name_node() {return m_ClassNameNode;}

	///	returns if class name is forward declared
		static bool forward_declared() {return m_bForwardDeclared;}

	///	returns if the class has been named by user
		static bool named() {return !m_bForwardDeclared && !m_ClassNameNode.empty();}

	protected:
	///	sets a temporary name to the class
		static void set_forward_declared();

	private:
	///	vector of parent class names (depreciated)
		static std::vector<const char*> m_names;

	///	Name of group, we're this class is sorted
		static std::string m_group;

	///	set to true if class has not been named by user, but a default name given
		static bool m_bForwardDeclared;

	///	class name node holding own name and pointers to base class name nodes
		static ClassNameNode m_ClassNameNode;
};

template <typename TClass>
const char* GetClassName(){
	return ClassNameProvider<TClass>::name().c_str();
}

///	static cast function for two classes
template <typename TBase, typename TDerived>
void* StaticVoidCast(void* DerivVoidPtr);


struct UGError_ClassCastFailed : UGError{
	UGError_ClassCastFailed(const std::string& from, const std::string& to) :
		UGError("Class cast failed"), m_from(from), m_to(to)	{}

	std::string m_from;
	std::string m_to;
};

///	provides castings from derived classes to base classes
class ClassCastProvider
{
	public:
	///	add a cast function to the registry: Casts: Derived -> Base
		template <typename TBase, typename TDerived>
		static void add_cast_func();

	///	cast a pointer to the desired base class
	/**
	 * This method casts a void pointer to a given derived class to the void
	 * pointer of a base class. If conversion fails, an exception of type
	 * UGError_ClassCastFailed is thrown.
	 *
	 * \param[in]	pDerivVoid		void pointer to Derived object
	 * \param[in,out]	node		on entry: class name node corresponding to pDerivVoid
	 * 								on exit:  class name node corresponding to baseName
	 * \param[in]		baseName	name of base class the pointer should be casted to
	 * \returns		void* to base class
	 */
	/// \{
		static void* cast_to_base_class(void* pDerivVoid,
		                                const ClassNameNode*& node,
		                                const std::string& baseName);
		static const void* cast_to_base_class(const void* pDerivVoid,
		                                      const ClassNameNode*& node,
		                                      const std::string& baseName);
	//// \}

	///	casts a void pointer to a concrete class
	/**
	 * This method casts a void pointer to a given derived classed and returns it
	 * as a reinterpreted cast to the type specified by the template argument.
	 * If conversion fails, an exception of type UGError_ClassCastFailed is thrown.
	 *
	 * \param[in]	ptr		void pointer to Derived object
	 * \param[in,out]	node		on entry: class name node corresponding to pDerivVoid
	 * 								on exit:  class name node corresponding to baseName
	 * \returns		void* to base class
	 */
	/// \{
		template <typename T>
		static T* cast_to(void* ptr, const ClassNameNode*& node);
		template <typename T>
		static const T* cast_to(const void* ptr, const ClassNameNode*& node);
		template <typename T>
		static SmartPtr<T> cast_to(SmartPtr<void> ptr, const ClassNameNode*& node);
		template <typename T>
		static ConstSmartPtr<T> cast_to(ConstSmartPtr<void> ptr, const ClassNameNode*& node);
	///	\}

	protected:
	//	type of cast pointer
		using CastFunc = void* (*)(void*);

	//	cast map
		static std::map<std::pair<const ClassNameNode*, const ClassNameNode*>, CastFunc> m_mmCast;
};


/// returns the vector containing all names in the name tree for node and its base classes
void ExtractClassNameVec(std::vector<const char*>& names,
                         const ClassNameNode& node,
                         bool clearVec = true);

///	returns if a name is contained in the name vector
bool ClassNameVecContains(const std::vector<const char*>& names, const std::string& name);

///	returns if a name is contained in the name tree at node or in base classes
bool ClassNameTreeContains(const ClassNameNode& node, const std::string& name);

/// returns an std::vector that contains in reverse order the base class that
///	must be used in the Class Hierarchy to get to the base class
bool ClassNameTreeWay(std::vector<size_t>& vWay, const ClassNameNode& node, const std::string& name);

// end group registry
/// \}

} // end namespace
} // end namespace

// include implementation
#include "class_name_provider_impl.h"

#endif