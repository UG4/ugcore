/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG_BRIDGE__REGISTRY__
#define __H__UG_BRIDGE__REGISTRY__

#include <vector>
#include <string>
#include <cstring>
#include <typeinfo>
#include <iostream>
#include <boost/function.hpp>
#include <boost/type_traits.hpp>


#include "global_function.h"
#include "class.h"
#include "param_to_type_value_list.h"
#include "parameter_stack.h"
#include "common/ug_config.h"


namespace ug
{
namespace bridge
{

/**
 * \defgroup registry Registry
 * \ingroup bridge
 * \{
 */

//	PREDECLARATIONS
class Registry;

///	declaration of registry callback function.
/**	Allows to notify listeners if the registry changes.
 * Since FuncRegistryChanged is a functor, you can either
 * pass a normal function or a member function of a class
 * (Have a look at boost::bind in the second case).
 */
typedef boost::function<void (Registry* pReg)> FuncRegistryChanged;


///	groups classes. One of the members is the default member.
class UG_API ClassGroupDesc
{
	public:
		ClassGroupDesc() : m_defaultClass(NULL)	{}

	///	sets name of group
		void set_name(const std::string& name) 	{m_name = name;}

	///	returns name of group
		const std::string& name() const			{return m_name;}

	///	adds a class to group
		void add_class(IExportedClass* c, const std::string& tag)
			{m_classes.push_back(c); m_classTags.push_back(tag);}

	///	returns number of classes in group
		size_t num_classes() const			{return m_classes.size();}

	///	returns if classes in group
		bool empty() const					{return num_classes() == 0;}

	///	returns a class of the group
		IExportedClass* get_class(size_t i)	{return m_classes[i];}

	///	returns a class of the group
		const IExportedClass* get_class(size_t i) const	{return m_classes[i];}

	///	returns the class group tag for a class
		const std::string& get_class_tag(size_t i) const{return m_classTags[i];}

	///	sets the i'th class as default
		void set_default_class(size_t i)	{m_defaultClass = m_classes[i];}

	///	if no default class is set, this method returns NULL.
		IExportedClass* get_default_class()	const {return m_defaultClass;}

	private:
	///	name of class group
		std::string						m_name;

	///	classes registered to the class group
		std::vector<IExportedClass*>	m_classes;

	/// tags can be used to describe classes. One tag for each class.
		std::vector<std::string>		m_classTags;

	///	the current default class
		IExportedClass* 				m_defaultClass;
};


/// Registry for functions and classes that are exported to scripts and visualizations
/**
 * It also allows to register callbacks that are called if the registry changes.
 *
 * Please note that once a class or method is registered at the registry, it
 * can not be removed (This is important for the implementation of callbacks).
 *
 * Please note, that some of the methods take std::string by copy instead of
 * using const std::string&. This is on purpose in order to allow a call
 * of the method using a const char* as well.
 * @sa \ref pageUG4Registry
 */
class UG_API Registry {
	public:
	///	constructor
		Registry();

	///	destructor
		~Registry();
		
	///	requires every constructable class to be constructed via Smart-Pointer
		void set_force_construct_via_smart_pointer(bool bForceConstructionWithSmartPtr);

	////////////////////////
	//	callbacks
	////////////////////////

	///	adds a callback which is triggered whenever Registry::registry_changed is called.
		void add_callback(FuncRegistryChanged callback);

	///	call this method if to forward changes of the registry to its listeners
		bool registry_changed();

	//////////////////////
	// global functions
	//////////////////////
		
	/**
	 * @brief adds a function to the registry
	 * @param funcName the name of the function
	 * @param func function pointer of the function
	 * @param group registry group. use / for subgroups e.g. ug4/mygroup/mysubgroup (optional)
	 * @param retValInfos string documenting what the function returns (optional)
	 * @param paramInfos string documenting the parameters of the function
	 * seperate parameters with an # e.g. "x#y#z" (don't specify the type of the values)  (optional)
	 * @param toolTip small documentation for the function (optional)
	 * @param help help string for the function
	 * @sa \ref pageUG4Registry
	 * @sa \ref secSTHowToSpecifyParameterInformation
	 *
	 * References the template function proxy_function<TFunc> and stores
	 * it with the FuntionWrapper.
	 */
		template<typename TFunc>
		Registry& add_function(std::string funcName, TFunc func, std::string group = "",
		                       std::string retValInfos = "", std::string paramInfos = "",
		                       std::string tooltip = "", std::string help = "");

	///	Similar to add_function but returns the added function instead of the registry.
		template<typename TFunc>
		ExportedFunction*
		add_and_get_function(std::string funcName, TFunc func, std::string group = "",
		                     std::string retValInfos = "", std::string paramInfos = "",
		                     std::string tooltip = "", std::string help = "");

	/// number of functions registered at the Registry (overloads are not counted)
		size_t num_functions() const;

	/// returns the first overload of an exported function
		ExportedFunction& get_function(size_t ind);
		const ExportedFunction& get_function(size_t ind) const;

	///	returns the number of overloads of a function
		size_t num_overloads(size_t ind) const;

	///	returns the i-th overload of a function
		ExportedFunction& get_overload(size_t funcInd, size_t oInd);

	///	returns a group which contains all overloads of a function
		ExportedFunctionGroup& get_function_group(size_t ind);

	///	returns an exported function group by name
		ExportedFunctionGroup* get_exported_function_group(const std::string& name);

	///////////////////
	// classes
	///////////////////

	/**
	 * @brief Register a class at this registry
	 * @param className name of the class to appear in the registry
	 * @param group registry group. use / for subgroups e.g. ug4/mygroup/mysubgroup (optional)
	 * @param toolTip describing text for the class (optional)
	 */
		template <typename TClass>
		ExportedClass<TClass>& add_class_(std::string className,
		                                   std::string group = "",
		                                   std::string tooltip = "");

	/**
	 * @brief Register a class at this registry together with its base class
	 * @param className name of the class to appear in the registry
	 * @param group registry group. use / for subgroups e.g. ug4/mygroup/mysubgroup (optional)
	 * @param toolTip describing text for the class (optional)
	 */
		template <typename TClass, typename TBaseClass>
		ExportedClass<TClass>& add_class_(std::string className,
		                                   std::string group = "",
		                                   std::string tooltip = "");

	/**
	 * @brief Register a class at this registry together with two base classes
	 * @param className name of the class to appear in the registry
	 * @param group registry group. use / for subgroups e.g. ug4/mygroup/mysubgroup (optional)
	 * @param toolTip describing text for the class (optional)
	 */
		template <typename TClass, typename TBaseClass1, typename TBaseClass2>
		ExportedClass<TClass>& add_class_(std::string className,
		                                   std::string group = "",
		                                   std::string tooltip = "");

	/// Get Reference to already registered class
		template <typename TClass>
		ExportedClass<TClass>& get_class_();

	/// number of classes registered at the Registry
		size_t num_classes() const;

	/// returns an exported class
		const IExportedClass& get_class(size_t ind)	const;

	/// returns an exported class
		IExportedClass* get_class(const std::string& name);

	/// returns an exported class
		const IExportedClass* get_class(const std::string& name) const;

	///	returns true if everything well-declared, else false
		bool check_consistency();


	///////////////////
	// class-groups
	///////////////////

	///	returns the number of available class groups
		size_t num_class_groups() const;

	///	returns a const pointer to the i-th class group
		const ClassGroupDesc* get_class_group(size_t i) const;

    ///	returns a pointer to the i-th class group
		ClassGroupDesc* get_class_group(size_t i);

	///	Returns the class-group with the given name.
	/**	If no such group exists at the time of calling, it will be created.*/
		ClassGroupDesc* get_class_group(const std::string& name);

	///	Returns the class-group with the given name.
	/**	If no such group exists at the time of calling, NULL is returned.*/
		const ClassGroupDesc* get_class_group(const std::string& name) const;

	///	adds the given class to the given group.
	/**	Groups are constructed automatically if required.
	 * This method is just for convenience. It is effectively the same as:
	 * get_class_group(groupName).add_class(reg.get_class(className), classTag).*/
		void add_class_to_group(std::string className, std::string groupName,
		                        std::string classTag = "");


	/// returns true if classname is already used by a class in this registry
		bool classname_registered(const std::string& name);

	/// returns true if groupname is already used by a class in this registry
		bool groupname_registered(const std::string& name);

	/// returns true if functionname is already used by a function in this registry
		bool functionname_registered(const std::string& name);

	protected:
	///	performs some checks, throws error if something wrong
		template <typename TClass, typename TBaseClass>
		void check_base_class(const std::string& className);

	private:
	//	disallow copy
		Registry(const Registry& reg);
		
	///	registered functions
		std::vector<ExportedFunctionGroup*>	m_vFunction;

	///	registered classes
		std::vector<IExportedClass*> m_vClass;

	///	registered class groups
		std::vector<ClassGroupDesc*> m_vClassGroups;

	///	Callback, that are called when registry changed is invoked
		std::vector<FuncRegistryChanged> m_callbacksRegChanged;

	///	flag if classes must be constructed via smart-pointer
		bool m_bForceConstructionWithSmartPtr;
};

// end group registry
/// \}

} // end namespace registry

} // end namespace ug




// Auxiliary function.
#define UG_REGISTRY_CONCAT(first, second) first##second

#ifdef UG_USE_PYBIND11

// Use this function to declare.
#define UG_REGISTRY_DECL(fname, ...)\
		namespace bridge { void fname (__VA_ARGS__); }\
		namespace pybind { void fname (__VA_ARGS__); }

#define UG_REGISTRY_DEFINE(fname)\
		namespace bridge {\
			void fname (ug::bridge::Registry& r, std::string grp)\
			{ UG_REGISTRY_CONCAT(fname,_) (r, grp); }\
		}\
		namespace pybind {\
			void fname (ug::pybind::RegistryAdapter& r, std::string grp)\
			{ ug::bridge::UG_REGISTRY_CONCAT(fname,_) (r, grp); }\
		}
#else

#define UG_REGISTRY_DECL(fname, ...)\
		namespace bridge { void fname (__VA_ARGS__); }

#define UG_REGISTRY_DEFINE(fname)\
		namespace bridge {\
			void fname (ug::bridge::Registry& r, std::string grp)\
			{ UG_REGISTRY_CONCAT(fname,_) (r, grp); }\
		}

#endif

////////////////////////////////
//	include implementation
#include "registry_impl.h"

#endif /* __H__UG_BRIDGE__REGISTRY__ */
