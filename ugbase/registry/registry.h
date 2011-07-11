//	Authors: Andreas Vogel, Sebastian Reiter

#ifndef __H__UG_BRIDGE__REGISTRY__
#define __H__UG_BRIDGE__REGISTRY__

#include <vector>
#include <string>
#include <map>
#include <cstring>
#include <typeinfo>
#include <iostream>
#include <boost/function.hpp>
#include <boost/type_traits.hpp>


#include "global_function.h"
#include "class.h"
#include "param_to_type_value_list.h"
#include "parameter_stack.h"

namespace ug
{
namespace bridge
{

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
class ClassGroupDesc
{
	public:
		ClassGroupDesc() : m_defaultClass(NULL)	{}

		void set_name(const char* name)		{m_name = name;}
		const char* name() const			{return m_name.c_str();}

		void add_class(IExportedClass* c)	{m_classes.push_back(c);}
		size_t num_classes() const			{return m_classes.size();}
		bool empty() const					{return num_classes() != 0;}
		IExportedClass* get_class(size_t i)	{return m_classes[i];}
		const IExportedClass* get_class(size_t i) const	{return m_classes[i];}

		void set_default_class(size_t i)	{m_defaultClass = m_classes[i];}

	///	if no default class is set, this method returns NULL.�
		IExportedClass* get_default_class()	const {return m_defaultClass;}

	private:
		std::string						m_name;
		std::vector<IExportedClass*>	m_classes;
		IExportedClass* 				m_defaultClass;
};


// Registry
/** registers functions and classes that are exported to scripts and visualizations
 * It also allows to register callbacks that are called if the
 * registry changes.
 *
 * Please note that once a class or method is registered at the
 * registry, it will can not be removed (This is important for the
 * implementation of callbacks).
 */
class Registry {
	public:
		Registry();
		~Registry();
		
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
		
	/**	References the template function proxy_function<TFunc> and stores
	 * it with the FuntionWrapper.
	 */
		template<class TFunc>
		Registry& add_function(const char* funcName, TFunc func, const char* group = "",
								const char* retValInfos = "", const char* paramInfos = "",
								const char* tooltip = "", const char* help = "");

	/// number of functions registered at the Registry (overloads are not counted)
		size_t num_functions() const;

	/// returns the first overload of an exported function
		ExportedFunction& get_function(size_t ind);

	///	returns the number of overloads of a function
		size_t num_overloads(size_t ind);

	///	returns the i-th overload of a function
		ExportedFunction& get_overload(size_t funcInd, size_t oInd);

	///	returns a group which contains all overloads of a function
		ExportedFunctionGroup& get_function_group(size_t ind);

	///////////////////
	// classes
	///////////////////

	/** Register a class at this registry
	 * This function registers any class
	 */
		template <typename TClass>
		ExportedClass_<TClass>& add_class_(const char* className,
											const char* group = "",
											const char *tooltip="");

	/** Register a class at this registry
	 * This function registers any class together with its base class
	 */
		template <typename TClass, typename TBaseClass>
		ExportedClass_<TClass>& add_class_(const char* className,
											const char* group = "",
											const char *tooltip = "");

	/** Register a class at this registry
	 * This function registers any class together with its base class
	 */
		template <typename TClass, typename TBaseClass1, typename TBaseClass2>
		ExportedClass_<TClass>& add_class_(const char* className,
											const char* group = "");

	/**
	 * Get Reference to already registered class
	 */
		template <typename TClass>
		ExportedClass_<TClass>& get_class_();

	/// number of classes registered at the Registry
		size_t num_classes() const;

	/// returns an exported class
		const IExportedClass& get_class(size_t ind)	const;

	/// returns an exported class
		IExportedClass* get_class(const char* name);

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
		ClassGroupDesc* get_class_group(const char* name);

	///	Returns the class-group with the given name.
	/**	If no such group exists at the time of calling, NULL is returned.*/
		const ClassGroupDesc* get_class_group(const char* name) const;

	///	adds the given class to the given group.
	/**	Groups are constructed automatically if required.
	 * This method is just for conveniance. It is effectively the same as:
	 * get_class_group(groupName).add_class(reg.get_class(className))*/
		void add_class_to_group(const char* className, const char* groupName);

	protected:
	///	performs some checks, throws error if something wrong
		template <typename TClass, typename TBaseClass>
		void check_base_class(const char* className);

	// returns true if classname is already used by a class in this registry
		bool classname_registered(const char* name);

		// returns true if functionname is already used by a function in this registry
		bool functionname_registered(const char* name);

		ExportedFunctionGroup* get_exported_function_group(const char* name);

	private:
	//	disallow copy
		Registry(const Registry& reg);
		
	//	registered functions
		std::vector<ExportedFunctionGroup*>	m_vFunction;

	//	registered classes
		std::vector<IExportedClass*> m_vClass;

	//	class groups
		std::vector<ClassGroupDesc*> m_vClassGroups;

	//	Callback, that are called when registry changed is invoked
		std::vector<FuncRegistryChanged> m_callbacksRegChanged;
};

} // end namespace registry

} // end namespace ug

////////////////////////////////
//	include implementation
#include "registry_impl.h"

#endif /* __H__UG_BRIDGE__REGISTRY__ */
