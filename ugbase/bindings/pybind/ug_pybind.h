#pragma once

// pybind11 lib.
#include <pybind11/pybind11.h>


// UG4 lib.
#include "registry/registry.h"
#include "registry/class.h"
#include "registry/function_traits.h"

// #include "bridge/util.h"
//#include "bridge/util_domain_dependent.h"

// #include "bridge/util_domain_algebra_dependent.h"


// #define UNROLL_TYPE_LIST(TL) ((TL::length>0)? TL::head, UNROLL_TYPE_LIST(TL:tail) :  TL::head )


//! Use SmartPtr
//PYBIND11_DECLARE_HOLDER_TYPE(T, SmartPtr<T>, true);

PYBIND11_DECLARE_HOLDER_TYPE(T, SmartPtr<T>);

namespace ug{

namespace pybind {

//! Adapter for exported classes.
/*! Maps to pybind11::class_.
 * 	May be defined using variadic templates. */


template <typename ...TClass>
struct ExportedClass : public pybind11::class_<TClass...>
{

	typedef pybind11::class_<TClass...> base_type;
	typedef ExportedClass<TClass...> exported_class_type;

	//! Construct from pybind class
	ExportedClass(const base_type &py)
	: base_type(py) {}



public:
	// Constructor (w/o arguments).
	exported_class_type &add_constructor()
	{
		base_type::def( pybind11::init<>() ); // no argument
		return (*this);
	}

protected:
	// Use function specialization (instead of partial templates)
	template <typename T> struct type{};

	template <typename TRet, class... Ts>
	void variadic_def(type <TRet (*) (Ts...)>)
	{ base_type::def( pybind11::init<Ts...>() ); }; // template arguments

	template <typename TRet, typename TNew, class... Ts>
	void variadic_def(type <TRet (TNew::*) (Ts...)>)
	{ base_type::def( pybind11::init<Ts...>() ); }; // template arguments

public:
	// Constructor (w/ template arguments).
	template<typename TFunc>
	exported_class_type &add_constructor(std::string paramInfos = "",
            std::string tooltip = "", std::string help = "",
            std::string options = "")
	{
		variadic_def(type<TFunc>{}); // extract arguments from signature
		return (*this);
	}

	//! Add a method to class.
	template<typename TMethod>
	exported_class_type &add_method (std::string methodName, TMethod func,
			                     std::string retValInfos = "", std::string paramInfos = "",
			                     std::string tooltip = "", std::string help = "")
	{
		base_type::def(methodName.c_str(), func);
		return (*this);
	}

	//! Always use smart pointers.
	void set_construct_as_smart_pointer (bool flag)
	{};

	void construct_as_smart_pointer ()
	{ set_construct_as_smart_pointer(true); }


};

//! This adapter provides access to a pybind11::module as a UG registry conforming object.
struct RegistryAdapter : public pybind11::module
{
	//! CTOR (wraps 'module')
	RegistryAdapter(pybind11::module &py)
	: pybind11::module(py) {}

	//! Add a function.
	template <typename TFunc>
	RegistryAdapter &add_function(std::string name, TFunc f, std::string grp="",
						std::string retValInfos = "", std::string paramInfos = "",
						std::string tooltip = "", std::string help = "")
	{
		pybind11::module::def(name.c_str(), f, tooltip.c_str());
		return (*this);
	}

	//! Add a class.
	template <typename TClass>
	ExportedClass<TClass> add_class_(std::string className,
			                                   std::string group = "",
			                                   std::string tooltip = "")
	{
		typedef pybind11::class_<TClass> pyclass;
		return pyclass(*this, className.c_str());
	}

	template <typename TClass, typename TBaseClass1>
	ExportedClass<TClass, TBaseClass1> add_class_(std::string className,
			                                   std::string group = "",
			                                   std::string tooltip = "")
	{
		typedef pybind11::class_<TClass,TBaseClass1> pyclass;
		ExportedClass<TClass,TBaseClass1> myclass(pyclass(*this, className.c_str()));
		return myclass;

	}

	template <typename TClass, typename TBaseClass1, typename TBaseClass2>
	ExportedClass<TClass, TBaseClass1, TBaseClass2> add_class_(std::string className,
			                                   std::string group = "",
			                                   std::string tooltip = "")
	{
		typedef pybind11::class_<TClass,TBaseClass1,TBaseClass2> pyclass;
		return pyclass(*this, className.c_str());
	}

	//! TODO: Does this make sense?
	void add_class_to_group(std::string group, std::string className, std::string tag)
	{

	}

	void registry_changed()
	{};

};




} // namespace pybind
} // namespace ug
