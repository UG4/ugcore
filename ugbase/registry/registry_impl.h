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
	bool success = funcGrp->add_overload(func, &FunctionProxy<TFunc>::apply,
	                                     methodOptions, group,
	                                     retValInfos, paramInfos,
	                                     tooltip, help);

	if(!success){
		UG_THROW_REGISTRY_ERROR(strippedMethodName,
		"Trying to register function name '"<<funcName
		<< "', that is already used by another function in this registry.");
	}

	return *this;
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
