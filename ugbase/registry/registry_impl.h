// created by Andreas Vogel, Sebastian Reiter
// s.b.reiter@googlemail.com
// 08.07.2011 (m,d,y)

#ifndef __H__UG_BRIDGE__REGISTRY_IMPL__
#define __H__UG_BRIDGE__REGISTRY_IMPL__

#include "common/util/string_util.h"

namespace ug{
namespace bridge
{

//////////////////////
// global functions
//////////////////////

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
		UG_LOG(strippedMethodName << " ... | ... " << methodOptions << std::endl);
	}

//	trim whitespaces
	strippedMethodName = TrimString(strippedMethodName);
	methodOptions = TrimString(methodOptions);

// 	check that name is not empty
	if(strippedMethodName.empty())
	{
		UG_LOG("### Registry ERROR: Trying to register empty function name."
				<< "\n### Please change register process. Aborting ...\n");
		throw(UG_REGISTRY_ERROR_RegistrationFailed(strippedMethodName));
	}
	
	// check that name does not contain illegal characters
	if (!IdentifierIsValid(strippedMethodName)) {
		UG_LOG("### Registry ERROR: Trying to register function '" 
				<< strippedMethodName << "' that"
				<< " contains illegal characters.\n"
				<< GetIdentifierMessage()
				<< "\n### Please change register process. Aborting ...\n");
		throw(UG_REGISTRY_ERROR_RegistrationFailed(strippedMethodName));
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
		UG_LOG("### Registry ERROR: Trying to register function name '"<<funcName
				<< "', that is already used by another function in this registry."
				<< "\n### Please change register process. Aborting ...\n");
		throw(UG_REGISTRY_ERROR_RegistrationFailed(strippedMethodName));
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
		UG_LOG("### Registry ERROR: Trying to register class name '"<<className
				<< "', that is already used by another class in this registry."
				<< "\n### Please change register process. Aborting ...\n");
		throw(UG_REGISTRY_ERROR_RegistrationFailed(className));
	}
// 	check that name is not empty
	if(className.empty())
	{
		UG_LOG("### Registry ERROR: Trying to register empty class name."
				<< "\n### Please change register process. Aborting ...\n");
		throw(UG_REGISTRY_ERROR_RegistrationFailed(className));
	}

// 	check that base class is not same type as class
	if(typeid(TClass) == typeid(TBaseClass))
	{
		UG_LOG("### Registry ERROR: Trying to register class "<<className
				<< "\n### that derives from itself. Aborting ...\n");
		throw(UG_REGISTRY_ERROR_RegistrationFailed(className));
	}

// 	check that class derives from base class
	if(boost::is_base_of<TBaseClass, TClass>::value == false)
	{
		UG_LOG("### Registry ERROR: Trying to register class "<<className
				<< "\n### with base class that is no base class. Aborting ...\n");
		throw(UG_REGISTRY_ERROR_RegistrationFailed(className));
	}
}

template <typename TClass>
ExportedClass<TClass>& Registry::
add_class_(std::string className, std::string group, std::string tooltip)
{
//	check that className is not already used
	if(classname_registered(className))
	{
		UG_LOG("### Registry ERROR: Trying to register class name '"<<className
				<< "', that is already used by another class in this registry."
				<< "\n### Please change register process. Aborting ...\n");
		throw(UG_REGISTRY_ERROR_RegistrationFailed(className));
	}
//	check that className is not already used as a group name
	if(groupname_registered(className))
	{
		UG_LOG("### Registry ERROR: Trying to register class name '"<<className
				<< "', that is already used by another group in this registry."
				<< "\n### Please change register process. Aborting ...\n");
		throw(UG_REGISTRY_ERROR_RegistrationFailed(className));
	}
// 	check that name is not empty
	if(className.empty())
	{
		UG_LOG("### Registry ERROR: Trying to register empty class name."
				<< "\n### Please change register process. Aborting ...\n");
		throw(UG_REGISTRY_ERROR_RegistrationFailed(className));
	}
	
	// check that name does not contain illegal characters
	if (!IdentifierIsValid(className)) {
		UG_LOG("### Registry ERROR: Trying to register class '" 
				<< className << "' that"
				<< " contains illegal characters.\n"
				<< GetIdentifierMessage()
				<< "\n### Please change register process. Aborting ...\n");
		throw(UG_REGISTRY_ERROR_RegistrationFailed(className));
	}

//	new class pointer
	ExportedClass<TClass>* newClass = NULL;

//	try creation
	try
	{
		newClass = new ExportedClass<TClass>(className, group, tooltip);
	}
	catch(ug::bridge::REGISTRY_ERROR_ClassAlreadyNamed ex)
	{
		UG_LOG("### Registry ERROR: Trying to register class with name '"<<className
				<< "', that has already been named. This is not allowed. "
				<< "\n### Please change register process. Aborting ...\n");
		throw(UG_REGISTRY_ERROR_RegistrationFailed(className));
	}
	catch(ug::bridge::REGISTRY_ERROR_Message ex)
	{
		UG_LOG("### Registry ERROR: " << ex.msg << "\n");
		throw(UG_REGISTRY_ERROR_RegistrationFailed("(unknown)"));
	}

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
		UG_LOG("### Registry ERROR: Trying to register class name '"<<className
				<< "', that is already used by another class in this registry."
				<< "\n### Please change register process. Aborting ...\n");
		throw(UG_REGISTRY_ERROR_RegistrationFailed(className));
	}
//	check that className is not already used as a group name
	if(groupname_registered(className))
	{
		UG_LOG("### Registry ERROR: Trying to register class name '"<<className
				<< "', that is already used by another group in this registry."
				<< "\n### Please change register process. Aborting ...\n");
		throw(UG_REGISTRY_ERROR_RegistrationFailed(className));
	}
// 	check that name is not empty
	if(className.empty())
	{
		UG_LOG("### Registry ERROR: Trying to register empty class name."
				<< "\n### Please change register process. Aborting ...\n");
		throw(UG_REGISTRY_ERROR_RegistrationFailed(className));
	}
	
	// check that name does not contain illegal characters
	if (!IdentifierIsValid(className)) {
		UG_LOG("### Registry ERROR: Trying to register class '" 
				<< className << "' that"
				<< " contains illegal characters.\n"
				<< GetIdentifierMessage()
				<< "\n### Please change register process. Aborting ...\n");
		throw(UG_REGISTRY_ERROR_RegistrationFailed(className));
	}

//	check
	check_base_class<TClass, TBaseClass>(className);

//	new class pointer
	ExportedClass<TClass>* newClass = NULL;

//	try creation of new class
	try { newClass = new ExportedClass<TClass>(className, group, tooltip);}
	catch(ug::bridge::REGISTRY_ERROR_ClassAlreadyNamed ex)
	{
		UG_LOG("### Registry ERROR: Trying to register class with name '"<<className
				<< "', that has already been named. This is not allowed. "
				<< "\n### Please change register process. Aborting ...\n");
		throw(UG_REGISTRY_ERROR_RegistrationFailed(className));
	}
	catch(ug::bridge::REGISTRY_ERROR_Message ex)
	{
		UG_LOG("### Registry ERROR: " << ex.msg << "\n");
		throw(UG_REGISTRY_ERROR_RegistrationFailed("(unknown)"));
	}


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
		UG_LOG("### Registry ERROR: Trying to register class name '"<<className
				<< "', that is already used by another class in this registry."
				<< "\n### Please change register process. Aborting ...\n");
		throw(UG_REGISTRY_ERROR_RegistrationFailed(className));
	}
//	check that className is not already used as a group name
	if(groupname_registered(className))
	{
		UG_LOG("### Registry ERROR: Trying to register class name '"<<className
				<< "', that is already used by another group in this registry."
				<< "\n### Please change register process. Aborting ...\n");
		throw(UG_REGISTRY_ERROR_RegistrationFailed(className));
	}
// 	check that name is not empty
	if(className.empty())
	{
		UG_LOG("### Registry ERROR: Trying to register empty class name."
				<< "\n### Please change register process. Aborting ...\n");
		throw(UG_REGISTRY_ERROR_RegistrationFailed(className));
	}
	
	// check that name does not contain illegal characters
	if (!IdentifierIsValid(className)) {
		UG_LOG("### Registry ERROR: Trying to register class '" 
				<< className << "' that"
				<< " contains illegal characters.\n"
				<< GetIdentifierMessage()
				<< "\n### Please change register process. Aborting ...\n");
		throw(UG_REGISTRY_ERROR_RegistrationFailed(className));
	}

//	check
	check_base_class<TClass, TBaseClass1>(className);
	check_base_class<TClass, TBaseClass2>(className);

//	new class pointer
	ExportedClass<TClass>* newClass = NULL;

//	try creation of new class
	try { newClass = new ExportedClass<TClass>(className, group, tooltip);}
	catch(ug::bridge::REGISTRY_ERROR_ClassAlreadyNamed ex)
	{
		UG_LOG("### Registry ERROR: Trying to register class with name '"<<className
				<< "', that has already been named. This is not allowed. "
				<< "\n### Please change register process. Aborting ...\n");
		throw(UG_REGISTRY_ERROR_RegistrationFailed(className));
	}
	catch(ug::bridge::REGISTRY_ERROR_Message ex)
	{
		UG_LOG("### Registry ERROR: " << ex.msg << "\n");
		throw(UG_REGISTRY_ERROR_RegistrationFailed("(unknown)"));
	}

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
	UG_LOG("### Registry ERROR: Trying to get class with name '" << name
			<< "', that has not yet been registered to this Registry."
			<< "\n### Please change register process. Aborting ...\n");
	throw(UG_REGISTRY_ERROR_RegistrationFailed(name));
}

}//	end of namespace
}//	end of namespace

#endif
