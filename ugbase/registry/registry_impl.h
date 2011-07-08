// created by Andreas Vogel, Sebastian Reiter
// s.b.reiter@googlemail.com
// 08.07.2011 (m,d,y)

#ifndef __H__UG_BRIDGE__REGISTRY_IMPL__
#define __H__UG_BRIDGE__REGISTRY_IMPL__

namespace ug{
namespace bridge
{
//////////////////////
// global functions
//////////////////////

template<class TFunc>
Registry& Registry::
add_function(const char* funcName, TFunc func, const char* group,
			 const char* retValInfos, const char* paramInfos,
			 const char* tooltip, const char* help)
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
	{
		const size_t start = strippedMethodName.find_first_not_of(" \t");
		const size_t end = strippedMethodName.find_last_not_of(" \t");
		if(start != std::string::npos && end != std::string::npos)
			strippedMethodName = strippedMethodName.substr(start, end - start + 1);
	}
	{
		const size_t start = methodOptions.find_first_not_of(" \t");
		const size_t end = methodOptions.find_last_not_of(" \t");
		if(start != std::string::npos && end != std::string::npos)
			methodOptions = methodOptions.substr(start, end - start + 1);
	}

// 	check that name is not empty
	if(strippedMethodName.empty())
	{
		std::cout << "### Registry ERROR: Trying to register empty function name."
				<< "\n### Please change register process. Aborting ..." << std::endl;
		throw(UG_REGISTRY_ERROR_RegistrationFailed(strippedMethodName.c_str()));
	}

//	if the function is already in use, we have to add an overload
	ExportedFunctionGroup* funcGrp = get_exported_function_group(strippedMethodName.c_str());
	if(!funcGrp)
	{
	//	we have to create a new function group
		funcGrp = new ExportedFunctionGroup(strippedMethodName.c_str());
		m_vFunction.push_back(funcGrp);
	}

//  add an overload to the function group
	bool success = funcGrp->add_overload(func, &FunctionProxy<TFunc>::apply,
										methodOptions.c_str(), group,
										retValInfos, paramInfos,
										tooltip, help);

	if(!success){
		std::cout << "### Registry ERROR: Trying to register function name '" << funcName
				<< "', that is already used by another function in this registry."
				<< "\n### Please change register process. Aborting ..." << std::endl;
		throw(UG_REGISTRY_ERROR_RegistrationFailed(strippedMethodName.c_str()));
	}

	return *this;
}


///////////////////
// classes
///////////////////

template <typename TClass, typename TBaseClass>
void Registry::
check_base_class(const char* className)
{
//	check that className is not already used
	if(classname_registered(className))
	{
		std::cout << "### Registry ERROR: Trying to register class name '" << className
				<< "', that is already used by another class in this registry."
				<< "\n### Please change register process. Aborting ..." << std::endl;
		throw(UG_REGISTRY_ERROR_RegistrationFailed(className));
	}
// 	check that name is not empty
	if(strlen(className) == 0)
	{
		std::cout << "### Registry ERROR: Trying to register empty class name."
				<< "\n### Please change register process. Aborting ..." << std::endl;
		throw(UG_REGISTRY_ERROR_RegistrationFailed(className));
	}

// 	check that base class is not same type as class
	if(typeid(TClass) == typeid(TBaseClass))
	{
		std::cout << "### Registry ERROR: Trying to register class " << className
				<< "\n### that derives from itself. Aborting ..." << std::endl;
		throw(UG_REGISTRY_ERROR_RegistrationFailed(className));
	}

// 	check that class derives from base class
	if(boost::is_base_of<TBaseClass, TClass>::value == false)
	{
		std::cout << "### Registry ERROR: Trying to register class " << className
				<< "\n### with base class that is no base class. Aborting ..." << std::endl;
		throw(UG_REGISTRY_ERROR_RegistrationFailed(className));
	}
}

template <typename TClass>
ExportedClass_<TClass>& Registry::
add_class_(const char* className, const char* group, const char *tooltip)
{
//	check that className is not already used
	if(classname_registered(className))
	{
		std::cout << "### Registry ERROR: Trying to register class name '" << className
				<< "', that is already used by another class in this registry."
				<< "\n### Please change register process. Aborting ..." << std::endl;
		throw(UG_REGISTRY_ERROR_RegistrationFailed(className));
	}
// 	check that name is not empty
	if(strlen(className) == 0)
	{
		std::cout << "### Registry ERROR: Trying to register empty class name."
				<< "\n### Please change register process. Aborting ..." << std::endl;
		throw(UG_REGISTRY_ERROR_RegistrationFailed(className));
	}

//	new class pointer
	ExportedClass_<TClass>* newClass = NULL;

//	try creation
	try
	{
		newClass = new ExportedClass_<TClass>(className, group, tooltip);
	}
	catch(ug::bridge::UG_REGISTRY_ERROR_ClassAlreadyNamed ex)
	{
		std::cout << "### Registry ERROR: Trying to register class with name '" << className
				<< "', that has already been named. This is not allowed. "
				<< "\n### Please change register process. Aborting ..." << std::endl;
		throw(UG_REGISTRY_ERROR_RegistrationFailed(className));
	}

//	add new class to list of classes
	m_vClass.push_back(newClass);

	return *newClass;
}

template <typename TClass, typename TBaseClass>
ExportedClass_<TClass>& Registry::
add_class_(const char* className, const char* group, const char *tooltip)
{
//	check
	check_base_class<TClass, TBaseClass>(className);

//	new class pointer
	ExportedClass_<TClass>* newClass = NULL;

//	try creation of new class
	try { newClass = new ExportedClass_<TClass>(className, group, tooltip);}
	catch(ug::bridge::UG_REGISTRY_ERROR_ClassAlreadyNamed ex)
	{
		std::cout << "### Registry ERROR: Trying to register class with name '" << className
				<< "', that has already been named. This is not allowed. "
				<< "\n### Please change register process. Aborting ..." << std::endl;
		throw(UG_REGISTRY_ERROR_RegistrationFailed(className));
	}

// 	set base class names
	try
	{
		ClassNameProvider<TClass>::template set_name<TBaseClass>(className, group);
	}
	catch(ug::bridge::UG_REGISTRY_ERROR_ClassUnknownToRegistry ex)
	{
		std::cout <<"### Registry ERROR: Trying to register class with name '" << className
				<< "', that derives from class, that has not yet been registered to this Registry."
				<< "\n### Please change register process. Aborting ..." << std::endl;
		throw(UG_REGISTRY_ERROR_RegistrationFailed(className));
	}

//	add cast function
	ClassCastProvider::add_cast_func<TBaseClass, TClass>();

//	add new class to list of classes
	m_vClass.push_back(newClass);
	return *newClass;
}

template <typename TClass, typename TBaseClass1, typename TBaseClass2>
ExportedClass_<TClass>& Registry::
add_class_(const char* className, const char* group)
{
//	check
	check_base_class<TClass, TBaseClass1>(className);
	check_base_class<TClass, TBaseClass2>(className);

//	new class pointer
	ExportedClass_<TClass>* newClass = NULL;

//	try creation of new class
	try { newClass = new ExportedClass_<TClass>(className, group);}
	catch(ug::bridge::UG_REGISTRY_ERROR_ClassAlreadyNamed ex)
	{
		std::cout << "### Registry ERROR: Trying to register class with name '" << className
				<< "', that has already been named. This is not allowed. "
				<< "\n### Please change register process. Aborting ..." << std::endl;
		throw(UG_REGISTRY_ERROR_RegistrationFailed(className));
	}

// 	set base class names
	try
	{
		ClassNameProvider<TClass>::template set_name<TBaseClass1, TBaseClass2>(className, group);
	}
	catch(ug::bridge::UG_REGISTRY_ERROR_ClassUnknownToRegistry ex)
	{
		std::cout <<"### Registry ERROR: Trying to register class with name '" << className
				<< "', that derives from class, that has not yet been registered to this Registry."
				<< "\n### Please change register process. Aborting ..." << std::endl;
		throw(UG_REGISTRY_ERROR_RegistrationFailed(className));
	}

//	add cast function
	ClassCastProvider::add_cast_func<TBaseClass1, TClass>();
	ClassCastProvider::add_cast_func<TBaseClass2, TClass>();

//	add new class to list of classes
	m_vClass.push_back(newClass);
	return *newClass;
}

template <typename TClass>
ExportedClass_<TClass>& Registry::
get_class_()
{
	const char* name;
// get class names
	try
	{
		name = ClassNameProvider<TClass>::name();
	}
	catch(ug::bridge::UG_REGISTRY_ERROR_ClassUnknownToRegistry ex)
	{
		std::cout <<"### Registry ERROR: Trying to get class "
				<< "that has not yet been named."
				<< "\n### Please change register process. Aborting ..." << std::endl;
		throw(UG_REGISTRY_ERROR_RegistrationFailed(""));
	}

//	look for class in this registry
	for(size_t i = 0; i < m_vClass.size(); ++i)
	{
	//  compare strings
		if(strcmp(name, m_vClass[i]->name()) == 0)
			return *dynamic_cast<ExportedClass_<TClass>* >(m_vClass[i]);
	}

//	not found
	std::cout <<"### Registry ERROR: Trying to get class with name '" << name
			<< "', that has not yet been registered to this Registry."
			<< "\n### Please change register process. Aborting ..." << std::endl;
	throw(UG_REGISTRY_ERROR_RegistrationFailed(name));
}

}//	end of namespace
}//	end of namespace

#endif
