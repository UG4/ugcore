// created by Andreas Vogel, Sebastian Reiter
// s.b.reiter@googlemail.com
// 08.07.2011 (m,d,y)
 
#include <string>

#include "registry.h"

namespace ug{
namespace bridge
{

Registry::Registry()
{
//	register native types as provided in ParameterStack
	add_class_<bool>("bool");
	add_class_<int>("int");
	add_class_<size_t>("size_t");
	add_class_<float>("float");
	add_class_<double>("double");
	add_class_<char>("char");
	add_class_<std::string>("string");
}

Registry::Registry(const Registry& reg)
{
}

Registry::~Registry()
{
//  delete registered functions
	for(size_t i = 0; i < m_vFunction.size(); ++i)
	{
		if(m_vFunction[i] != NULL)
			delete m_vFunction[i];
	}
//  delete registered classes
	for(size_t i = 0; i < m_vClass.size(); ++i)
	{
		if(m_vClass[i] != NULL)
			delete m_vClass[i];
	}
//	delete registered class groups
	for(size_t i = 0; i < m_vClassGroups.size(); ++i){
		delete m_vClassGroups[i];
	}
}


////////////////////////
//	callbacks
////////////////////////

void Registry::add_callback(FuncRegistryChanged callback)
{
	m_callbacksRegChanged.push_back(callback);
}

bool Registry::registry_changed()
{
//	check that the registered components are correct
	if(!check_consistency()) return false;

//	iterate through all callbacks and call them
	for(size_t i = 0; i < m_callbacksRegChanged.size(); ++i){
		m_callbacksRegChanged[i](this);
	}

//	ok
	return true;
}

//////////////////////
// global functions
//////////////////////

size_t Registry::num_functions() const
{
	return m_vFunction.size();
}

ExportedFunction& Registry::get_function(size_t ind)
{
	return *m_vFunction.at(ind)->get_overload(0);
}

size_t Registry::num_overloads(size_t ind)
{
	return m_vFunction.at(ind)->num_overloads();
}

ExportedFunction& Registry::get_overload(size_t funcInd, size_t oInd)
{
	return *m_vFunction.at(funcInd)->get_overload(oInd);
}

ExportedFunctionGroup& Registry::get_function_group(size_t ind)
{
	return *m_vFunction.at(ind);
}


///////////////////
// classes
///////////////////

size_t Registry::num_classes() const
{
	return m_vClass.size();
}

/// returns an exported function
const IExportedClass& Registry::get_class(size_t ind) const
{
	return *m_vClass.at(ind);
}

IExportedClass* Registry::get_class(const char* name)
{
//todo:	use a map to access classes by name.
	for(size_t i = 0; i < m_vClass.size(); ++i)
	{
	//  compare strings
		if(strcmp(name, m_vClass[i]->name()) == 0)
			return m_vClass[i];
	}
	return NULL;
}

bool Registry::check_consistency()
{
	size_t globFctUndef = 0;
//	check global functions
	for(size_t i=0; i<num_functions(); i++)
	{
	//	get function
		ExportedFunctionGroup& funcGrp = get_function_group(i);

	//	check all overloads
		for(size_t j = 0; j < funcGrp.num_overloads(); ++j){
			if(!funcGrp.get_overload(j)->check_consistency())
				globFctUndef++;
		}
	}

// 	check classes and their methods
	size_t baseClassUndef = 0;
	size_t methodUndef = 0;
	for(size_t i=0; i<num_classes(); i++)
	{
	//	get class
		const bridge::IExportedClass &c = get_class(i);

	//	check class (e.g. that base classes have been named)
		if(!c.check_consistency())
			baseClassUndef++;

	//	check methods
		for(size_t j=0; j<c.num_methods(); j++)
			if(!c.get_method(j).check_consistency(c.name()))
				methodUndef++;

		for(size_t j=0; j<c.num_const_methods(); j++)
			if(!c.get_const_method(j).check_consistency(c.name()))
				methodUndef++;
	}

//	log error messages
	if(globFctUndef > 0)
		UG_LOG("#### ERROR in 'Registry::check_consistency': "<<globFctUndef<<
		       " global Functions are using undeclared Classes.\n");
	if(methodUndef > 0)
		UG_LOG("#### ERROR in 'Registry::check_consistency': "<<methodUndef<<
		       " Methods are using undeclared Classes.\n");
	if(baseClassUndef > 0)
		UG_LOG("#### ERROR in 'Registry::check_consistency': "<<baseClassUndef<<
		       " Base-Classes are using undeclared Classes.\n");

//	return false if undefined classes have been found
	if(globFctUndef > 0) return false;
	if(methodUndef > 0) return false;
	if(baseClassUndef > 0) return false;

//	everything fine
	return true;
}


size_t Registry::num_class_groups() const
{
	return m_vClassGroups.size();
}

const ClassGroupDesc* Registry::get_class_group(size_t i) const
{
	return m_vClassGroups[i];
}

ClassGroupDesc* Registry::get_class_group(size_t i)
{
	return m_vClassGroups[i];
}

ClassGroupDesc* Registry::get_class_group(const char* name)
{
//todo:	use a map to quickly access classGroups by name
	for(size_t i = 0; i < m_vClassGroups.size(); ++i){
		if(strcmp(m_vClassGroups[i]->name(), name) == 0)
			return m_vClassGroups[i];
	}

//	since we reached this point, no class-group with the given name exists.
	ClassGroupDesc* classGroup = new ClassGroupDesc();
	classGroup->set_name(name);
	m_vClassGroups.push_back(classGroup);

	return classGroup;
}

const ClassGroupDesc* Registry::get_class_group(const char* name) const
{
//todo:	use a map to quickly access classGroups by name
	for(size_t i = 0; i < m_vClassGroups.size(); ++i){
		if(strcmp(m_vClassGroups[i]->name(), name) == 0)
			return m_vClassGroups[i];
	}
//	since we reached this point, no class-group with the given name exists.
	return NULL;
}

void Registry::add_class_to_group(const char* className, const char* groupName)
{
//todo: make sure that no class with groupName exists.
	ClassGroupDesc* groupDesc = get_class_group(groupName);
//todo:	make sure that groupDesc does not already contain className.
	IExportedClass* expClass = get_class(className);
	UG_ASSERT(expClass, "The given class has to be registered before "
						"adding it to a group: " << className);
	if(expClass)
		groupDesc->add_class(expClass);
}


bool Registry::classname_registered(const char* name)
{
	return get_class(name) != NULL;
}

// returns true if functionname is already used by a function in this registry
bool Registry::functionname_registered(const char* name)
{
	for(size_t i = 0; i < m_vFunction.size(); ++i)
	{
	//  compare strings
		if(strcmp(name, (m_vFunction[i]->name()).c_str()) == 0)
			return true;
	}
	return false;
}

ExportedFunctionGroup* Registry::get_exported_function_group(const char* name)
{
	for(size_t i = 0; i < m_vFunction.size(); ++i)
	{
	//  compare strings
		if(strcmp(name, (m_vFunction[i]->name()).c_str()) == 0)
			return m_vFunction[i];
	}
	return NULL;
}

}// end of namespace
}// end of namespace
