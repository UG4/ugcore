/**
 * \file class_helper.cpp
 *
 * \author Martin Rupp
 *
 * \date 20.10.2010
 *
 * ClassHierarchy implementation, GetClassHierarchy, FindClass
 *
 * Goethe-Center for Scientific Computing 2009-2010.
 */

#include <string>
#include <vector>
#include <sstream>
#include "parameter_stack.h"
#include "function_traits.h"
#include "param_to_type_value_list.h"
#include <iostream>

#include "class_helper.h"
#include "registry.h"
//#include "ug_script/ug_script.h"

namespace ug
{
namespace bridge
{
extern Registry& GetUGRegistry();

static std::string trim(const std::string& str)
{
	const size_t start = str.find_first_not_of(" \t");
	const size_t end = str.find_last_not_of(" \t");
	if(start == std::string::npos || end == std::string::npos) return "";
	return str.substr(start, end - start + 1);
}

static void tokenize(const std::string& str, std::vector<std::string>& tokens, const char delimiter)
{
	tokens.clear();
	std::stringstream tokenstream;
	tokenstream << str;
	std::string token;

	while ( std::getline (tokenstream, token, delimiter ) )
	{
		tokens.push_back(trim(token));
	}
}

void ClassHierarchy::insert_class(const IExportedClass &c)
{
	//	get name and visualization options of function
	std::vector<std::string> vGroups;
	tokenize(c.group(), vGroups, '/');

	ClassHierarchy *base = this;
	for(vector<std::string>::const_iterator it = vGroups.begin(); it != vGroups.end(); ++it)
	{
		const std::string &thename = (*it);
		if(thename.length() <= 0) continue;
		bool bFound = false;
		for(size_t j=0; j<base->subclasses.size(); j++)
		{
			if(base->subclasses.at(j).name.compare(thename) == 0)
			{
				base = &base->subclasses.at(j);
				bFound = true;
				break;
			}
		}

		if(!bFound)
		{
			ClassHierarchy newclass;
			newclass.name = thename;
			newclass.bGroup = true;
			base->subclasses.push_back(newclass);
			base = &base->subclasses.back();
		}
	}

	const vector<const char *> *pNames = c.class_names();

	if(pNames == NULL) return;

	for(vector<const char*>::const_reverse_iterator rit = pNames->rbegin(); rit != pNames->rend(); ++rit)
	{
		const char *thename = (*rit);
		bool bFound = false;
		for(size_t j=0; j<base->subclasses.size(); j++)
		{
			if(base->subclasses.at(j).name.compare(thename) == 0)
			{
				base = &base->subclasses.at(j);
				bFound = true;
				break;
			}
		}

		if(!bFound)
		{
			ClassHierarchy newclass;
			newclass.name = thename;
			newclass.bGroup = false;
			base->subclasses.push_back(newclass);
			base = &base->subclasses.back();
		}
	}
}


ClassHierarchy *ClassHierarchy::find_class(const char *classname)
{
	if(name.compare(classname) == 0)
		return this;
	for(size_t i=0; i < subclasses.size(); i++)
	{
		ClassHierarchy *c = subclasses[i].find_class(classname);
		if(c) return c;
	}
	return NULL;
}

void GetClassHierarchy(ClassHierarchy &hierarchy, const bridge::Registry &reg)
{
	hierarchy.subclasses.clear();
	hierarchy.name = "UGBase";
	for(size_t i=0; i<reg.num_classes(); ++i)
		hierarchy.insert_class(reg.get_class(i));
}

const IExportedClass *FindClass(const char* classname)
{
	bridge::Registry &reg = ug::bridge::GetUGRegistry();
	for(size_t j=0; j<reg.num_classes(); ++j)
		if(strcmp(classname, reg.get_class(j).name()) == 0)
		{
			return &reg.get_class(j);
			break;
		}
	return NULL;
}

} // namespace bridge
} // namespace ug

