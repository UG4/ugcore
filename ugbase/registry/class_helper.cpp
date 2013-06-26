/**
 * \file class_helper.cpp
 *
 * \author Martin Rupp
 *
 * \date 20.10.2010
 *
 * ClassHierarchy implementation, GetClassHierarchy
 *
 * Goethe-Center for Scientific Computing 2009-2010.
 */

#include <string>
#include <vector>
#include <algorithm>

#include "class_helper.h"
#include "registry.h"
#include "common/util/string_util.h"

using namespace std;

namespace ug
{
namespace bridge
{
extern Registry& GetUGRegistry();


void ClassHierarchy::insert_class(const IExportedClass &c)
{

	//	get name and visualization options of function
	std::vector<std::string> vGroups;
	TokenizeString(c.group(), vGroups, '/');

	ClassHierarchy *base = this;

	for(vector<std::string>::const_iterator it = vGroups.begin(); it != vGroups.end(); ++it)
	{
		const std::string thename = TrimString(*it);
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

void ClassHierarchy::sort()
{
	std::sort(subclasses.begin(), subclasses.end());
	for(size_t i=0; i<subclasses.size(); i++)
		subclasses[i].sort();
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

void GetClassHierarchy(ClassHierarchy &hierarchy, const Registry &reg)
{
	hierarchy.subclasses.clear();
	hierarchy.name = "UGBase";
	for(size_t i=0; i<reg.num_classes(); ++i)
		hierarchy.insert_class(reg.get_class(i));
	hierarchy.sort();
}





void PrintClassSubHierarchy(ClassHierarchy &c, int level)
{
	for(int j=0; j<level; j++) UG_LOG("  ");
	UG_LOG(c.name << endl);
	if(c.subclasses.size())
	{
		for(size_t i=0; i<c.subclasses.size(); i++)
			PrintClassSubHierarchy(c.subclasses[i], level+1);
	}
}

string ClassHierarchyString(const Registry &reg, const char *classname)
{
	std::stringstream ss;
	const IExportedClass *c = reg.get_class(classname);
	if(c == NULL)
	{
		ss << "Class name " << classname << " not found\n";
		return ss.str();
	}

	ss << "\nClass Hierarchy of " << classname << "\n";

	int level = 0;
	const std::vector<const char*> *names = c->class_names();
	if(names != NULL && !names->empty())
	{
		for(int i = names->size()-1; i>0; i--)
		{
			for(int j=0; j<level; j++) ss << "  ";
			ss << names->at(i) << endl;
			level++;
		}
	}

	ClassHierarchy hierarchy;
	GetClassHierarchy(hierarchy, reg);
	ClassHierarchy *ch = hierarchy.find_class(classname);
	if(ch)
		PrintClassSubHierarchy(*ch, level);
	else
	{
		for(int j=0; j<level; j++) ss << "  ";
		ss << classname;
	}

	return ss.str();
}



/**
 *
 * \brief Gets a description of the i-th parameter of a ParameterStack
 * todo: perhaps this function would be better somewhere else like in parameter_stack.cpp
  */
string ParameterToString(const ParameterInfo &par, int i)
{
	string res = string("");
	if(par.is_vector(i)) res.append("std::vector<");
	switch(par.type(i))
	{
	default:
	case Variant::VT_INVALID:
		res.append("unknown");
		break;
	case Variant::VT_BOOL:
		res.append("bool");
		break;
	case Variant::VT_INT:
		res.append("integer");
		break;
	case Variant::VT_FLOAT:
		res.append("number");
		break;
	case Variant::VT_DOUBLE:
		res.append("number");
		break;
	case Variant::VT_CSTRING:
		res.append("string");
		break;
	case Variant::VT_STDSTRING:
		res.append("string");
		break;
	case Variant::VT_POINTER:
		res.append(par.class_name(i)).append("*");
		break;
	case Variant::VT_CONST_POINTER:
		res.append("const ").append(par.class_name(i)).append("*");
		break;
	case Variant::VT_SMART_POINTER:
		res.append("SmartPtr<").append(par.class_name(i)).append(">");
		break;
	case Variant::VT_CONST_SMART_POINTER:
		res.append("ConstSmartPtr<").append(par.class_name(i)).append(">");
		break;
	}
	if(par.is_vector(i)) res.append(">");
	return res;
}

template<typename T>
string PrintParametersIn(const T &thefunc, const char*highlightclassname)
{
	std::stringstream ss;
	ss << "(";
	for(size_t i=0; i < (size_t)thefunc.params_in().size(); ++i)
	{
		if(i>0) ss << ", ";
		bool b=false;
		if(highlightclassname != NULL && thefunc.params_in().class_name(i) != NULL &&
				strcmp(thefunc.params_in().class_name(i), highlightclassname)==0)
			b = true;
		if(b) ss << "[";
		ss << ParameterToString(thefunc.params_in(), i);
		if(i < thefunc.num_parameter())
			ss << " " << thefunc.parameter_name(i);
		if(b) ss << "]";
	}
	ss << ")";
	return ss.str();
}

//bool PrintParametersIn(const ExportedFunctionBase &thefunc, const char*highlightclassname = NULL);
//bool PrintParametersIn(const ExportedConstructor &thefunc, const char*highlightclassname = NULL);

string PrintParametersOut(const ExportedFunctionBase &thefunc)
{
	std::stringstream ss;
	if(thefunc.params_out().size() == 1)
	{
		ss << ParameterToString(thefunc.params_out(), 0);
		//file << " " << thefunc.return_name();
		ss << " ";
	}
	else if(thefunc.params_out().size() > 1)
	{
		ss << "(";
		for(int i=0; i < thefunc.params_out().size(); ++i)
		{
			if(i>0) ss << ", ";
			ss << ParameterToString(thefunc.params_out(), i);

		}
		ss << ") ";
	}
	return ss.str();
}

/**
 *
 * Prints parameters of the function thefunc.
 * If highlightclassname != NULL, it highlights parameters which implement the highlightclassname class.
 */
string FunctionInfo(const ExportedFunctionBase &thefunc, bool isConst, const char *classname, const char *highlightclassname, bool bPrintHelp)
{
	std::stringstream ss;
	PrintParametersOut(thefunc);
	if(classname)
		ss << classname << ":";

	ss << thefunc.name();

	PrintParametersIn<ExportedFunctionBase>(thefunc, highlightclassname);


	if(isConst) ss << " const";
	if(bPrintHelp)
	{
		if(thefunc.tooltip().length() != 0)
			ss << "\n     tooltip: " << thefunc.tooltip();
		if(thefunc.help().length() != 0)
			ss << "\n     help: " << thefunc.help();
	}
	return ss.str();
}

/**
 *
 * Prints parameters of the constructor constr.
 * If highlightclassname != NULL, it highlights parameters which implement the highlightclassname class.
 */
string ConstructorInfo(const ExportedConstructor &constr,
		const char *classname, const char *highlightclassname)
{
	UG_LOG(classname << ":" << classname);
	PrintParametersIn<ExportedConstructor>(constr, highlightclassname);
}


const ExportedFunction *FindFunction(const Registry &reg, const char *functionname)
{
	for(size_t i=0; i<reg.num_functions(); i++)
	{
		if(strcmp(functionname, reg.get_function(i).name().c_str()) == 0)
			return &reg.get_function(i);
	}
	return NULL;
}

/**
 *
 * searches for a function named functionname in the registry and prints it
 */
string FunctionInfo(const Registry &reg, const char *functionname)
{
	const ExportedFunction *f = FindFunction(reg, functionname);
	if(f)
	{
		return FunctionInfo(*f, false, NULL, NULL, true);
	}
	else
		return "";
}

/**
 *
 * \brief Prints the (const) method of one class
 */
string ClassInfo(const IExportedClass &c)
{
	std::stringstream ss;
	ss << "class " << c.name() << "\n";
	if(c.is_instantiable())
		ss << " has constructor\n";
	else
		ss << " has no constructor\n";
	ss << " " << c.num_methods() << " method(s):" << endl;

	for(size_t k=0; k<c.num_methods(); ++k)
	{
		ss << " - ";
		ss << FunctionInfo(c.get_method(k), false, NULL, NULL, true);
		ss << endl;
	}
	ss << " " << c.num_const_methods() << " const method(s):" << endl;
	for(size_t k=0; k<c.num_const_methods(); ++k)
	{
		ss << " - ";
		ss << FunctionInfo(c.get_const_method(k), true);
		ss << endl;
	}
}


/**
 *
 * Searches the classname in the Registry and prints info of the class
 */
string ClassInfo(const Registry &reg, const char *classname)
{
	// search registry for that class
	const IExportedClass *c = reg.get_class(classname);
	if(c)
	{
		return ClassInfo(*c);
	}
	else
		return "";
}


/**
 *
 * \return true, if the class classname is in a parameter in the ParameterStack par
 */
bool IsClassInParameters(const ParameterInfo &par, const char *classname)
{
	int i;
	for(i=0; i<par.size(); ++i)
	{
		if(par.type(i) != Variant::VT_POINTER && par.type(i) != Variant::VT_CONST_POINTER)
			continue;
		if(par.class_name_node(i) != NULL && strcmp(par.class_name(i), classname)==0)
			break;
	}

	if(i==par.size()) return false;
	else return true;
}


/**
 *
 * \param reg			Registry
 * \param classname 	the class (and only this class) to print usage in functions/member functions of.
 * \param OutParameters
 */
string ClassUsageExact(const Registry &reg, const char *classname, bool OutParameters)
{
	std::stringstream ss;
	// check functions
	for(size_t i=0; i<reg.num_functions(); i++)
	{
		const ExportedFunctionBase &thefunc = reg.get_function(i);
		if((!OutParameters && IsClassInParameters(thefunc.params_in(), classname)) ||
				(OutParameters && IsClassInParameters(thefunc.params_out(), classname)))
		{
			ss << " " << FunctionInfo(thefunc, false, classname) << endl;
		}
	}

	// check classes
	for(size_t i=0; i<reg.num_classes(); i++)
	{
		const IExportedClass &c = reg.get_class(i);
		for(size_t i=0; i<c.num_methods(); i++)
		{
			const ExportedFunctionBase &thefunc = c.get_method(i);
			if((!OutParameters && IsClassInParameters(thefunc.params_in(), classname)) ||
					(OutParameters && IsClassInParameters(thefunc.params_out(), classname)))
			{
				ss << " " << FunctionInfo(thefunc, false, c.name().c_str(), classname)  << endl;
			}
		}

		for(size_t i=0; i<c.num_const_methods(); i++)
		{
			const ExportedFunctionBase &thefunc = c.get_const_method(i);
			if((!OutParameters && IsClassInParameters(thefunc.params_in(), classname)) ||
					(OutParameters && IsClassInParameters(thefunc.params_out(), classname)))
			{
				ss << " " << FunctionInfo(thefunc, false, c.name().c_str(), classname) << endl;
			}
		}
	}
	return ss.str();
}


} // namespace bridge
} // namespace ug

