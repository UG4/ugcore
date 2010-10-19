/**
 * \file typeinfo.cpp
 *
 * \author Martin Rupp
 *
 * \date 15.10.2010
 *
 * Goethe-Center for Scientific Computing 2010.
 *
 * Comfort functions for the lua ug shell.
 * example: TypeInfo("
 */

#include <iostream>
#include <sstream>

#include <dirent.h>
#include <sys/stat.h>

#include "ug.h"

#include "ug_script/ug_script.h"
extern "C"
{
#include "ug_script/externals/lua/lstate.h"
}

namespace ug
{
namespace bridge
{


const IExportedClass *FindClass(const char* classname)
{
	bridge::Registry &reg = GetUGRegistry();
	for(size_t j=0; j<reg.num_classes(); ++j)
		if(strcmp(classname, reg.get_class(j).name()) == 0)
		{
			return &reg.get_class(j);
			break;
		}
	return NULL;
}

/**
 *
 * \brief Gets a description of the i-th parameter of a ParameterStack
 * todo: perhaps this function would be better somewhere else like in parameter_stack.cpp
  */
string ParameterToString(const bridge::ParameterStack &par, int i)
{
	switch(par.get_type(i))
	{
	default:
	case PT_UNKNOWN:
		return string("unknown");
	case PT_BOOL:
		return string("bool");

	case PT_INTEGER:
		return string("integer");

	case PT_NUMBER:
		return string("number");

	case PT_STRING:
		return string("string");

	case PT_POINTER:
		return string(par.class_name(i)).append("*");
		break;

	case PT_CONST_POINTER:
		return string("const ").append(par.class_name(i)).append("*");
		break;
	}
}

/**
 *
 * Prints parameters of the ParameterStack par.
 * If highlightclassname != NULL, it highlights parameters which implement the highlightclassname class.
 */
bool PrintParameters(const bridge::ParameterStack &par, const char *highlightclassname=NULL)
{
	for(int i=0; i<par.size(); ++i)
	{
		if(i>0) UG_LOG(", ");
		bool b=false;
		if(highlightclassname != NULL && par.class_names(i) != NULL && strcmp(par.class_name(i), highlightclassname)==0)
			b = true;
		if(b) UG_LOG("[");
		UG_LOG(ParameterToString(par, i));;
		if(b) UG_LOG("]");
	}
	return true;
}

/**
 *
 * Prints parameters of the function thefunc.
 * If highlightclassname != NULL, it highlights parameters which implement the highlightclassname class.
 */
void PrintFunctionInfo(const bridge::ExportedFunctionBase &thefunc, bool isConst, const char *highlightclassname=NULL)
{
	UG_LOG(" " << thefunc.name() << " (");

	PrintParameters(thefunc.params_in(), highlightclassname);
	UG_LOG(")");
	if(isConst) { UG_LOG(" const"); }

	if(thefunc.params_out().size() > 0)
	{
		UG_LOG(" => (");
		PrintParameters(thefunc.params_out(), highlightclassname);
		UG_LOG(")\n");
	}
	else
		UG_LOG("\n");
}

/**
 *
 * searches for a function named functionname in the registry and prints it
 */
bool PrintFunctionInfo(const char *functionname)
{
	bridge::Registry &reg = GetUGRegistry();
	for(size_t i=0; i<reg.num_functions(); i++)
	{
		if(strcmp(functionname, reg.get_function(i).name().c_str()) == 0)
		{
			PrintFunctionInfo(reg.get_function(i), false);
			return true;
		}
	}
	return false;
}

/**
 *
 * \brief Prints the (const) method of one class
 */
void PrintClassInfo(const IExportedClass &c)
{
	UG_LOG("class " << c.name() << ", " << c.num_methods() << " method(s), " <<
		c.num_const_methods() << " const method(s):" << endl);
	for(size_t k=0; k<c.num_methods(); ++k)
		PrintFunctionInfo(c.get_method(k), false);
	for(size_t k=0; k<c.num_const_methods(); ++k)
		PrintFunctionInfo(c.get_const_method(k), true);
}


/**
 *
 * Searches the classname in the Registry and prints info of the class
 */
bool PrintClassInfo(const char *classname)
{
	// search registry for that class
	const IExportedClass *c = FindClass(classname);
	if(c)
	{
		PrintClassInfo(*c);
		return true;
	}
	else
		return false;
}

namespace lua
{
	string GetLuaTypeString(lua_State* L, int index);
}

/**
 * \brief Prints info to a lua type
 * \param 	p			the name of the object in lua.
 * you can use class names, function names or the names of an object
 * - TypeInfo("class") prints all member functions+parameters of this class and its parents
 * - TypeInfo("Function") prints all member functions+parameters
 * - TypeInfo("variable") prints class information if variable is a object of a ug class, otherwise what type in lua it is
 */
bool UGTypeInfo(const char *p)
{
	UG_LOG("\n");
	const IExportedClass *c = FindClass(p);
	if(c)
	{
		const std::vector<const char*> *names = c->class_names();
		for(size_t i=0; i < names->size(); ++i)
			PrintClassInfo(names->at(i));
		UG_LOG(endl);
		return true;
	}

	lua_State* L = script::GetDefaultLuaState();

	struct UserDataWrapper
	{
		bool	is_const;
		void*	obj;
	};

	lua_getglobal(L, p);
	if(lua_isnil(L, -1))
	{
		lua_pop(L, 1);
		UG_LOG(p << " is neither a global variable nor a class name." << endl);
		return false;
	}

	if(lua_isfunction(L, -1))
	{
		UG_LOG(p << " is a function\n");
		PrintFunctionInfo(p);
	}
	else if(lua_iscfunction(L, -1))
	{
		UG_LOG(p << " is a cfunction\n");
		PrintFunctionInfo(p);
	}
	else if(lua_isuserdata(L, -1))
	{
		if(lua_getmetatable(L, -1) == 0)
		{
			UG_LOG("global variable " << p << " has light user data, but no metatable." << endl);
			lua_pop(L, 1); // pop globals
			return false;
		}

		// get names
		lua_pushstring(L, "names");
		lua_rawget(L, -2);
		if(lua_isnil(L, -1) || !lua_isuserdata(L, -1))
		{
			UG_LOG("global variable " << p << " has metatable, but cannot access names." << endl);
			lua_pop(L, 3); // pop metatable, userdata, globals
			return false;
		}
		const std::vector<const char*> *names =
				(const std::vector<const char*>*) lua_touserdata(L, -1);
		lua_pop(L, 2); // pop metatable, userdata
		UG_LOG("Typeinfo for " << p << ": " << endl);
		for(size_t i=0; i < names->size(); ++i)
			PrintClassInfo(names->at(i));
	}
	else
		UG_LOG(p << " is a " << lua::GetLuaTypeString(L, -1) << endl);

	lua_pop(L, 1); // pop globals

	UG_LOG("\n");
	return true;
}


/**
 *
 * \return true, if the class classname is in a parameter in the ParameterStack par
 */
bool IsClassInParameters(const bridge::ParameterStack &par, const char *classname)
{
	int i;
	for(i=0; i<par.size(); ++i)
	{
		if(par.get_type(i) != PT_POINTER && par.get_type(i) != PT_CONST_POINTER)
			continue;
		if(par.class_names(i) != NULL && strcmp(par.class_name(i), classname)==0)
			break;
	}

	if(i==par.size()) return false;
	else return true;
}

/**
 *
 * \param classname the class (and only this class) to print usage in functions/member functions of.
 */
bool ClassUsageExact(const char *classname, bool OutParameters)
{
	bridge::Registry &reg = GetUGRegistry();
	// check functions
	for(size_t i=0; i<reg.num_functions(); i++)
	{
		const bridge::ExportedFunctionBase &thefunc = reg.get_function(i);
		if(!OutParameters && IsClassInParameters(thefunc.params_in(), classname) ||
				OutParameters && IsClassInParameters(thefunc.params_out(), classname))
			PrintFunctionInfo(thefunc, false, classname);
	}

	// check classes
	for(size_t i=0; i<reg.num_classes(); i++)
	{
		const IExportedClass &c = reg.get_class(i);
		for(size_t i=0; i<c.num_methods(); i++)
		{
			const bridge::ExportedFunctionBase &thefunc = c.get_method(i);
			if(!OutParameters && IsClassInParameters(thefunc.params_in(), classname) ||
					OutParameters && IsClassInParameters(thefunc.params_out(), classname))
			{
				UG_LOG(c.name() << " ::");
				PrintFunctionInfo(thefunc, false, classname);
			}
		}

		for(size_t i=0; i<c.num_const_methods(); i++)
		{
			const bridge::ExportedFunctionBase &thefunc = c.get_const_method(i);
			if(!OutParameters && IsClassInParameters(thefunc.params_in(), classname) ||
					OutParameters && IsClassInParameters(thefunc.params_out(), classname))
			{
				UG_LOG(c.name() << " ::");
				PrintFunctionInfo(thefunc, false, classname);
			}
		}
	}
	return true;
}

bool SubClasses(const char *classname)
{
	bridge::Registry &reg = GetUGRegistry();
	const IExportedClass *c = FindClass(classname);
	if(c == NULL)
	{
		UG_LOG("Class name " << classname << " not found\n");
		return false;
	}

	UG_LOG("\n");
	for(size_t j=0; j<reg.num_classes(); ++j)
	{
		const std::vector<const char*> *names = reg.get_class(j).class_names();
		if(names == NULL) continue;

		for(size_t k=0; k<names->size(); k++)
		{
			if(strcmp(names->at(k), classname) == 0)
			{
				UG_LOG(reg.get_class(j).name() << " (");
				for(size_t l=0; l<names->size(); l++)
				{
					if(l>0) UG_LOG(" :: ");
					UG_LOG(names->at(l));
				}
				UG_LOG(")\n");
				break;
			}
		}
	}
	UG_LOG("\n");
	return true;
}

/**
 *
 * \param classname the class to print usage in functions/member functions of (and all its subclasses) .
 * class in in/out parameters is highlighted with [class].
 */
bool ClassUsage(const char *classname)
{
	UG_LOG("\n");

	// find class
	const IExportedClass *c = FindClass(classname);
	if(c == NULL)
	{
		UG_LOG("Class name " << classname << " not found\n");
		return false;
	}

	// print usages in functions

	UG_LOG("--- Functions returning " << classname << ": ---\n");
	ClassUsageExact(classname, true);

	const std::vector<const char*> *names = c->class_names();
	if(names != NULL && names->size() > 1)
	{
		for(size_t i = 1; i<names->size(); i++)
		{
			UG_LOG("--- Functions using " << names->at(i) << ": ---\n");
			ClassUsageExact(names->at(i), false);
		}
	}

	UG_LOG("\n");
	return true;
}
bool ClassNameVecContains(const std::vector<const char*>& names, const char* name);
bool ClassInstantiations(const char *classname)
{
	UG_LOG(endl);
	const IExportedClass *c = FindClass(classname);
	if(c == NULL)
	{
		UG_LOG("Class " << classname << " not found\n");
		return false;
	}

	lua_State* L = script::GetDefaultLuaState();
	bool bFound = false;

	// iterate through all of lua's global string table
	for(int i=0; i<G(L)->strt.size; i++)
	{
		GCObject *obj;
		for (obj = G(L)->strt.hash[i]; obj != NULL; obj = obj->gch.next)
		{
			// get the string
			TString *ts = rawgco2ts(obj);
			if(ts == NULL) continue;

			const char *luastr = getstr(ts);
			// check is of a global variable
			lua_getglobal(L, luastr);
			if(lua_isnil(L, -1) || !lua_isuserdata(L, -1) || lua_getmetatable(L, -1) == 0)
			{
				lua_pop(L, 1); // remove global from stack
				continue; // nope
			}

			// get names
			lua_pushstring(L, "names");
			lua_rawget(L, -2);
			if(lua_isnil(L, -1) || !lua_isuserdata(L, -1))
			{
				lua_pop(L, 3); // pop userdata, metatable, globals
				continue;
			}
			const std::vector<const char*> *names =
					(const std::vector<const char*>*) lua_touserdata(L, -1);
			lua_pop(L, 3); // pop userdata, metatable, globals
			if(ClassNameVecContains(*names, classname))
			{
				bFound = true;
				UG_LOG(setw(10) << left << luastr);
				UG_LOG(" (");
				for(size_t i=0; i<names->size(); i++)
				{
					if(i>0) UG_LOG(" :: ");
					UG_LOG(names->at(i));
				}
				UG_LOG(")\n");
			}
		}
	}
	if(!bFound) UG_LOG("No instantiations of " << classname << " or subclasses found.");
	UG_LOG(endl);
	return true;
}

void RegisterInfoCommands(Registry &reg)
{
	reg.add_function("TypeInfo", &UGTypeInfo);
	reg.add_function("ClassUsage", &ClassUsage);
	reg.add_function("ClassInstantiations" ,&ClassInstantiations);
	reg.add_function("SubClasses" ,&SubClasses);
}


} // namespace bridge
} // namespace ug
