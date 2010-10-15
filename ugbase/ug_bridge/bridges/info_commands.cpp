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
 * and its subclasses
 */
bool PrintClassInfo(const char *classname)
{
	bridge::Registry &reg = GetUGRegistry();

	// search registry for that class
	for(size_t j=0; j<reg.num_classes(); ++j)
	{
		if(strcmp(classname, reg.get_class(j).name()) == 0)
		{
			PrintClassInfo(reg.get_class(j));
			return true;
		}
	}
	return false;
}


/**
 * \brief Prints info to a lua type
 * \param 	p			the name of the object in lua.
 * you can use class names, function names or the names of an object
 */
bool UGTypeInfo(const char *p)
{
	UG_LOG("\n");
	if(PrintClassInfo(p))
		return true;

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
		if(PrintFunctionInfo(p) == false)
			UG_LOG(p << " is a Lua function\n");
	}
	else if(lua_iscfunction(L, -1))
	{
		if(PrintFunctionInfo(p) == false)
			UG_LOG(p << " is a cfunction\n");
	}
	else if(lua_islightuserdata(L, -1))
	{
		UG_LOG(p << " is userdata\n");
	}
	else if(lua_isstring(L, -1))
	{ UG_LOG(p << " is a string: \"" << lua_tostring(L, -1) << "\"\n"); }
	else if(lua_isuserdata(L, -1))
	{
		if(lua_getmetatable(L, -1) == 0)
		{
			UG_LOG("global variable " << p << " has light user data, but no metatable." << endl);
			return false;
		}

		// get names
		lua_pushstring(L, "names");
		lua_rawget(L, -2);
		if(lua_isnil(L, -1) || !lua_isuserdata(L, -1))
		{
			UG_LOG("global variable " << p << " has metatable, but cannot access names." << endl);
			lua_pop(L, 1);
			return false;
		}
		const std::vector<const char*> *names =
				(const std::vector<const char*>*) lua_touserdata(L, -1);
		lua_pop(L, 2);
		UG_LOG("Typeinfo for " << p << ": " << endl);
		for(size_t i=0; i < names->size(); ++i)
			PrintClassInfo((*names)[i]);
	}
	else
	{
		UG_LOG(lua_tostring(L, -1));
	}
	lua_pop(L, 1);

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
bool ClassUsageExact(const char *classname)
{
	bridge::Registry &reg = GetUGRegistry();
	// check functions
	for(size_t i=0; i<reg.num_functions(); i++)
	{
		const bridge::ExportedFunctionBase &thefunc = reg.get_function(i);
		if(IsClassInParameters(thefunc.params_in(), classname) || IsClassInParameters(thefunc.params_out(), classname))
			PrintFunctionInfo(thefunc, false, classname);
	}

	// check classes
	for(size_t i=0; i<reg.num_classes(); i++)
	{
		const IExportedClass &c = reg.get_class(i);
		for(size_t i=0; i<c.num_methods(); i++)
		{
			const bridge::ExportedFunctionBase &thefunc = c.get_method(i);
			if(IsClassInParameters(thefunc.params_in(), classname) || IsClassInParameters(thefunc.params_out(), classname))
			{
				UG_LOG(c.name() << " ::");
				PrintFunctionInfo(thefunc, false, classname);
			}
		}

		for(size_t i=0; i<c.num_const_methods(); i++)
		{
			const bridge::ExportedFunctionBase &thefunc = c.get_const_method(i);
			if(IsClassInParameters(thefunc.params_in(), classname) || IsClassInParameters(thefunc.params_out(), classname))
			{
				UG_LOG(c.name() << " ::");
				PrintFunctionInfo(thefunc, false, classname);
			}
		}
	}
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
	bridge::Registry &reg = GetUGRegistry();
	const IExportedClass *c = NULL;
	// find class
	for(size_t j=0; j<reg.num_classes(); ++j)
		if(strcmp(classname, reg.get_class(j).name()) == 0)
		{
			c = &reg.get_class(j);
			break;
		}
	if(c == NULL)
	{
		UG_LOG("Class name " << classname << " not found\n");
		return false;
	}

	const std::vector<const char*> *names = c->class_names();
	if(names == 0 || names->size() < 2)
		ClassUsageExact(classname);
	else
	{
		for(size_t i = 0; i<names->size(); i++)
		{
			UG_LOG("--- as " << names->at(i) << ": ---\n");
			ClassUsageExact(names->at(i));
		}
	}

	UG_LOG("\n");
	return true;
}

void RegisterInfoCommands(Registry &reg)
{
	reg.add_function("TypeInfo", &UGTypeInfo);
	reg.add_function("ClassUsage", &ClassUsage);
}


} // namespace bridge
} // namespace ug
