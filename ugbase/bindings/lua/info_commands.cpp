/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

/**
 * \file info_commands.cpp
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

#include <iomanip>
#include "bindings/lua/lua_util.h"
#include "bindings/lua/bindings_lua.h"
#include "bridge/bridge.h"
#include "bridge/util.h"
#include "registry/class_helper.h"
#include "common/util/sort_util.h"
#include "common/util/string_util.h"
#include "lua_stack_check.h"
#include "registry/class_name_provider.h"
#include "common/util/string_util.h"
#include "common/util/stringify.h"
#include "common/util/string_table_stream.h"

#include "common/profiler/profiler.h"
#ifdef UG_PLUGINS
	#include "common/util/plugin_util.h"
#endif

#ifdef UG_POSIX
#include <signal.h>
#endif

#ifndef USE_LUAJIT
extern "C" // default lua
{
#include "bindings/lua/externals/lua/lstate.h"
#include "externals/lua/lua.h"
}
#else
// luajit
#include <lua.hpp>
#endif


#include "info_commands.h"


using namespace std;

std::string get_gcc_backtrace();
void lua_backtrace();
void shiny_backtrace();
void ug_backtrace();

namespace ug
{
bool useLua2VM=false;
bool useLuaCompiler=false;

namespace bridge
{

#ifdef UG_POSIX
void (*oldSIGSEGVhandler)(int);
void (*oldSIGINThandler)(int);

void signal_callback_handler(int signum)
{
	cout << "\n##############################################################################################################\n";
	printf("CAUGHT SIGNAL SIGSEV = Segmentation Violation (Invalid access to storage). Possibly an uninitialized pointer.");
	cout << "\n##############################################################################################################\n\n";
	ug_backtrace();
	UG_THROW("CAUGHT SIGNAL SIGSEV = Segmentation Violation (Invalid access to storage). Possibly an uninitialized pointer.");
	exit(signum);
}

void sigint_handler(int signum)
{
	cout << "\n##############################################################################################################\n";
	cout << "--- Received SIGINT (Ctrl+C) --- \n";
	cout << "\n##############################################################################################################\n\n";

	lua_backtrace();
	shiny_backtrace();
	exit(signum);
}

void InitSignals()
{
	oldSIGSEGVhandler = signal(SIGSEGV, signal_callback_handler);
	oldSIGINThandler = signal(SIGINT, sigint_handler);
}
#else
void InitSignals()
{
	UG_LOG("Signals not available on this system (non-posix).\n");
}
#endif

void SetLuaNamespaceInTable(string name, string value)
{
	UG_LOGN("SetLuaNamespaceInTable " << name << " str " << value);
	lua_State* L = script::GetDefaultLuaState();
	LUA_STACK_CHECK(L, 0);
	int dotPos = name.find(".");

	if(dotPos == -1)
	{
		lua_pushstring(L, name.c_str()); // key
		lua_pushstring(L, value.c_str()); // value
		lua_settable(L, -3);
	}
	else
	{
		string subname = name.substr(0, dotPos);
		string restname = name.substr(dotPos+1, name.size());

		lua_pushstring(L, subname.c_str());
		lua_rawget(L, -2);
		if(lua_isnil(L, -1))
		{
			lua_pop(L, 1);

			lua_newtable(L);
			SetLuaNamespaceInTable(restname, value);
			lua_pushstring(L, subname.c_str());
			lua_insert(L, -2); /* swap uppertable and 0 */
			lua_settable(L, -3);
		}
		else
		{
			SetLuaNamespaceInTable(restname.c_str(), value);
			lua_pop(L, 1);
		}
	}
}
void SetLuaNamespace(string name, string value)
{
	UG_LOGN("SetLuaNamespace " << name << " str " << value);
	lua_State* L = script::GetDefaultLuaState();
	LUA_STACK_CHECK(L, 0);

	int dotPos = name.find(".");

	if(dotPos == -1)
	{
		lua_pushstring(L, value.c_str());
		lua_setglobal(L, name.c_str());
	}
	else
	{
		string subname = name.substr(0, dotPos);
		string restname = name.substr(dotPos+1, name.size());
		lua_getglobal(L, subname.c_str());
		if(lua_isnil(L, -1))
		{
			lua_pop(L, 1);
			lua_newtable(L);
			lua_setglobal(L, subname.c_str());
			lua_getglobal(L, subname.c_str());
		}
		SetLuaNamespaceInTable(restname, value);
		lua_pop(L, 1);

	}
}


/**
 * this function also includes all lines before fromline which begin with -- (lua comments)
 * and adds line numbers
 * @param filename
 * @param fromline
 * @param toline
 * @return string with lines
 */
string GetFileLinesLUA(const char *filename, size_t fromline, size_t toline)
{
	if(GetLogAssistant().is_output_process()) return "";
	char buf[512];
	fstream file(filename, ios::in);
	if(file.is_open() == false) return string("");
	stringstream *pss = nullptr;
	for(size_t i=0; i<fromline-1 && !file.eof(); i++)
	{
		file.getline(buf, 512);
		if(strncmp(buf+strspn(buf, "\t "), "--", 2)==0)
		{
			if(pss == nullptr) pss = new stringstream;
			*pss << "   \t" << buf+strspn(buf, "-\t ") << '\n';
		}
		else if(pss) { delete pss; pss = nullptr; }

	}
	stringstream ss;
	if(pss != nullptr) ss << "\n" << pss->str() << "\n";
	for(; fromline <= toline && !file.eof(); fromline++)
	{
		file.getline(buf, 512);
		ss << fromline << "\t" << buf << "\n";
	}
	return ss.str();
}


double LuaGetNumber(lua_State *L, const char *name, double notAvailable)
{
	LUA_STACK_CHECK(L, 0);
	if(GetLuaNamespace(L, name)==false || !lua_isnumber(L, -1))
	{
		lua_pop(L, 1);
		return notAvailable;
	}
	double d = lua_tonumber(L, -1);
	lua_pop(L, 1);
	return d;
}

string LuaGetString(lua_State *L, const char *name, const char *notAvailable)
{
	LUA_STACK_CHECK(L, 0);
	if(GetLuaNamespace(L, name)==false || !lua_isstring(L, -1))
	{
		lua_pop(L, 1);
		return notAvailable;
	}
	string str = lua_tostring(L, -1);
	lua_pop(L, 1);
	return str;
}

bool LuaGetBoolean(lua_State *L, const char *name, bool notAvailable)
{
	LUA_STACK_CHECK(L, 0);
	if(GetLuaNamespace(L, name)==false || !lua_isboolean(L, -1))
	{
		lua_pop(L, 1);
		return notAvailable;
	}
	bool b = lua_toboolean(L, -1);
	lua_pop(L, 1);
	return b;
}




/**
 * \brief searches for a namespaces (lists) and pushes it onto the stack
 * \param L				the lua state
 * \param name			name of the namespace. also nested namespaces are allowed (like struc1.struc2)
 * \return true if namespace found
 */
bool GetLuaNamespace(lua_State* L, string name)
{
	LUA_STACK_CHECK(L, 1);
	vector<string> tokens;
	TokenizeString(name, tokens, '.');
	if(tokens.empty())
	{
		lua_pushnil(L);
		return false;
	}

	lua_getglobal(L, tokens[0].c_str());
	if(lua_isnil(L, -1))
	{
		return 0; 	// global name not found
	}

	size_t i=1;
	for(; i<tokens.size(); i++)
	{
		lua_pushstring(L, tokens[i].c_str());
		lua_rawget(L, -2);
		lua_remove(L, -2); // remove parent table from stack
		if(lua_isnil(L, -1))
			return false;
	}

	return true;
}


const ClassNameNode* GetClassNameNode(lua_State *L, int index)
{
	LUA_STACK_CHECK(L, 0);
	const ClassNameNode* classNameNode = nullptr;
	if(lua_getmetatable(L, index) != 0)
	{
		// get names
		lua_pushstring(L, "class_name_node");
		lua_rawget(L, -2);
		if(!lua_isnil(L, -1) && lua_isuserdata(L, -1))
			classNameNode = (const ClassNameNode*) lua_touserdata(L, -1);
		lua_pop(L, 2); // pop userdata, metatable
	}
	return classNameNode;
}


const std::vector<const char*> *GetClassNames(lua_State *L, int index)
{
	LUA_STACK_CHECK(L, 0);
	const std::vector<const char*> *p = nullptr;
	if(lua_getmetatable(L, index) != 0)
	{
		// get names
		lua_pushstring(L, "__names");
		lua_rawget(L, -2);
		if(!lua_isnil(L, -1) && lua_isuserdata(L, -1))
			p = (const std::vector<const char*>*) lua_touserdata(L, -1);
		lua_pop(L, 2); // pop userdata, metatable
	}
	return p;
}

const std::vector<const char*> *GetClassNames(lua_State* L, const char *name)
{
	// get the lua object with that name
	lua_getglobal(L, name);
	if(lua_isnil(L, -1))
	{
		lua_pop(L, 1);	// remove global from stack
		return nullptr; 	// global name not found
	}

	const std::vector<const char*> *p = GetClassNames(L, -1);
	lua_pop(L, 1); // remove global from stack;
	return p;
}

string GetLUAScriptFunctionDefined(const char *functionName)
{

	lua_State* L = script::GetDefaultLuaState();
	LUA_STACK_CHECK(L, 0);

	string str = functionName;
	GetLuaNamespace(L, str);

	if(lua_isnil(L, -1) || lua_isfunction(L, -1)==false)
	{
		// it is nil
		lua_pop(L, 1);
		return "?";
	}

	lua_Debug ar;

	if(lua_getinfo(L, ">Snlu", &ar) != 0 && ar.source)
	{
		const char *src = ar.source[0]=='@' ? ar.source+1 : ar.source;
		std::stringstream ss;
		ss << src << ":" << ar.linedefined << "-" << ar.lastlinedefined;
		return ss.str();
	}
	else
	{
		return "?";
	}
}

string FunctionInfo(lua_State *L, bool bComplete, const char *functionName)
{
	LUA_STACK_CHECK(L, 0);
	lua_pushvalue(L, -1);
	lua_Debug ar;

	std::stringstream ss;
	if(lua_getinfo(L, ">Snlu", &ar) != 0 && ar.source)
	{
		if(bComplete)
		{
			const char *src = ar.source[0]=='@' ? ar.source+1 : ar.source;
			ss << src << ":" << ar.linedefined << "-" << ar.lastlinedefined << "\n";
			ss << GetFileLinesLUA(src, ar.linedefined, ar.lastlinedefined) << "\n";
		}
		else
		{
			ss << GetFileLine(ar.source+1, ar.linedefined) << "\n";
		}
	}
	else
		ss << functionName;

	return ss.str();
}

string LuaClassMethodInfo(lua_State *L, int index, const ExportedMethod &thefunc)
{
	const std::vector<const char*> *names = GetClassNames(L, index);
	const char *classname = "(unknown class)";
	if(names != nullptr)
		classname = names->at(0);
	return FunctionInfo(thefunc, false, classname);
}


/**
 * \brief Prints info to a lua type
 * \param 	p			the name of the object in lua.
 * you can use class names, function names or the names of an object
 * - TypeInfo("class") prints all member functions+parameters of this class and its parents
 * - TypeInfo("Function") prints all member functions+parameters
 * - TypeInfo("variable") prints class information if variable is a object of a ug class, otherwise what type in lua it is
 */
int UGTypeInfo(const char *p)
{
	UG_LOG("\n");
	const bridge::Registry &reg = GetUGRegistry();

	// check if it is a class
	const ClassGroupDesc *cg = reg.get_class_group(p);
	const IExportedClass *c;
	if(cg != nullptr)
	{
		UG_LOG("ClassGroup " << p << " consisting of classes\n");
		for(size_t i=0; i<cg->num_classes(); i++)
		{
			if(i != 0) UG_LOG(", ");
			UG_LOG(cg->get_class(i)->name());
		}
		c = cg->get_default_class();
		if(c)
		{	UG_LOG("\nDefault class is " << c->name() << "\n");	}
		else
		{
			UG_LOG("\nDefault class (not yet?) set, using first class in list.\n");
			if(cg->num_classes() > 0)
				c = cg->get_class(0);
		}
		UG_LOG("\n\n");
	}
	else
		c = reg.get_class(p);

	if(c)
	{
		const std::vector<const char*> *names = c->class_names();
		for(size_t i=0; i < names->size(); ++i)
			UG_LOG(ClassInfo(reg, names->at(i)));
		UG_LOG(endl);
		UG_LOG(ClassHierarchyString(reg, c->name().c_str()));
		ClassInstantiations(c->name().c_str());
		UG_LOG("\n");
		return true;
	}


	lua_State* L = script::GetDefaultLuaState();
	LUA_STACK_CHECK(L, 0);

	string str = p;
	GetLuaNamespace(L, str);

	if(lua_isnil(L, -1))
	{
		// it is nil
		lua_pop(L, 1);
		UG_LOG(p << " is neither a global variable nor a class name." << endl);
		return false;
	}

	if(lua_iscfunction(L, -1))
	{
		// it is a cfunction
		UG_LOG(p << " is a cfunction\n " << FunctionInfo(reg, p) << "\n");
	}
	else if(lua_isfunction(L, -1))
	{
		// it is a lua function
		UG_LOG(p << " is a function\n " << FunctionInfo(L, true, p) << "\n");
	}
	else if(lua_isuserdata(L, -1))
	{
		// it is user data, that is a class object
		// names = GetClassNames(L, -1))
		if(lua_getmetatable(L, -1) == 0)
		{
			UG_LOG(p << " is a global variable which has light user data, but no metatable." << endl);
			lua_pop(L, 1); // pop globals
			return false;
		}

		// get names
		lua_pushstring(L, "__names");
		lua_rawget(L, -2);
		if(lua_isnil(L, -1) || !lua_isuserdata(L, -1))
		{
			UG_LOG(p << " is a global variable which has a metatable, but cannot access names." << endl);
			lua_pop(L, 3); // pop metatable, userdata, globals
			return false;
		}
		const std::vector<const char*> *names =
				(const std::vector<const char*>*) lua_touserdata(L, -1);
		lua_pop(L, 2); // pop metatable, userdata
		UG_LOG("Typeinfo for " << p << ": " << endl);
		for(size_t i=0; i < names->size(); ++i)
			UG_LOG(ClassInfo(reg, names->at(i)));
		if(!names->empty())
			UG_LOG(ClassHierarchyString(reg, names->at(0)));
	}
	else if (lua_istable(L, -1))
	{
		// it is a table
		LuaPrintTable(L, 1);
		UG_LOG(p << " is a table" << "\n");
	}
	else
	{
		// it is something else like number string or boolean
		UG_LOG(p << ": type is " << GetLuaTypeString(L, -1) << ": " << lua_tostring(L, -1));
	}
	lua_pop(L, 1);

	UG_LOG("\n");
	return 0;
}


/**
 * \brief this function prints all objects which are of a certain class
 * \param 	classname	the class name to be searched
 * \return true if class found, otherwise false
 */
bool ClassInstantiations(const char *classname)
{
	bridge::Registry &reg = GetUGRegistry();
	// search for the class
	const IExportedClass *c = reg.get_class(classname);
	if(c == nullptr)
	{
		UG_LOG("Class " << classname << " not found\n");
		return false;
	}

	UG_LOG(endl);
	UG_LOG("Instantiations of Class " << classname << ":\n");

	lua_State* L = script::GetDefaultLuaState();
	LUA_STACK_CHECK(L, 0);
	bool bFound = false;

#ifndef USE_LUAJIT
	// iterate through all of lua's global string table
	for(int i=0; i<G(L)->strt.size; i++)
	{
		TString *obj;
		for (obj = G(L)->strt.hash[i]; obj != nullptr; obj = obj->u.hnext)
		{
			// get the string
			TString *ts = obj;
			if(ts == nullptr) continue;

			const char *luastr = getstr(ts);
			// check is of a global variable

			const std::vector<const char*> *names = GetClassNames(L, luastr);
			if(names == nullptr)
				continue;

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
#else
	// traversal using API
	// cf. http://stackoverflow.com/questions/20527659/how-to-filter-out-user-defined-globals-in-lua-from-c

	lua_pushvalue(L, LUA_GLOBALSINDEX); 	//lua_pushglobaltable(L);
	lua_pushnil(L);
	while (lua_next(L,-2) != 0)
	{
	  const char* luastr = lua_tostring(L,-2);

	  if (luastr) {
		  std::cerr << "Found global: " << luastr << std::endl;
	  }

	  lua_pop(L,1); // pop value
	}
	lua_pop(L,1); // pop global table

	// Check required!
	UG_ASSERT(0, "ERROR: Implement for LuaJit!");

#endif
	if(!bFound) UG_LOG("No instantiations of " << classname << " or subclasses found.");
	UG_LOG(endl);
	return true;
}

/**
 *
 * \param classname the class to print usage in functions/member functions of (and all its subclasses) .
 * class in in/out parameters is highlighted with [class].
 * \return true if class found, otherwise fals
 */
string ClassUsage(const char *classname)
{
	bridge::Registry &reg = GetUGRegistry();
	std::stringstream ss;
	ss << "\n";

	// find class
	const IExportedClass *c = reg.get_class(classname);
	if(c == nullptr)
	{
		ss << "Class name " << classname << " not found\n";
		return ss.str();
	}

	// print usages in functions

	ss << "--- Functions returning " << classname << ": ---\n";
	ss << ClassUsageExact(reg, classname, true);

	const std::vector<const char*> *names = c->class_names();
	if(names != nullptr && !names->empty())
	{
		for(size_t i = 0; i<names->size(); i++)
		{
			ss << "--- Functions using " << names->at(i) << ": ---\n";
			ss << ClassUsageExact(reg, names->at(i), false);
		}
	}

	ss << ClassInstantiations(classname);

	ss << "\n";
	return ss.str();
}

string LuaGetScriptFunctionString(lua_State *L, int index)
{
	LUA_STACK_CHECK(L, 0);
	lua_pushvalue(L, index);
	lua_Debug ar;

	if(lua_getinfo(L, ">S", &ar) != 0 && ar.linedefined != -1)
	{

		const char *p=GetFileLine(ar.source[0] == '@' ? ar.source+1 : ar.source,
								  ar.linedefined).c_str();
		p+=strspn(p, " \t");
		return p;
	}
	return "";
}

/**
 * \brief prints the source of a lua script function which is on top of the stack
 * \param L		the lua state
 */
void LuaPrintTable(lua_State *L, size_t iSpace, int index)
{
	LUA_STACK_CHECK(L, 0);

	if(index != -1)
		lua_pushvalue(L, index);

	UG_LOG(repeat(' ', iSpace));
	UG_LOG("{\n");
	//lua_getglobal(L, "ugargv"); // -2
	int len=0;
	lua_pushnil( L ); // -1
	while( lua_next( L, (-1-len)*2) )
	{
		// key -2
		// val -1
		lua_pushvalue(L, -2);
		len++;
	}

	std::vector<SortStruct<int, string> > sorted;
	sorted.resize(len);

	for(int i=0; i<len; i++)
	{
		sorted[i].index = -2*len+2*i+1; // valueIndex
		size_t indexIndex = -2*len+2*i;
		const char *name = lua_tostring(L, indexIndex);
		if(name) sorted[i].value = name;
		else sorted[i].value = GetLuaTypeString(L, indexIndex);
	}

	sort(sorted.begin(), sorted.end());

	for(int i=0; i<len; i++)
	{
		int index2 = sorted[i].index;

		UG_LOG(repeat(' ', iSpace));
		UG_LOG(sorted[i].value);
		UG_LOG(" (" << GetLuaTypeString(L, index2) << ")");
		if(lua_isfunction(L, index2))
		{
			UG_LOG(": " << LuaGetScriptFunctionString(L, index2));
		}
		else if(lua_istable(L, index2))
		{
			UG_LOG(" = \n");
			LuaPrintTable(L, iSpace+1, index2);
		}
		else
		{
			const char * value = lua_tostring(L, index2);
			if(value) { UG_LOG(" = \"" << value << "\"") };
		}
		UG_LOG("\n");
	}

	lua_pop(L, 2*len);
	UG_LOG(repeat(' ', iSpace));
	UG_LOG("}\n");
	if(index != -1)
		lua_pop(L, 1);
}


/**
 * @brief function to list all objects of the default LUA state
 * @param functions
 * @param internalFunctions
 * @param scriptFunctions
 * @param luaObjects
 * @param classInstantiations
 */
void GetLuaGlobals(std::vector<std::string>  *functions,
		std::vector<std::string> *internalFunctions,
		std::vector<std::string> *scriptFunctions,
		std::vector<std::string> *luaObjects,
		std::vector<std::string> *classInstantiations)
{
	bridge::Registry &reg = GetUGRegistry();
	lua_State* L = script::GetDefaultLuaState();
	LUA_STACK_CHECK(L, 0);

	// iterate through all of lua's globals

	lua_getglobal(L, "_G");
	lua_pushnil( L ); // -1
	while( lua_next( L, -2 ))
	{
		const char *luastr = lua_tostring(L, -2);
		if(luastr && strcmp(luastr, "_G") != 0 && strcmp(luastr, "package") != 0)
		{
			if(!reg.get_class(luastr))
			{
				if(FindFunction(reg, luastr))
				{
					if(functions) functions->push_back(luastr);
				}
				else if(lua_isfunction(L, -1) || lua_iscfunction(L, -1))
				{
					lua_Debug ar;
					lua_pushvalue(L, -1);
					if(lua_getinfo(L, ">S", &ar) != 0 && ar.linedefined != -1) {
						if(scriptFunctions) scriptFunctions->push_back(luastr);
					}
					else {
						if(internalFunctions) internalFunctions->push_back(luastr);
					}
				}
				else if(lua_isuserdata(L, -1)) {
					if(classInstantiations) classInstantiations->push_back(luastr);
				}
				else {
					if(luaObjects) luaObjects->push_back(luastr);
				}
			}
		}
		lua_pop(L, 1); // remove table entry from stack
	}
	lua_pop(L, 1); // remove global _G from stack

	if(functions) sort(functions->begin(), functions->end());
	if(internalFunctions) sort(internalFunctions->begin(), internalFunctions->end());
	if(scriptFunctions) sort(scriptFunctions->begin(), scriptFunctions->end());
	if(luaObjects) sort(luaObjects->begin(), luaObjects->end());
	if(classInstantiations) sort(classInstantiations->begin(), classInstantiations->end());
}

void GetLuaGlobal_functions(std::vector<std::string> &functions)
{
	GetLuaGlobals(&functions, nullptr, nullptr, nullptr, nullptr);
}

void GetLuaGlobal_internalFunctions(std::vector<std::string> &internalFunctions)
{
	GetLuaGlobals(nullptr, &internalFunctions, nullptr, nullptr, nullptr);
}

void GetLuaGlobal_scriptFunctions(std::vector<std::string> &scriptFunctions)
{
	GetLuaGlobals(nullptr, nullptr, &scriptFunctions, nullptr, nullptr);
}

void GetLuaGlobal_luaObjects(std::vector<std::string> &luaObjects)
{
	GetLuaGlobals(nullptr, nullptr, nullptr, &luaObjects, nullptr);
}

void GetLuaGlobal_classInstantiations(std::vector<std::string> &classInstantiations)
{
	GetLuaGlobals(nullptr, nullptr, nullptr, nullptr, &classInstantiations);
}

void LuaList_classes()
{
	bridge::Registry &reg = GetUGRegistry();
	std::vector<std::string> classes;
	classes.reserve(reg.num_classes());
	for(size_t j=0; j<reg.num_classes(); ++j)
		classes.push_back(reg.get_class(j).name());
	sort(classes.begin(), classes.end());

	UG_LOG(endl << "--- Classes: --------------------" << endl)
	for(size_t i=0; i<classes.size(); i++)
		UG_LOG(classes[i] << endl);
}


void LuaList_cfunctions()
{
	bridge::Registry &reg = GetUGRegistry();
	std::vector<std::string> functions;
	GetLuaGlobal_functions(functions);
	UG_LOG(endl << "--- C Functions: ------------------" << endl)
	for(size_t i=0; i<functions.size(); i++)
		UG_LOG(FunctionInfo(reg, functions[i].c_str()) << "\n");
}

void LuaList_scriptFunctions()
{
	lua_State* L = script::GetDefaultLuaState();

	std::vector<std::string> scriptFunctions;
	GetLuaGlobal_scriptFunctions(scriptFunctions);
	UG_LOG(endl << "--- Script Functions: ---" << endl)

	if(scriptFunctions.empty())	return;
	int maxLength = (*max_element(scriptFunctions.begin(), scriptFunctions.end(), IsLonger)).size();
	for(size_t i=0; i<scriptFunctions.size(); i++)
	{
		lua_getglobal(L, scriptFunctions[i].c_str());  /* get global 'f' */
		UG_LOG(left << setw(maxLength) << scriptFunctions[i] << ": ");
		UG_LOG(LuaGetScriptFunctionString(L));
		lua_pop(L, 1);
		UG_LOG("\n");
	}
}

void LuaList_internalFunctions()
{
	std::vector<std::string> internalFunctions;
	GetLuaGlobal_internalFunctions(internalFunctions);

	UG_LOG(endl << "--- Internal Functions: ---" << endl)
	for(size_t i=0; i<internalFunctions.size(); i++)
		UG_LOG(internalFunctions[i] << endl);
}


void LuaList_luaObjects()
{
	lua_State* L = script::GetDefaultLuaState();
	std::vector<std::string> luaObjects;
	GetLuaGlobal_luaObjects(luaObjects);

	UG_LOG(endl << "--- Lua Objects: ----------------" << endl)
	if(luaObjects.empty()) return;
	int maxLength = (*max_element(luaObjects.begin(), luaObjects.end(), IsLonger)).size();
	for(size_t i=0; i<luaObjects.size(); i++)
	{
		if(luaObjects[i].compare("_G") == 0) continue;
		if(luaObjects[i] == "package") continue;
		lua_getglobal(L, luaObjects[i].c_str());
		UG_LOG(left << setw(maxLength) << luaObjects[i]);
		UG_LOG(" (" << GetLuaTypeString(L, -1) << ")");
		/*if(lua_istable(L, -1))
		{
			UG_LOG(" = \n");
			LuaPrintTable(L, 1);
		}*/
		const char *p = lua_tostring(L, -1);
		if(p) UG_LOG(": " << p)
		lua_pop(L, 1);
		UG_LOG(endl);
	}
}

void LuaList_classInstantiations()
{
	lua_State* L = script::GetDefaultLuaState();
	std::vector<std::string> instantiations;
	GetLuaGlobal_classInstantiations(instantiations);

	UG_LOG(endl << "--- Class Instantiations: ---------" << endl)
	if(instantiations.empty()) return;
	int maxLength = (*max_element(instantiations.begin(), instantiations.end(), IsLonger)).size();
	for(size_t i=0; i<instantiations.size(); i++)
	{
		lua_getglobal(L, instantiations[i].c_str());

		const char * type = "UNKNOWN";
		int ref = 0;
		void* ptr = lua_touserdata(L, 1);
		
	//	we perform delete if the user-data is a raw pointer
		if(((lua::UserDataWrapper*)ptr)->is_raw_ptr()){
			type = "raw ptr ";
		}
		else if(((lua::UserDataWrapper*)ptr)->is_smart_ptr()){

		//	invalidate the associated smart-pointer
			if(((lua::UserDataWrapper*)ptr)->is_const())
			{
				type = "ConstSmartPtr ";
				ref = ((lua::ConstSmartUserDataWrapper*)ptr)->smartPtr.refcount();
			}
			else
			{
				type = "SmartPtr ";
				ref =  ((lua::SmartUserDataWrapper*)ptr)->smartPtr.refcount();
			}
		}

		if(!lua_isuserdata(L, -1)) continue;
		const std::vector<const char*> *n  = GetClassNames(L, -1);
		if(n && !n->empty())
		{
			UG_LOG(left << setw(maxLength) << instantiations[i]);
			if(ref != 0)
			{	UG_LOG(" refcount = " << ref << "\t"); }
			UG_LOG(type);
//			for(size_t j = 0; j < n->size(); j++)
//			{
//				if(j > 0) UG_LOG(", ");
//				UG_LOG(n->at(j));
//			}
//			UG_LOG(endl);
			UG_LOG(n->at(0));

			UG_LOG("\n");
		}
		lua_pop(L, 1);
	}
}


void LuaList()
{
	LuaList_cfunctions();
	LuaList_classes();
	LuaList_internalFunctions();
	LuaList_luaObjects();
	LuaList_scriptFunctions();
}

string GetLuaTypeString(lua_State* L, int index)
{
	if(lua_isnil(L, index))
		return string("nil");
	string str("");
	// somehow lua_typeinfo always prints userdata
	if(lua_isboolean(L, index)) str.append("boolean/");
	if(lua_iscfunction(L, index)) str.append("cfunction/");
	if(lua_isfunction(L, index)) str.append("function/");
	if(lua_islightuserdata(L, index)) str.append("lightuserdata/");
	if(lua_isnil(L, index)) str.append("nil/");
	if(lua_isnone(L, index)) str.append("none/");
	if(lua_isnumber(L, index)) 	str.append("number/");
	if(lua_isstring(L, index)) str.append("string/");

	if(lua_istable(L, index)) {
		//lua_pushnil(L);
		//str.append("{");
		/*bool bFirst = true;
		while (lua_next(L, index) != 0) {
			if(bFirst) {bFirst = false;} else {str.append(", ");};
			str.append(GetLuaTypeString(L, -1));
			lua_pop(L, 1);
	   }*/
		str.append("table/");
	}
	if(lua_isthread(L, index)) str.append("thread/");
	if(lua_isuserdata(L, index))
	{
		if(((lua::UserDataWrapper*)lua_touserdata(L, index))->is_const()){
			str.append("const ");
		}
		const ClassNameNode* classNameNode = GetClassNameNode(L, index);
		if(classNameNode == nullptr || classNameNode->empty()) str.append("userdata/");
		else str.append(classNameNode->name());
		str.append("*/");
	}

	if(lua_type(L, index) == LUA_TNONE)	str.append("none/");

	if(str.size() == 0)
		return string("unknown type");
	else
		return str.substr(0, str.size()-1);
}


/**
 * @param L
 * @param entry (out) returns the lua_Debug entry associated with the deepest call stack level (the "current line")
 */
void LuaGetLastLine(lua_State* L, lua_Debug entry)
{
    for(int depth = 0; lua_getstack(L, depth, &entry); depth++)
	{
    	int status = lua_getinfo(L, "Sln", &entry);
    	if(!status || entry.currentline < 0) continue;
    }
}

/**
 * @param L
 * @return the current line with content
 * example:
 * @/Users/mrupp/Documents/workspace/ug4svn/apps/amg//setup.lua:576        local dim = p.approxSpace:get_dim()
 */
string LuaCurrentLine(lua_State* L)
{
	std::stringstream ss;
	lua_Debug entry;
	LuaGetLastLine(L, entry);
	ss << entry.short_src << ":" << entry.currentline;
	ss << " " << GetFileLine(entry.short_src, entry.currentline) << "\n";
	return ss.str();

}

/**
 * this function returns e.g.
 * 1  @/Users/mrupp/Documents/workspace/ug4svn/apps/amg//setup.lua:576        local dim = p.approxSpace:get_dim()
 * 2  @/Users/mrupp/Documents/workspace/ug4svn/apps/amg//setup.lua:649        amgsetup.CreateApproxSpaceAndDisc()
 * 3  @/Users/mrupp/Documents/workspace/ug4svn/apps/amg/famg_laplace.lua:21   amgsetup.LoadStdDomAndApproxSpace()
 * @param L
 * @param fromLevel  where to start
 * @param toLevel  how far we want to go back
 * @return a list list
 */
string LuaStackTraceString(lua_State* L, int fromLevel, int toLevel)
{
	StringTableStream sts;
    lua_Debug entry;

    std::vector<string> filenames;

    //sts.table().set_col_alignments("ll");
    int luaLevel=0;
    for(int depth = 0; lua_getstack(L, depth, &entry); depth++)
	{
		int status = lua_getinfo(L, "Sln", &entry);
		if(entry.currentline <0) continue;
		if(!status || !entry.source || entry.currentline < 0) continue;
		luaLevel++;
		if(luaLevel < fromLevel) continue;
		// first col
		sts << (Stringify() << " " << luaLevel << "  " << entry.source << ":" << entry.currentline).str();
		// second col
		sts << TrimString(GetFileLine(entry.source, entry.currentline));
		// next row
		sts << "\n";
		if(luaLevel == toLevel) break;
	}

    return sts.to_string();
}

/// returns the current file and line ( \sa LuaStackTrace ).
bool GetLuaFileAndLine(lua_State* L, std::string &file, size_t &line)
{
	file = "unknown";
	line = -1;
	lua_Debug entry;
	for(int depth = 0; lua_getstack(L, depth, &entry); depth++)
	{
	   	int status = lua_getinfo(L, "Sln", &entry);
	   	if(!status || !entry.source || entry.currentline <0) continue;
		file = entry.source;
		line = entry.currentline;
		return true;
	}
	return false;
}

/// returns the current file and line ( \sa LuaStackTrace ).
std::string GetLuaFileAndLine(lua_State* L)
{
	PROFILE_FUNC();
	string file; size_t line;
	if(GetLuaFileAndLine(L, file, line) == false) return "[LUA File could not be determined]";
	std::stringstream ss;
	if(GetLogAssistant().is_output_process())
		ss << file << ":" << line << " " << GetFileLine(file.c_str(), line);
	else
		ss << file << ":" << line;

	return ss.str();
}

std::string GetLuaFileAndLineNumber(lua_State* L)
{
	string file; size_t line;
	if(GetLuaFileAndLine(L, file, line) == false) return "[ --unknown file -- ]";
	std::stringstream ss;
	ss << file << ":" << line;
	return ss.str();
}

std::string GetLuaLine(lua_State* L)
{
	PROFILE_FUNC();
	string file; size_t line;
	if(GetLuaFileAndLine(L, file, line) == false) return "[ --unknown script line -- ]";
	if(GetLogAssistant().is_output_process())
		return GetFileLine(file.c_str(), line);
	else
		return Stringify() << file << ":" << line;
}

std::string LuaStackTraceString()
{
	return LuaStackTraceString(script::GetDefaultLuaState(), 0, -1);
}

void LuaStackTrace(int fromLevel)
{
	UG_LOG(LuaStackTraceString(script::GetDefaultLuaState(), fromLevel, -1));
}

void ScriptPrintClassHierarchy(const char *classname)
{
	UG_LOG(ClassHierarchyString(GetUGRegistry(), classname));
}

	
bool ScriptHasClass(const char *classname)
{
	const Registry &reg = GetUGRegistry();
	return reg.get_class(classname) != nullptr;
}

bool ScriptHasClassGroup(const char *classname)
{
	const Registry &reg = GetUGRegistry();
	return reg.get_class_group(classname) != nullptr;
}

void ScriptPrintClassUsage(const char *classname)
{
	UG_LOG(ClassUsage(classname));
}

#ifdef UG_PLUGINS
bool PluginRequired(const char *name)
{
	if(PluginLoaded(name) == false)
	{
		string msg = string("plugin ") + name + string(" not loaded. Please use 'cmake -D") + name + string("=ON .' in your build directory.");		
		std::string file; size_t line;
		if(GetLuaFileAndLine(script::GetDefaultLuaState(), file, line))
			throw UGError(msg.c_str(), file.c_str(), line);
		else
			throw UGError(msg.c_str());
		return false;
	}
	return true;
}
#endif

void EnableLUA2C(bool b)
{
#ifndef USE_LUA2C
	UG_LOG("Warning: LUA2C not enabled. Enable with \"cmake -DUSE_LUA2C=ON ..\"\n")
#else
	useLuaCompiler=b;
	useLua2VM=false;
#endif
}

void EnableLUA2VM(bool b)
{
	useLuaCompiler=b;
	useLua2VM=b;
}

bool RegisterSerializationCommands(Registry &reg, const char* parentGroup);

bool RegisterInfoCommands(Registry &reg, const char* parentGroup)
{
	RegisterSerializationCommands(reg, parentGroup);
	stringstream grpSS; grpSS << parentGroup << "/Info";
	std::string grp = grpSS.str();

	try
	{
		reg.add_function("ls", &LuaList, grp.c_str(), 
		                 "", "", "list all objects");
		reg.add_function("list_cfunctions", &LuaList_cfunctions, grp.c_str(), 
		                 "", "", "list all cfunctions");
		reg.add_function("list_classes", &LuaList_classes, grp.c_str(), 
		                 "", "", "list all classes");
		reg.add_function("list_internalFunctions", &LuaList_internalFunctions, grp.c_str(), 
		                 "", "", "list all of LUAs internal functions");
		reg.add_function("list_luaObjects", &LuaList_luaObjects, grp.c_str(), 
		                 "", "", "list all created LUA objects");
		reg.add_function("list_scriptFunctions", LuaList_scriptFunctions, grp.c_str(), 
		                 "", "", "list all LUA script functions");
		reg.add_function("list_objects", LuaList_classInstantiations, grp.c_str(),
				                 "", "", "list all LUA class objects");

		reg.add_function("TypeInfo", &UGTypeInfo, grp.c_str(), 
		                 "", "typeName", "print information about a type");
		reg.add_function("ClassUsage", &ScriptPrintClassUsage, grp.c_str(),
		                 "", "typeName", "print information about the usage of a type");
		reg.add_function("ClassInstantiations" ,&ClassInstantiations, grp.c_str(), 
		                 "", "typeName", "print all objects of the type");
		reg.add_function("ClassHierarchy" ,&ScriptPrintClassHierarchy, grp.c_str(), 
		                 "", "typeName", "print the class hierachy of type");
		reg.add_function("Stacktrace", &LuaStackTrace, grp.c_str(),
		                 "", "", "prints the LUA function stack, that is which functions are called up to this point");
		reg.add_function("HasClass", &ScriptHasClass, grp.c_str(), 
		                 "true if class exists", "className", "use only if you know that you're not using a class group, otherwise HasClassGroup");
		reg.add_function("HasClassGroup", &ScriptHasClassGroup, grp.c_str(), 
		                 "true if class oder classGroup exists", "classGroupName", "can be used before instantiating a class");
#ifdef UG_PLUGINS
		reg.add_function("PluginLoaded", &PluginLoaded, grp.c_str(), 
		                 "true if plugin loaded", "pluginName", "pluginName as listed when using cmake ..");

		reg.add_function("PluginRequired", &PluginRequired, grp.c_str(),
		                 "true if plugin loaded", "pluginName", "throws an error if plugin not loaded, displays help string how to enable plugins via cmake -DpluginName=ON ..");
		reg.add_function("GetLoadedPlugins", &GetLoadedPlugins, grp.c_str(), 
		                 "list of loaded plugins names", "", "");
#endif
		reg.add_function("EnableLUA2C", &EnableLUA2C, grp.c_str(), 
		                 "", "bEnable", "");
		reg.add_function("EnableLUA2VM", &EnableLUA2VM, grp.c_str(),
				"", "bEnable", "");
		reg.add_function("InitSignals", &InitSignals, grp.c_str());
	}
	UG_REGISTRY_CATCH_THROW(grp);

	return true;
}


} // namespace bridge

} // namespace ug
