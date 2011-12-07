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
#include "registry/class_helper.h"
#include "common/util/sort_util.h"
#include "common/util/string_util.h"
#include "lua_stack_check.h"
#include "registry/class_name_provider.h"

extern "C"
{
#include "bindings/lua/externals/lua/lstate.h"
}

#include "info_commands.h"


using namespace std;

namespace ug
{
namespace bridge
{

/**
 * \brief returns an integer to a lua-variable.
 * \param L				the lua state
 * \param name			name of the variable. namespaces are possible (like math.pi)
 * \param notAvailable	return value if variable was not found
  */
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

/**
 * \brief returns the string to a lua-variable.
 * \param L				the lua state
 * \param name			name of the variable. namespaces are possible (like math.pi)
 * \param notAvailable	return value if variable was not found
  */
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

/**
 * \brief returns a boolean to a lua-variable.
 * \param L				the lua state
 * \param name			name of the variable. namespaces are possible (like math.pi)
 * \param notAvailable	return value if variable was not found
  */
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



string GetFileLine(const char *filename, size_t line);
string GetFileLines(const char *filename, size_t fromline, size_t toline, bool includeLineNumbers=false);
void LuaPrintTable(lua_State *L, size_t iSpace, int index=-1);
bool ClassInstantiations(const char *classname);

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
	if(tokens.size() == 0)
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
		lua_remove(L, -2);
		if(lua_isnil(L, -1))
			return false;
	}

	return true;
}

/**
 * \param L 		the lua state
 * \param index		the index of the lua object
 * \return the corresponing ClassNameNode
 */
const ClassNameNode* GetClassNameNode(lua_State *L, int index)
{
	LUA_STACK_CHECK(L, 0);
	const ClassNameNode* classNameNode = NULL;
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


/**
 * \param L 		the lua state
 * \param index		the index of the lua object
 * \return a std::vector of class names
 */
const std::vector<const char*> *GetClassNames(lua_State *L, int index)
{
	LUA_STACK_CHECK(L, 0);
	const std::vector<const char*> *p = NULL;
	if(lua_getmetatable(L, index) != 0)
	{
		// get names
		lua_pushstring(L, "names");
		lua_rawget(L, -2);
		if(!lua_isnil(L, -1) && lua_isuserdata(L, -1))
			p = (const std::vector<const char*>*) lua_touserdata(L, -1);
		lua_pop(L, 2); // pop userdata, metatable
	}
	return p;
}

/**
 * \param L 		the lua state
 * \param name		the name of the lua object
 * \return the class names of the object or NULL if not found
 */
const std::vector<const char*> *GetClassNames(lua_State* L, const char *name)
{
	// get the lua object with that name
	lua_getglobal(L, name);
	if(lua_isnil(L, -1))
	{
		lua_pop(L, 1);	// remove global from stack
		return NULL; 	// global name not found
	}

	const std::vector<const char*> *p = GetClassNames(L, -1);
	lua_pop(L, 1); // remove global from stack;
	return p;
}

/**
 * \brief prints out information for a lua function (a function defined in lua script)
 * \param L 			the lua state
 * \param bComplete		if complete, we print the source code of the function
 * \return 0
 */
int PrintFunctionInfo(lua_State *L, bool bComplete)
{
	LUA_STACK_CHECK(L, 0);
	lua_pushvalue(L, -1);
	lua_Debug ar;
	lua_getinfo(L, ">Snlu", &ar);
	if(ar.source)
	{
		if(bComplete)
		{
			UG_LOG(ar.source+1 << ":" << ar.linedefined << "-" << ar.lastlinedefined << "\n");
			UG_LOG(GetFileLines(ar.source+1, ar.linedefined, ar.lastlinedefined, true) << "\n");
		}
		else
		{
			UG_LOG(GetFileLine(ar.source+1, ar.linedefined) << "\n");
		}
	}

	return 0;
}

/**
 * \brief prints out information for a method of a class
 * \param L 			the lua state
 * \param index			index of the class object on the lua stack
 * \param thefunc		method to print
 * \return 0
 */
void PrintLuaClassMethodInfo(lua_State *L, int index, const ExportedMethod &thefunc)
{
	const std::vector<const char*> *names = GetClassNames(L, index);
	const char *classname = "(unknown class)";
	if(names != NULL)
		classname = names->at(0);
	PrintFunctionInfo(thefunc, false, classname);
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
	bridge::Registry &reg = GetUGRegistry();

	// check if it is a class
	const IExportedClass *c = FindClass(reg, p);
	if(c)
	{
		const std::vector<const char*> *names = c->class_names();
		for(size_t i=0; i < names->size(); ++i)
			PrintClassInfo(reg, names->at(i));
		UG_LOG(endl);
		PrintClassHierarchy(reg, c->name().c_str());
		ClassInstantiations(c->name().c_str());
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
		UG_LOG(p << " is a cfunction\n ");
		PrintFunctionInfo(reg, p);
		UG_LOG(endl);
	}
	else if(lua_isfunction(L, -1))
	{
		// it is a lua function
		UG_LOG(p << " is a function\n ");
		PrintFunctionInfo(L, true);
		UG_LOG(endl);
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
		lua_pushstring(L, "names");
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
			PrintClassInfo(reg, names->at(i));
		if(names->size() > 0)
			PrintClassHierarchy(reg, names->at(0));
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
	const IExportedClass *c = FindClass(reg, classname);
	if(c == NULL)
	{
		UG_LOG("Class " << classname << " not found\n");
		return false;
	}

	UG_LOG(endl);
	UG_LOG("Instantiations of Class " << classname << ":\n");

	lua_State* L = script::GetDefaultLuaState();
	LUA_STACK_CHECK(L, 0);
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

			const std::vector<const char*> *names = GetClassNames(L, luastr);
			if(names == NULL)
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
bool ClassUsage(const char *classname)
{
	bridge::Registry &reg = GetUGRegistry();
	UG_LOG("\n");

	// find class
	const IExportedClass *c = FindClass(reg, classname);
	if(c == NULL)
	{
		UG_LOG("Class name " << classname << " not found\n");
		return false;
	}

	// print usages in functions

	UG_LOG("--- Functions returning " << classname << ": ---\n");
	ClassUsageExact(reg, classname, true);

	const std::vector<const char*> *names = c->class_names();
	if(names != NULL && names->size() > 0)
	{
		for(size_t i = 0; i<names->size(); i++)
		{
			UG_LOG("--- Functions using " << names->at(i) << ": ---\n");
			ClassUsageExact(reg, names->at(i), false);
		}
	}

	ClassInstantiations(classname);

	UG_LOG("\n");
	return true;
}

static void LOGFillSpace(size_t n)
{
	for(size_t i=0; i<n; i++)
		UG_LOG(" ");
}

bool IsLonger(const std::string &a, const std::string &b)
{
	return b.size() > a.size();
}


/**
 * \brief function to get some lines of a file
 * \param filename				file name
 * \param fromline
 * \param toline
 * \param includeLineNumbers	if true, add the line number in front of each line and a tab.
 *
 * \return lines fromline to toline of file filename.
 */
string GetFileLines(const char *filename, size_t fromline, size_t toline, bool includeLineNumbers)
{
	char buf[512];
	fstream file(filename, ios::in);
	if(file.is_open() == false) return string("");
	for(size_t i=0; i<fromline; i++)
		file.getline(buf, 512);
	stringstream ss;
	if(includeLineNumbers)
		ss << fromline << "\t";
	ss << buf;
	for(; fromline < toline; fromline++)
	{
		file.getline(buf, 512);
		ss << "\n";
		if(includeLineNumbers)
			ss << fromline << "\t";
		ss << buf;
	}
	return ss.str();
}

/**
 * \brief function to get a line of a file
 * \param filename				file name
 * \param line
 * \return the line of the file.
 */
string GetFileLine(const char *filename, size_t line)
{
	return GetFileLines(filename, line, line, false);
}


/**
 * \brief prints the source of a lua script function which is on top of the stack
 * \param L		the lua state
  */
void PrintLuaScriptFunction(lua_State *L, int index=-1)
{
	LUA_STACK_CHECK(L, 0);
	lua_pushvalue(L, index);
	lua_Debug ar;
	lua_getinfo(L, ">S", &ar);
	if(ar.linedefined != -1)
		UG_LOG(GetFileLine(ar.source+1, ar.linedefined)); // skip '@'
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

	LOGFillSpace(iSpace);
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
		int index = sorted[i].index;

		LOGFillSpace(iSpace);
		UG_LOG(sorted[i].value);
		UG_LOG(" (" << GetLuaTypeString(L, index) << ")");
		if(lua_isfunction(L, index))
		{
			UG_LOG(": ");
			PrintLuaScriptFunction(L, index);
		}
		else if(lua_istable(L, index))
		{
			UG_LOG(" = \n");
			LuaPrintTable(L, iSpace+1, index);
		}
		else
		{
			const char * value = lua_tostring(L, index);
			if(value) { UG_LOG(" = \"" << value << "\"") };
		}
		UG_LOG("\n");
	}

	lua_pop(L, 2*len);
	LOGFillSpace(iSpace);
	UG_LOG("}\n");
	if(index != -1)
		lua_pop(L, 1);
}





/**
 * \brief prints information about lua variables and classes
  */
void LuaList()
{
	bridge::Registry &reg = GetUGRegistry();
	lua_State* L = script::GetDefaultLuaState();
	LUA_STACK_CHECK(L, 0);
	std::vector<std::string> classes, functions, internalFunctions, scriptFunctions, names, instantiations;
	// iterate through all of lua's globals

	lua_getglobal(L, "_G");
	lua_pushnil( L ); // -1
	while( lua_next( L, -2 ))
	{
		const char *luastr = lua_tostring(L, -2);
		if(luastr && strcmp(luastr, "_G") != 0 && strcmp(luastr, "package") != 0)
		{
			if(FindClass(reg, luastr))
				classes.push_back(luastr);
			else if(FindFunction(reg, luastr))
				functions.push_back(luastr);
			else if(lua_isfunction(L, -1) || lua_iscfunction(L, -1))
			{
				lua_Debug ar;
				lua_pushvalue(L, -1);
				lua_getinfo(L, ">S", &ar);
				if(ar.linedefined != -1)
					scriptFunctions.push_back(luastr);
				else
					internalFunctions.push_back(luastr);
			}
			else if(lua_isuserdata(L, -1))
				instantiations.push_back(luastr);

				names.push_back(luastr);
		}
		lua_pop(L, 1); // remove table entry from stack
	}
	lua_pop(L, 1); // remove global _G from stack

	classes.clear();
	for(size_t j=0; j<reg.num_classes(); ++j)
		classes.push_back(reg.get_class(j).name());
	sort(classes.begin(), classes.end());
	sort(functions.begin(), functions.end());
	sort(names.begin(), names.end());
	sort(scriptFunctions.begin(), scriptFunctions.end());
	sort(internalFunctions.begin(), internalFunctions.end());
	sort(instantiations.begin(), instantiations.end());

	UG_LOG(endl << "--- Classes: --------------------" << endl)
	for(size_t i=0; i<classes.size(); i++)
		UG_LOG(classes[i] << endl);


	UG_LOG(endl << "--- Functions: ------------------" << endl)
	for(size_t i=0; i<functions.size(); i++)
	{
		if(PrintFunctionInfo(reg, functions[i].c_str()) == false)
			UG_LOG(functions[i]);
		UG_LOG(endl);
	}

	UG_LOG(endl << "--- Script Functions: ---" << endl)

	int maxLength = (*max_element(scriptFunctions.begin(), scriptFunctions.end(), IsLonger)).size();
	for(size_t i=0; i<scriptFunctions.size(); i++)
	{
		lua_getglobal(L, scriptFunctions[i].c_str());  /* get global 'f' */
		UG_LOG(left << setw(maxLength) << scriptFunctions[i] << ": ");
		PrintLuaScriptFunction(L);
		lua_pop(L, 1);
		UG_LOG("\n");
	}

	UG_LOG(endl << "--- Internal Functions: ---" << endl)

	for(size_t i=0; i<internalFunctions.size(); i++)
		UG_LOG(internalFunctions[i] << endl);



	UG_LOG(endl << "--- Lua Objects: ----------------" << endl)
	maxLength = (*max_element(names.begin(), names.end(), IsLonger)).size();
	for(size_t i=0; i<names.size(); i++)
	{
		if(names[i].compare("_G") == 0) continue;
		if(names[i].compare("package") == 0) continue;
		lua_getglobal(L, names[i].c_str());
		UG_LOG(left << setw(maxLength) << names[i]);
		UG_LOG(" (" << GetLuaTypeString(L, -1) << ")");
		if(lua_istable(L, -1))
		{
			UG_LOG(" = \n");
			LuaPrintTable(L, 1);
		}
		const char *p = lua_tostring(L, -1);
		if(p) UG_LOG(": " << p)
		lua_pop(L, 1);
		UG_LOG(endl);
	}

	UG_LOG(endl << "--- Class Instantiations: ---------" << endl)
	maxLength = (*max_element(instantiations.begin(), instantiations.end(), IsLonger)).size();
	for(size_t i=0; i<instantiations.size(); i++)
	{
		lua_getglobal(L, instantiations[i].c_str());
		if(!lua_isuserdata(L, -1)) continue;
		const std::vector<const char*> *n  = GetClassNames(L, -1);
		if(n && n->size() > 0)
		{
			UG_LOG(left << setw(maxLength) << instantiations[i] << ": class ");
			for(size_t j = 0; j < n->size(); j++)
			{
				if(j > 0) UG_LOG(", ");
				UG_LOG(n->at(j));
			}
			UG_LOG(endl);
		}
		lua_pop(L, 1);
	}
}


/**
 * \brief returns a String describing the content of the lua stack at a given index
 * and all types it is compatible with (for ex. "2" is string and number)
 */
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

	if(lua_istable(L, index)) str.append("table/");
	if(lua_isthread(L, index)) str.append("thread/");
	if(lua_isuserdata(L, index))
	{
		if(((lua::UserDataWrapper*)lua_touserdata(L, index))->is_const()){
			str.append("const ");
		}
		const ClassNameNode* classNameNode = GetClassNameNode(L, index);
		if(classNameNode == NULL || classNameNode->empty()) str.append("userdata/");
		else str.append(classNameNode->name()); str.append("*/");
	}

	if(lua_type(L, index) == LUA_TNONE)	str.append("none/");

	if(str.size() == 0)
		return string("unknown type");
	else
		return str.substr(0, str.size()-1);
}

/// prints information about lua's call stack (file:line source).
void lua_stacktrace(lua_State* L)
{
    lua_Debug entry;
    for(int depth = 1; lua_getstack(L, depth, &entry); depth++)
	{
    	int status = lua_getinfo(L, "Sln", &entry);
    	if(!status || !entry.short_src || entry.currentline < 0) return;
		UG_LOG(entry.short_src << ":" << entry.currentline);
		UG_LOG(" " << GetFileLine(entry.short_src, entry.currentline));
		UG_LOG("\n");
    }
}

/// returns the current file and line ( \sa lua_stacktrace ).
std::string GetLuaFileAndLine(lua_State* L)
{
	lua_Debug entry;
	lua_getstack(L, 1, &entry);
	int status = lua_getinfo(L, "Sln", &entry);
	if(!status || !entry.short_src || entry.currentline < 0) return std::string("");
	std::stringstream ss;
	ss << entry.short_src << ":" << entry.currentline;
	return ss.str();
}

void ScriptStacktrace()
{
	lua_stacktrace(script::GetDefaultLuaState());
}

bool ScriptPrintClassHierarchy(const char *classname)
{
	return PrintClassHierarchy(GetUGRegistry(), classname);
}

bool RegisterInfoCommands(Registry &reg, const char* parentGroup)
{
	try
	{
		stringstream grpSS; grpSS << parentGroup << "/Info";
		std::string grp = grpSS.str();

		reg.add_function("ls", &LuaList, grp.c_str());
		reg.add_function("TypeInfo", &UGTypeInfo, grp.c_str());
		reg.add_function("ClassUsage", &ClassUsage, grp.c_str());
		reg.add_function("ClassInstantiations" ,&ClassInstantiations, grp.c_str());
		reg.add_function("ClassHierarchy" ,&ScriptPrintClassHierarchy, grp.c_str());
		reg.add_function("Stacktrace", &ScriptStacktrace, grp.c_str());
	}
	catch(UG_REGISTRY_ERROR_RegistrationFailed ex)
	{
		UG_LOG("### ERROR in RegisterInfoCommands: "
				"Registration failed (using name " << ex.name << ").\n");
		return false;
	}

	return true;
}


} // namespace bridge

} // namespace ug
