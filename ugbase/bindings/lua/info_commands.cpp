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
#ifdef UG_PLUGINS
	#ifndef UG_EMBEDDED_PLUGINS
		#include "common/util/plugin_util.h"
	#endif
#endif

extern "C"
{
#include "bindings/lua/externals/lua/lstate.h"
}

#include "info_commands.h"


using namespace std;

namespace ug
{
	extern bool useLua2C;

namespace bridge
{


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
	char buf[512];
	fstream file(filename, ios::in);
	if(file.is_open() == false) return string("");
	stringstream *pss = NULL;
	for(size_t i=0; i<fromline-1 && !file.eof(); i++)
	{
		file.getline(buf, 512);
		if(strncmp(buf+strspn(buf, "\t "), "--", 2)==0)
		{
			if(pss == NULL) pss = new stringstream;
			*pss << "   \t" << buf+strspn(buf, "-\t ") << '\n';
		}
		else if(pss) { delete pss; pss = NULL; }

	}
	stringstream ss;
	if(pss != NULL) ss << "\n" << pss->str() << "\n";
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
		lua_remove(L, -2);
		if(lua_isnil(L, -1))
			return false;
	}

	return true;
}


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
			const char *src = ar.source[0]=='@' ? ar.source+1 : ar.source;
			UG_LOG(src << ":" << ar.linedefined << "-" << ar.lastlinedefined << "\n");
			UG_LOG(GetFileLinesLUA(src, ar.linedefined, ar.lastlinedefined) << "\n");
		}
		else
		{
			UG_LOG(GetFileLine(ar.source+1, ar.linedefined) << "\n");
		}
	}

	return 0;
}

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
	UG_LOG("\n");
	const bridge::Registry &reg = GetUGRegistry();

	// check if it is a class
	const ClassGroupDesc *cg = reg.get_class_group(p);
	const IExportedClass *c=NULL;
	if(cg != NULL)
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
			PrintClassInfo(reg, names->at(i));
		UG_LOG(endl);
		PrintClassHierarchy(reg, c->name().c_str());
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
		if(!names->empty())
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
	const IExportedClass *c = reg.get_class(classname);
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
	const IExportedClass *c = reg.get_class(classname);
	if(c == NULL)
	{
		UG_LOG("Class name " << classname << " not found\n");
		return false;
	}

	// print usages in functions

	UG_LOG("--- Functions returning " << classname << ": ---\n");
	ClassUsageExact(reg, classname, true);

	const std::vector<const char*> *names = c->class_names();
	if(names != NULL && !names->empty())
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

void PrintLuaScriptFunction(lua_State *L, int index)
{
	LUA_STACK_CHECK(L, 0);
	lua_pushvalue(L, index);
	lua_Debug ar;
	lua_getinfo(L, ">S", &ar);
	if(ar.linedefined != -1)
	{

		const char *p=GetFileLine(ar.source[0] == '@' ? ar.source+1 : ar.source,
				ar.linedefined).c_str();
		p+=strspn(p, " \t");
		UG_LOG(p);
	}
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
		int index = sorted[i].index;

		UG_LOG(repeat(' ', iSpace));
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
			if(reg.get_class(luastr)==false)
			{
				if(FindFunction(reg, luastr))
				{
					if(functions) functions->push_back(luastr);
				}
				else if(lua_isfunction(L, -1) || lua_iscfunction(L, -1))
				{
					lua_Debug ar;
					lua_pushvalue(L, -1);
					lua_getinfo(L, ">S", &ar);
					if(ar.linedefined != -1) {
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
	GetLuaGlobals(&functions, NULL, NULL, NULL, NULL);
}

void GetLuaGlobal_internalFunctions(std::vector<std::string> &internalFunctions)
{
	GetLuaGlobals(NULL, &internalFunctions, NULL, NULL, NULL);
}

void GetLuaGlobal_scriptFunctions(std::vector<std::string> &scriptFunctions)
{
	GetLuaGlobals(NULL, NULL, &scriptFunctions, NULL, NULL);
}

void GetLuaGlobal_luaObjects(std::vector<std::string> &luaObjects)
{
	GetLuaGlobals(NULL, NULL, NULL, &luaObjects, NULL);
}

void GetLuaGlobal_classInstantiations(std::vector<std::string> &classInstantiations)
{
	GetLuaGlobals(NULL, NULL, NULL, NULL, &classInstantiations);
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
	{
		if(PrintFunctionInfo(reg, functions[i].c_str()) == false)
			UG_LOG(functions[i]);
		UG_LOG(endl);
	}
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
		PrintLuaScriptFunction(L);
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
		if(luaObjects[i].compare("package") == 0) continue;
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
		if(!lua_isuserdata(L, -1)) continue;
		const std::vector<const char*> *n  = GetClassNames(L, -1);
		if(n && !n->empty())
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
		lua_pushnil(L);
		str.append("{");
		bool bFirst = true;
		while (lua_next(L, index) != 0) {
			if(bFirst) {bFirst = false;} else {str.append(", ");};
			str.append(GetLuaTypeString(L, -1));
			lua_pop(L, 1);
	   }
		str.append("} (table)/");
	}
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



void LuaGetLastLine(lua_State* L, lua_Debug entry)
{
    for(int depth = 0; lua_getstack(L, depth, &entry); depth++)
	{
    	int status = lua_getinfo(L, "Sln", &entry);
    	if(!status || entry.currentline < 0) continue;
    }
}


void LuaPrintCurrentLine(lua_State* L)
{
	lua_Debug entry;
	LuaGetLastLine(L, entry);
	UG_LOG(entry.short_src << ":" << entry.currentline);
	UG_LOG(" " << GetFileLine(entry.short_src, entry.currentline));
	UG_LOG("\n");

}
/// prints information about lua's call stack (file:line source).
void LuaStackTrace(lua_State* L)
{
    lua_Debug entry;
    for(int depth = 0; lua_getstack(L, depth, &entry); depth++)
	{
    	int status = lua_getinfo(L, "Sln", &entry);
    	if(entry.currentline <0) continue;
    	if(!status || !entry.short_src || entry.currentline < 0) continue;
    	UG_LOG(entry.short_src << ":" << entry.currentline);
    	UG_LOG(" " << GetFileLine(entry.short_src, entry.currentline));
    	UG_LOG("\n");
    }
}


/// returns the current file and line ( \sa LuaStackTrace ).
bool GetLuaFileAndLine(lua_State* L, std::string &file, size_t &line)
{
	lua_Debug entry;
	lua_getstack(L, 1, &entry);
	int status = lua_getinfo(L, "Sln", &entry);
	if(!status || !entry.short_src || entry.currentline < 0) return false;
	file = entry.short_src;
	line = entry.currentline;
	return true;
}

/// returns the current file and line ( \sa LuaStackTrace ).
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
	LuaStackTrace(script::GetDefaultLuaState());
}

bool ScriptPrintClassHierarchy(const char *classname)
{
	return PrintClassHierarchy(GetUGRegistry(), classname);
}

	
bool ScriptHasClass(const char *classname)
{
	const Registry &reg = GetUGRegistry();
	return reg.get_class(classname) != NULL;
}

bool ScriptHasClassGroup(const char *classname)
{
	const Registry &reg = GetUGRegistry();
	return reg.get_class_group(classname) != NULL;
}

#ifdef UG_PLUGINS
#ifndef UG_EMBEDDED_PLUGINS
bool AssertPluginLoaded(const char *name)
{
	if(PluginLoaded(name) == false)
	{
		string msg = string("plugin ") + name + string(" not loaded. Please use 'cmake -D") + name + string("=ON ..' in your build directory.");		
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
#endif

void EnableLUA2C(bool b)
{
#ifndef USE_LUA2C
	UG_LOG("Warning: LUA2C not enabled. Enable with \"cmake -DUSE_LUA2C=ON ..\"\n")
#endif
	useLua2C=b;
}

bool RegisterInfoCommands(Registry &reg, const char* parentGroup)
{
	stringstream grpSS; grpSS << parentGroup << "/Info";
	std::string grp = grpSS.str();

	try
	{
		reg.add_function("ls", &LuaList, grp.c_str(), "", "", "list all objects");
		reg.add_function("list_cfunctions", &LuaList_cfunctions, grp.c_str(), "", "", "list all cfunctions");
		reg.add_function("list_classes", &LuaList_classes, grp.c_str(), "", "", "list all classes");
		reg.add_function("list_internalFunctions", &LuaList_internalFunctions, grp.c_str(), "", "", "list all of LUAs internal functions");
		reg.add_function("list_luaObjects", &LuaList_luaObjects, grp.c_str(), "", "", "list all created LUA objects");
		reg.add_function("list_scriptFunctions", LuaList_scriptFunctions, grp.c_str(), "", "", "list all LUA script functions");
		reg.add_function("TypeInfo", &UGTypeInfo, grp.c_str(), "", "\"typeName\"", "print information about a type");
		reg.add_function("ClassUsage", &ClassUsage, grp.c_str(), "", "\"typeName\"", "print information about the usage of a type");
		reg.add_function("ClassInstantiations" ,&ClassInstantiations, grp.c_str(), "", "\"typeName\"", "print all objects of the type");
		reg.add_function("ClassHierarchy" ,&ScriptPrintClassHierarchy, grp.c_str(), "", "\"typeName\"", "print the class hierachy of type");
		reg.add_function("Stacktrace", &ScriptStacktrace, grp.c_str(), "", "", "prints the LUA function stack, that is which functions are called up to this point");
		reg.add_function("HasClass", &ScriptHasClass, grp.c_str(), "true if class exists", "\"className\"", "use only if you know that you're not using a class group, otherwise HasClassGroup");
		reg.add_function("HasClassGroup", &ScriptHasClassGroup, grp.c_str(), "true if class oder classGroup exists", "\"classGroupName\"", "can be used before instantiating a class");
#ifdef UG_PLUGINS
	#ifndef UG_EMBEDDED_PLUGINS
		reg.add_function("PluginLoaded", &PluginLoaded, grp.c_str(), "true if plugin loaded", "\"pluginName\"", "pluginName as listed when using cmake ..");
		reg.add_function("AssertPluginLoaded", &AssertPluginLoaded, grp.c_str(), "true if plugin loaded", "\"pluginName\"", "throws an error if plugin not loaded, displays help string how to enable plugins via cmake -DpluginName=ON ..");
	#endif
#endif
		reg.add_function("EnableLUA2C", &EnableLUA2C, grp.c_str(), "", "bEnable", "");
	}
	UG_REGISTRY_CATCH_THROW(grp);

	return true;
}


} // namespace bridge

} // namespace ug
