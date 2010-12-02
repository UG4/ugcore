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
#include "ug_script/ug_script.h"
#include "ug_bridge/ug_bridge.h"
#include "ug_bridge/class_helper.h"

extern "C"
{
#include "ug_script/externals/lua/lstate.h"
}

#include "info_commands.h"


using namespace std;

namespace ug
{
namespace bridge
{

namespace lua
{
	string GetLuaTypeString(lua_State* L, int index);
}



const std::vector<const char*> *GetClassNames(lua_State *L, int index)
{
	if(lua_getmetatable(L, index) != 0)
	{
		// get names
		lua_pushstring(L, "names");
		lua_rawget(L, -2);
		if(!lua_isnil(L, -1) && lua_isuserdata(L, -1))
		{
			return (const std::vector<const char*>*) lua_touserdata(L, -1);
		}
		lua_pop(L, 2); // pop userdata, metatable
	}
	return NULL;
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
bool UGTypeInfo(const char *p)
{
	bridge::Registry &reg = GetUGRegistry();

	UG_LOG("\n");
	const IExportedClass *c = FindClass(reg, p);
	if(c)
	{
		const std::vector<const char*> *names = c->class_names();
		for(size_t i=0; i < names->size(); ++i)
			PrintClassInfo(reg, names->at(i));
		UG_LOG(endl);
		PrintClassHierarchy(reg, c->name());
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
		UG_LOG(p << " is a function\n ");
		PrintFunctionInfo(reg, p);
		UG_LOG(endl);
	}
	else if(lua_iscfunction(L, -1))
	{
		UG_LOG(p << " is a cfunction\n ");
		PrintFunctionInfo(reg, p);
		UG_LOG(endl);
	}
	else if(lua_isuserdata(L, -1))
	{
		// names = GetClassNames(L, -1))
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
			PrintClassInfo(reg, names->at(i));
		if(names->size() > 0)
			PrintClassHierarchy(reg, names->at(0));
	}
	else
		UG_LOG(p << " is a " << lua::GetLuaTypeString(L, -1) << endl);

	lua_pop(L, 1); // pop globals

	UG_LOG("\n");
	return true;
}



bool ClassNameVecContains(const std::vector<const char*>& names, const char* name);
bool ClassInstantiations(const char *classname)
{
	bridge::Registry &reg = GetUGRegistry();
	const IExportedClass *c = FindClass(reg, classname);
	if(c == NULL)
	{
		UG_LOG("Class " << classname << " not found\n");
		return false;
	}

	UG_LOG(endl);
	UG_LOG("Instantiations of Class " << classname << ":\n");

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

/**
 *
 * \param classname the class to print usage in functions/member functions of (and all its subclasses) .
 * class in in/out parameters is highlighted with [class].
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


void LuaPrintTable(lua_State *L)
{
	UG_LOG("{ ");
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

	for(int i=0; i<len; i++)
	{
		if(i>0) UG_LOG(", ");
		const char * key = lua_tostring(L, -2*len+2*i);
		const char * value = lua_tostring(L, -2*len+2*i+1);
		UG_LOG(key);
		if(value) { UG_LOG(" = \"" << value << "\"") };
	}
	lua_pop(L, 2*len);
	UG_LOG(" }");
}



bool IsLonger(const std::string &a, const std::string &b)
{
	return b.size() > a.size();
}

void PrintFileLineFunction(const char *source, int linedefined)
{
	char buf[512];
	fstream file(source+1, ios::in);
	for(int i=0; i<linedefined; i++)
		file.getline(buf, 512);
	char *s = buf+strspn(buf, " \t");
	UG_LOG(s << " \t" << source << ":" << linedefined);

}

void LuaList()
{
	bridge::Registry &reg = GetUGRegistry();
	lua_State* L = script::GetDefaultLuaState();
	std::vector<std::string> classes, functions, nonregistered, names, instantiations;
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
			if(lua_isnil(L, -1))
			{
				lua_pop(L, 1); // remove global from stack
				continue; // nope
			}
			if(strcmp(luastr, "_G") == 0)
				;
			else if(FindClass(reg, luastr))
				classes.push_back(luastr);
			else if(FindFunction(reg, luastr))
				functions.push_back(luastr);
			else if(lua_isfunction(L, -1) || lua_iscfunction(L, -1))
				nonregistered.push_back(luastr);
			else if(lua_isuserdata(L, -1))
				instantiations.push_back(luastr);
			else
				names.push_back(luastr);

			lua_pop(L, 1); // remove global from stack
		}
	}
	sort(classes.begin(), classes.end());
	sort(functions.begin(), functions.end());
	sort(names.begin(), names.end());
	sort(nonregistered.begin(), nonregistered.end());
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

	UG_LOG(endl << "--- Not registered Functions: ---" << endl)

	std::vector<std::string>::const_iterator m = max_element(nonregistered.begin(), nonregistered.end(), IsLonger);
	for(size_t i=0; i<nonregistered.size(); i++)
	{
		lua_Debug ar;
		lua_getglobal(L, nonregistered[i].c_str());  /* get global 'f' */
		lua_getinfo(L, ">S", &ar);
		if(ar.linedefined != -1)
		{
			UG_LOG(left << setw((*m).size()) << nonregistered[i] << ": ");
			PrintFileLineFunction(ar.source, ar.linedefined);
			UG_LOG(endl);
		}
		else
		{ 	UG_LOG(nonregistered[i] << endl); }
		lua_pop(L, 1);
	}


	UG_LOG(endl << "--- Lua Objects: ----------------" << endl)
	m = max_element(names.begin(), names.end(), IsLonger);
	for(size_t i=0; i<names.size(); i++)
	{
		lua_getglobal(L, names[i].c_str());
		UG_LOG(left << setw((*m).size()) << names[i]);
		if(lua_istable(L, -1))
		{
			UG_LOG(" = ");
			LuaPrintTable(L);
		}
		const char *p = lua_tostring(L, -1);
		if(p) UG_LOG(": " << p)
		lua_pop(L, 1);
		UG_LOG(endl);
	}

	UG_LOG(endl << "--- Class Instantiations: ---------" << endl)
	m = max_element(instantiations.begin(), instantiations.end(), IsLonger);
	for(size_t i=0; i<instantiations.size(); i++)
	{
		lua_getglobal(L, instantiations[i].c_str());
		if(!lua_isuserdata(L, -1)) continue;
		const std::vector<const char*> *n  = GetClassNames(L, -1);
		if(n && n->size() > 0)
		{
			UG_LOG(left << setw((*m).size()) << instantiations[i] << ": class ");
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
