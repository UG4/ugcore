// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y10 m01 d20

//extern libraries
#include <cassert>
#include <cstring>
#include <string>
#include <stack>


// ug libraries
#include "ug.h"
#include "lua_util.h"
#include "common/util/path_provider.h"
#include "common/os_dependent/file_util.h"
#include "bindings_lua.h"
#include "bridge/bridge.h"
#include "registry/class_helper.h"
#include "info_commands.h"
#include "lua_user_data.h"
#include "registry/registry.h"
#include "lua_debug.h"

#ifdef UG_PARALLEL
#include "pcl/pcl.h"
#endif


using namespace std;

namespace ug
{

namespace script
{

stack<string> stkPathes;

///	loads a file relative to APP_PATH../scripts
bool LoadUGScript(const char* filename)
{
//	This stack is used to support relative pathes to a current script.
//	We do not use PathProvider here, since the current path not necessarily
//	holds the most current script-path.

//	get the path component of the current filename
	string strInFile = filename;
	string strInPath;
	string::size_type pos = strInFile.rfind("/");
	if(pos != string::npos){
		strInPath = strInFile.substr(0, pos);
	}
	else if((pos = strInFile.rfind("\\")) != string::npos){
		strInPath = strInFile.substr(0, pos);
	}

//	first we check whether the file can be found relative
//	to the location of the current script
	if(!stkPathes.empty()){
		string curPath = stkPathes.top();
		string file = curPath;
		file.append("/").append(filename);
		if(FileExists(file.c_str())){
			curPath.append("/").append(strInPath);
			stkPathes.push(curPath);
			PathProvider::push_current_path(curPath);
			bool success = ug::script::ParseFile(file.c_str());
			PathProvider::pop_current_path();
			stkPathes.pop();
			return success;
		}
	}

//	the file was not found relative to the current script (if there
//	was one at all). Check the default filename
	if(FileExists(filename)){
		stkPathes.push(strInPath);
		PathProvider::push_current_path(strInPath);
		bool success = ug::script::ParseFile(filename);
		PathProvider::pop_current_path();
		stkPathes.pop();
		return success;
	}

//	now we have to check relative to the script path
	string curPath = PathProvider::get_path(SCRIPT_PATH);
	string file = curPath;
	file.append("/").append(filename);

//	check script path
	if(FileExists(file.c_str())){
		curPath.append("/").append(strInPath);
		stkPathes.push(curPath);
		PathProvider::push_current_path(curPath);
		bool success = ug::script::ParseFile(file.c_str());
		PathProvider::pop_current_path();
		stkPathes.pop();
		return success;
	}


//	now we check relative to the apps path
	curPath = PathProvider::get_path(APPS_PATH);
	file = curPath;
	file.append("/").append(filename);

//	check script path
	if(FileExists(file.c_str())){
		curPath.append("/").append(strInPath);
		stkPathes.push(curPath);
		PathProvider::push_current_path(curPath);
		bool success = ug::script::ParseFile(file.c_str());
		PathProvider::pop_current_path();
		stkPathes.pop();
		return success;
	}

//todo: check path relative to UG4_ROOT
	UG_LOG("Couldn't find script " << filename << endl);
	return false;
}


static ug::bridge::Registry* g_pRegistry = NULL;

static void UpdateScriptAfterRegistryChange(ug::bridge::Registry* pReg)
{
	PROFILE_FUNC();
	UG_ASSERT(pReg == g_pRegistry, "static g_pRegistry does not match parameter pReg, someone messed up the registries!");
	
//	this can be called since CreateBindings automatically avoids
//	double registration
	ug::bridge::lua::CreateBindings_LUA(GetDefaultLuaState(),
										*pReg);
}

static lua_State* theLuaState = NULL;
lua_State* GetDefaultLuaState()
{
//	if the state has not already been opened then do it now.
	if(!theLuaState)
	{
		PROFILE_BEGIN(CreateLUARegistry);
		if(!g_pRegistry){
		//	store a pointer to the registry and avoid multiple callback registration
			g_pRegistry = &ug::bridge::GetUGRegistry();
			g_pRegistry->add_callback(UpdateScriptAfterRegistryChange);

			g_pRegistry->add_function("ug_load_script", &LoadUGScript, "/ug4/lua",
						"success", "", "Loads and parses a script and returns whether it succeeded.");

			RegisterLuaDebug(*g_pRegistry);

		//	this define makes sure that no methods are referenced that
		//	use the algebra, even if no algebra is included.
			#ifdef UG_ALGEBRA
		//	Register info commands
			RegisterInfoCommands(*g_pRegistry);

		//	Register user functions
			RegisterLuaUserData(*g_pRegistry, "/ug4");

			#endif

			if(!g_pRegistry->check_consistency())
				throw(UGFatalError("Script-Registry not ok."));

		}
		
	//	open a lua state
		theLuaState = lua_open();
	//	open standard libs
		luaL_openlibs(theLuaState);

	//	create lua bindings for registered functions and objects
		ug::bridge::lua::CreateBindings_LUA(theLuaState, *g_pRegistry);
		PROFILE_END();
	}
	
	return theLuaState;
}


void ReleaseDefaultLuaState()
{
	if(theLuaState != NULL)
	{
		lua_close(theLuaState);
		theLuaState = NULL;
	}
	FinalizeLUADebug();
	return;
}

int luaCallStackError( lua_State *L )
{
	UG_LOG("Error: " << lua_tostring(L, -1) << ". ");
    UG_LOG("call stack:\n"); ug::bridge::lua_stacktrace(L);
    return 1;
}

bool ParseBuffer(const char* buffer, const char *bufferName)
{
	lua_State* L = GetDefaultLuaState();
	lua_pushcfunction(L, luaCallStackError);
	int error = luaL_loadbuffer(L, buffer, strlen(buffer), bufferName) ||
				lua_pcall(L, 0, 0, -2);

	if(error)
	{
//		LOG("PARSE-ERROR: " << lua_tostring(L, -1) << endl);
		string msg = lua_tostring(L, -1);
		lua_pop(L, 1);
		ug::bridge::lua_stacktrace(L);
		throw(LuaError(msg.c_str()));
//		return false;
	}
	return true;

}

bool ParseFile(const char* filename)
{
	lua_State* L = GetDefaultLuaState();

	lua_pushcfunction(L, luaCallStackError);

	PROFILE_BEGIN(luaL_loadfile);
	int error = luaL_loadfile(L, filename);
	PROFILE_END();

	if(error == 0)
	{
		PROFILE_BEGIN(lua_pcall);
		error = lua_pcall(L, 0, 0, -2);
		PROFILE_END();
	}

	if(error)
	{
		//LOG("PARSE-ERROR in parse_file(" << filename << "): " << lua_tostring(L, -1) << endl);
		string msg = lua_tostring(L, -1);
		lua_pop(L, 1);
		ug::bridge::lua_stacktrace(L);
		throw(LuaError(msg.c_str()));
	}
	return true;
}

/// UGLuaPrint. Redirects LUA prints to UG_LOG
int UGLuaPrint(lua_State *L)
{
	int nArgs = lua_gettop(L);
	int i;
	lua_getglobal(L, "tostring");

	for(i=1; i<=nArgs; i++)
	{
		lua_pushvalue(L, -1);
		lua_pushvalue(L, i);
		lua_call(L, 1, 1);
		const char *s = lua_tostring(L, -1);
		if(s) UG_LOG(s);
		lua_pop(L, 1);
	}
	UG_LOG(endl);
	lua_pop(L,1);
	return 0;
}

/// UGLuaWrite. Redirects LUA prints to UG_LOG without adding newline at the end
int UGLuaWrite(lua_State *L)
{
	PROFILE_FUNC();
	int nArgs = lua_gettop(L);
	int i;
	lua_getglobal(L, "tostring");

	for(i=1; i<=nArgs; i++)
	{
		lua_pushvalue(L, -1);
		lua_pushvalue(L, i);
		lua_call(L, 1, 1);
		const char *s = lua_tostring(L, -1);
		if(s) UG_LOG(s);
		lua_pop(L, 1);
	}
	lua_pop(L,1);
	return 0;
}

}}//	end of namespace
