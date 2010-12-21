// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y10 m01 d20

#include <cassert>
#include <cstring>
#include "ug.h"
#include "ug_script.h"
#include "bindings/bindings_lua.h"
#include "ug_bridge/ug_bridge.h"
#include "ug_bridge/class_helper.h"
#include "info_commands.h"
#include "user_data/user_data.h"
using namespace std;

namespace ug
{

namespace script
{

///	loads a file relative to APP_PATH../scripts
bool LoadUGScript(const char* filename)
{
	string strFilename(UGGetScriptPath());
	strFilename.append("/").append(filename);

//	check absolute path
	if(FileExists(filename)){
		return ug::script::ParseFile(filename);
	}

//	check script path
	if(FileExists(strFilename.c_str())){
		return ug::script::ParseFile(strFilename.c_str());
	}

//todo: check path relative to UG4_ROOT

	return false;
}

/// checks if given file exists.
bool FileExists(const char* filename)
{
//todo: this could be improved.
	ifstream in(filename);
	if(in) {
		in.close();
		return true;
	}
	return false;
}

static ug::bridge::Registry* g_pRegistry = NULL;

static void UpdateScriptAfterRegistryChange(ug::bridge::Registry* pReg)
{
	UG_ASSERT(pReg == g_pRegistry, "static g_pRegistry does not match parameter pReg, someone messed up the registries!");
	
//	this can be called since CreateBindings automatically avoids
//	double registration
	ug::bridge::lua::CreateBindings_LUA(GetDefaultLuaState(),
										*pReg);
}

lua_State* GetDefaultLuaState()
{
	static lua_State* L = NULL;
	
//	if the state has not already been opened then do it now.
	if(!L)
	{
		if(!g_pRegistry){
		//	store a pointer to the registry and avoid multiple callback registration
			g_pRegistry = &ug::bridge::GetUGRegistry();
			g_pRegistry->add_callback(UpdateScriptAfterRegistryChange);
		}
		
	//	open a lua state
		L = lua_open();
	//	open standard libs
		luaL_openlibs(L);
	//	create lua bindings for registered functions and objects
		ug::bridge::lua::CreateBindings_LUA(L, *g_pRegistry);
		
	//	we use an extra registry to register some lua-only commands
		static ug::bridge::Registry scriptRegistry;

	//	this define makes sure that no methods are referenced that
	//	use the algebra, even if no algebra is included.
		#ifdef UG_ALGEBRA
		//	Register info commands
			RegisterInfoCommands(scriptRegistry);

		//	Register user functions
			RegisterLuaUserData(scriptRegistry, "/ug4");

		//	Register Boundary functions
			RegisterLuaBoundaryNumber(scriptRegistry, "/ug4");

			ug::bridge::lua::CreateBindings_LUA(L, scriptRegistry);
		#endif
			scriptRegistry.check_consistency();
	}
	
	return L;
}

bool ParseBuffer(const char* buffer)
{
	lua_State* L = GetDefaultLuaState();
	int error = luaL_loadbuffer(L, buffer, strlen(buffer), "buffer") ||
				lua_pcall(L, 0, 0, 0);

	if(error)
	{
		LOG("PARSE-ERROR: " << lua_tostring(L, -1) << endl);
		lua_pop(L, 1);
		return false;
	}
	return true;

}

bool ParseFile(const char* filename)
{
	lua_State* L = GetDefaultLuaState();

	UG_LOG("Parsing file: " << filename << std::endl);

	int error = luaL_loadfile(L, filename) || lua_pcall(L, 0, 0, 0);

	if(error)
	{
		LOG("PARSE-ERROR in parse_file(" << filename << "): " << lua_tostring(L, -1) << endl);
		lua_pop(L, 1);
		return false;
	}
	return true;
}

}}//	end of namespace
