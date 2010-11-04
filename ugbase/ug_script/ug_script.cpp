// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y10 m01 d20

#include <cassert>
#include <cstring>
#include "ug_script.h"
#include "bindings/bindings_lua.h"
#include "ug_bridge/ug_bridge.h"
#include "info_commands.h"

using namespace std;

namespace ug
{

namespace script
{

static ug::bridge::Registry* g_pRegistry = NULL;

static void UpdateScriptAfterRegistryChange(ug::bridge::Registry* pReg)
{
	assert((pReg == g_pRegistry) && "someone messed up the registries!");
	
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
		RegisterInfoCommands(scriptRegistry);
		ug::bridge::lua::CreateBindings_LUA(L, scriptRegistry);
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
