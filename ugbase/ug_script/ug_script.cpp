// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y10 m01 d20

#include <cassert>
#include "ug_script.h"
#include "bindings/bindings_lua.h"

using namespace std;

namespace ug{
namespace script
{

lua_State* GetDefaultLuaState()
{
	static lua_State* L = NULL;
	
//	if the state has not already been opened then do it now.
	if(!L)
	{
	//	open a lua state
		L = lua_open();

	//	open standard libs
		luaL_openlibs(L);
		
	//	create lua bindings for registered functions and objects
		ug::interface::lua::CreateBindings_LUA(L);
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
