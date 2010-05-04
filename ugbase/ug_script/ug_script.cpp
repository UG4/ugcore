// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y10 m01 d20

#include <cassert>

extern "C"{
#include "externals/lua/lua.h"
#include "externals/lua/lauxlib.h"
#include "externals/lua/lualib.h"
}

#include "luabind/luabind.hpp"
#include "common/common.h"

#include "ug_script.h"

//#include "lib_grid_script.h"
//#include "lib_grid/lib_grid.h"
//#include "lib_disc_script.h"

using namespace std;

#define ASSERT_LUA_IS_INITIALIZED assert(L && "ERROR in ugscript - ugscipt has not yet been initialized.");
namespace ug{
namespace script
{

static lua_State* L = NULL;
static bool bExit = false;


////////////////////////////////////////////////////////////////////////
//	some script-methods
void uglog(const char* str)
{
	UG_LOG(str << endl);
}

////////////////////////////////////////////////////////////////////////
///	initialises lib-grids script methods
void Initialize()
{
//	if the state has not already been opened then do it now.
	if(!L)
	{
	//	open a lua state
		L = lua_open();
	//	connect to luaBind
		luabind::open(L);
	//	open standard libs
		luaL_openlibs(L);

	//	register ugscript-base-methods
		luabind::module(L)[luabind::def("uglog", uglog)];
	}
}

void Finalize()
{
	ASSERT_LUA_IS_INITIALIZED;
	lua_close(L);
}


bool ExitTriggered()
{
	return bExit;
}

///	returns the lua-state that is used by ugscript
lua_State* GetLuaState()
{
	ASSERT_LUA_IS_INITIALIZED;
	return L;
}

///	parses and executes a buffer
bool ParseBuffer(const char* buffer)
{
	ASSERT_LUA_IS_INITIALIZED;
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

///	parses and executes a file
bool ParseFile(const char* filename)
{
	ASSERT_LUA_IS_INITIALIZED;
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
