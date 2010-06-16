// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y10 m01 d20

#ifndef __H__UG__UG_SCRIPT__
#define __H__UG__UG_SCRIPT__

#include <vector>
#include <string>
#include <utility>

extern "C" {
#include "externals/lua/lua.h"
#include "externals/lua/lauxlib.h"
#include "externals/lua/lualib.h"
}

#include "common/common.h"

namespace ug{
namespace script
{

////////////////////////////////////////////////////////////////////////
//	initialization routines of the various script-modules.
//	make sure to call them between Initialize and Finalize.
////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////
//	lib-grid script
///	registers lib_grids scripting routines at the given lua-instance.
bool InitLibGridScript(lua_State* L);
void FinalizeLibGridScript(lua_State* L);

//////////////////////////////////////////////////////
//	lib-disc script
///	registers lib_discretizations scripting routines at the given lua-instance.
bool InitLibDiscScript(lua_State* L);
void FinalizeLibDiscScript(lua_State* L);


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	- methods to initialize and finalize lua.
//	- access to the main lua-instance
//	- helpers to parse files and buffers.
////////////////////////////////////////////////////////////////////////

///	initializes lua
void Initialize();

///	exits lua and calls finalize on the libraries
void Finalize();

///	returns the lua-state that is used by ugscript
lua_State* GetLuaState();

///	returns true after a script called exit.
bool ExitTriggered();

///	parses and executes a buffer
bool ParseBuffer(const char* buffer);

///	parses and executes a file
bool ParseFile(const char* filename);

}}// end of namespace

#endif
