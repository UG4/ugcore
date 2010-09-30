// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y10 m01 d20

#ifndef __H__UG__UG_SCRIPT__
#define __H__UG__UG_SCRIPT__
#include <vector>

extern "C" {
#include "externals/lua/lua.h"
#include "externals/lua/lauxlib.h"
#include "externals/lua/lualib.h"
}

#include "common/common.h"
#include "../ug_bridge/registry.h"

namespace ug
{
namespace script
{

///	sets the registry from which methods and classes will be used
/**	As long as the script uses this registry, it may not be deleted!.*/
void SetScriptRegistry(ug::interface::InterfaceRegistry* pReg);

///	returns the default lua state
lua_State* GetDefaultLuaState();

///	parses and executes a buffer
bool ParseBuffer(const char* buffer);

///	parses and executes a file
bool ParseFile(const char* filename);

}//	end of namespace
}//	end of namespace


#endif
