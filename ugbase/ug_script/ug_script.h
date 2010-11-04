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

///	returns the default lua state
/**	When called for the first time, a new state is created and
 *	the methods and classes in ugs default registry
 *	(ug::bridge::GetUGRegistry) are registered. Furthermore a callback
 *	is registered, which registers new methods whenever
 *	Registry::registry_changed() is called on the default registry.*/
lua_State* GetDefaultLuaState();

///	parses and executes a buffer
bool ParseBuffer(const char* buffer);

///	parses and executes a file
bool ParseFile(const char* filename);

}//	end of namespace
}//	end of namespace


#endif
