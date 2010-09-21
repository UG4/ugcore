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

namespace ug
{
namespace script
{

///	returns the default lua state
lua_State* GetDefaultLuaState();

///	parses and executes a buffer
bool ParseBuffer(const char* buffer);

///	parses and executes a file
bool ParseFile(const char* filename);

}//	end of namespace
}//	end of namespace


#endif
