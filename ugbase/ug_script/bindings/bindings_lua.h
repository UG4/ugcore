//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m09 d16

#ifndef __H__UG__INTERFACE__BINDINGS_LUA__
#define __H__UG__INTERFACE__BINDINGS_LUA__

#include <vector>

extern "C" {
#include "../externals/lua/lua.h"
#include "../externals/lua/lauxlib.h"
#include "../externals/lua/lualib.h"
}

#include "common/common.h"

namespace ug
{
namespace interface
{
namespace lua
{

///	creates bindings for ug_interface and a given lua-state.
/**	If you use ug::script, this method will be invoked automatically.*/
bool CreateBindings_LUA(lua_State* L);

}//	end of namespace
}//	end of namespace
}//	end of namespace

#endif
