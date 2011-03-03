
#ifndef __H__UG_BRIDGE__BRIDGES__USER_DATA__USER_DATA__
#define __H__UG_BRIDGE__BRIDGES__USER_DATA__USER_DATA__

#include "ug_script/ug_script.h"
#include <stdarg.h>

namespace ug
{
namespace bridge
{

/// this class maps a scalar value an output scalar value using a lua callback
class LuaUserNumberNumberFunction
{
	public:
		LuaUserNumberNumberFunction();

		void set_lua_callback(const char* luaCallback);

		number operator() ( int numArgs, ... ) const;

	protected:
		const char* m_callbackName;
		int m_callbackRef;
		lua_State*	m_L;
};

void RegisterLuaUserData(Registry& reg, const char* parentGroup);

void RegisterLuaBoundaryNumber(Registry& reg, const char* parentGroup);

} // end namepace bridge
} // end namespace ug

#endif /* __H__UG_BRIDGE__BRIDGES__USER_DATA__USER_DATA__ */
