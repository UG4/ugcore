/*
 * lua_debug.h
 *
 *  Created on: 17.03.2012
 *      Author: mrupp
 */

#ifndef LUA_DEBUG_H_
#define LUA_DEBUG_H_

#include "registry/registry.h"

namespace ug
{
namespace script
{

/**
 * enum used to control execution flow in debug mode
 */
enum debug_return
{
	DEBUG_EXIT=0,  //!< exit ug
	DEBUG_CONTINUE,//!< continue execution
	DEBUG_NEXT,    //!< go to next line, but do not go deeper in stack
	DEBUG_STEP,    //!< go to next line, step into functions (deeper in stack)
	DEBUG_FINISH   //!< continue until we finish current function
};

/**
 * Register debug/profile functions
 * @param reg
 * @return true
 * @sa debug_return
 */
UG_API bool RegisterLuaDebug(ug::bridge::Registry &reg);

/**
 * function called when a breakpoint is reached
 * @param s debug shell function
 * @return 0
 * @sa DebugList, DebugBacktrace, DebugDown, DebugUp
 */
UG_API int SetDebugShell(debug_return (*s)());

UG_API void ProfileLUA(bool bProfile);

/// lists the current line in the script
UG_API void DebugList();

/// lists the function stack in lua
UG_API void DebugBacktrace(int fromLevel);

/// move down function stack
UG_API void DebugDown();

/// move down function stack
UG_API void DebugUp();

/// Free all memory associated with lua_debug
UG_API void FinalizeLUADebug();

UG_API void SetLuaDebugIDs(lua_State* L);

UG_API void DebugHold();
}
}

#endif /* LUA_DEBUG_H_ */
