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

enum debug_return
{
	DEBUG_EXIT=0,
	DEBUG_CONTINUE,
	DEBUG_NEXT,
	DEBUG_STEP,
	DEBUG_FINISH
};

int SetDebugShell(debug_return (*s)());

bool RegisterLuaDebug(ug::bridge::Registry &reg);
void ProfileLUA(bool bProfile);
void DebugList();
void DebugBacktrace();
void DebugDown();
void DebugUp();

}
}


#endif /* LUA_DEBUG_H_ */
