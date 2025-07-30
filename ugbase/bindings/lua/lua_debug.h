/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
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
