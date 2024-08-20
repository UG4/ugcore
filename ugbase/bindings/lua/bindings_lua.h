/*
 * Copyright (c) 2010-2014:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
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

#ifndef __H__UG__INTERFACE__BINDINGS_LUA__
#define __H__UG__INTERFACE__BINDINGS_LUA__

#include <vector>
#include <string>
#include "common/types.h"
extern "C" {
#include "externals/lua/lua.h"
#include "externals/lua/lauxlib.h"
#include "externals/lua/lualib.h"
}


#include "common/common.h"
#include "registry/registry.h"

namespace ug
{
namespace bridge
{
namespace lua
{

extern const bool IMLPICIT_SMART_PTR_TO_PTR_CONVERSION;

enum UserDataWrapperTypes{
	RAW_POINTER = 1,
	SMART_POINTER = 1 << 1,
	IS_CONST = 1 << 2
};

struct UserDataWrapper{
	byte_t type;

	bool is_const()		{return (type & IS_CONST) == IS_CONST;}
	bool is_raw_ptr()	{return (type & RAW_POINTER) == RAW_POINTER;}
	bool is_smart_ptr()	{return (type & SMART_POINTER) == SMART_POINTER;}
};

struct SmartUserDataWrapper : public UserDataWrapper
{
	SmartPtr<void>	smartPtr;
};

struct ConstSmartUserDataWrapper : public UserDataWrapper
{
	ConstSmartPtr<void>	smartPtr;
};

struct RawUserDataWrapper : public UserDataWrapper
{
	void*	obj;
	void (*deleteFunc)(const void*);
};


///	creates bindings for ug_interface and a given lua-state.
/**	If you use ug::script, this method will be invoked automatically.*/
bool CreateBindings_LUA(lua_State* L, Registry& reg);

SmartUserDataWrapper* CreateNewUserData(lua_State* L, const SmartPtr<void>& ptr,
										const char* metatableName);

ConstSmartUserDataWrapper* CreateNewUserData(lua_State* L, const ConstSmartPtr<void>& ptr,
											 const char* metatableName);

RawUserDataWrapper* CreateNewUserData(lua_State* L, void* ptr,
									  const char* metatableName,
									  void (*deleteFunc)(const void*),
									  bool is_const);


}//	end of namespace
}//	end of namespace
}//	end of namespace

#endif
