/*
 * Copyright (c) 2014-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef LUA_PARSING_H_
#define LUA_PARSING_H_

#ifdef USE_LUAJIT
#include <lua.h>
#else

extern "C" {
#include "externals/lua/lua.h"
#include "bindings/lua/externals/lua/lauxlib.h"
}

#endif

#include "bindings_lua.h"
#include "bindings/lua/lua_function_handle.h"
#include "bindings/lua/lua_table_handle.h"

namespace ug {
namespace bridge {
namespace lua {

template <typename T>
struct LuaParsing;

template <>
struct LuaParsing<bool>{
	static bool check(lua_State* L, int index){
		return lua_isboolean(L, index);
	}
	static bool get(lua_State* L, int index){
		return lua_toboolean(L, index);
	}
	static void push(lua_State* L, bool data){
		lua_pushboolean(L, (data ? 1 : 0));
	}
};

template <>
struct LuaParsing<int>{
	static bool check(lua_State* L, int index){
		return lua_isnumber(L, index);
	}
	static int get(lua_State* L, int index){
		return (int)lua_tointeger(L, index);
	}
	static void push(lua_State* L, int data){
		lua_pushnumber(L, data);
	}
};

template <>
struct LuaParsing<size_t>{
	static bool check(lua_State* L, int index){
		return lua_isnumber(L, index);
	}
	static size_t get(lua_State* L, int index){
		return lua_tointeger(L, index);
	}
	static void push(lua_State* L, size_t data){
		lua_pushnumber(L, (lua_Number)data);
	}
};

template <>
struct LuaParsing<float>{
	static bool check(lua_State* L, int index){
		return lua_isnumber(L, index);
	}
	static float get(lua_State* L, int index){
		return (float)lua_tonumber(L, index);
	}
	static void push(lua_State* L, float data){
		lua_pushnumber(L, data);
	}
};

template <>
struct LuaParsing<double>{
	static bool check(lua_State* L, int index){
		return lua_isnumber(L, index);
	}
	static double get(lua_State* L, int index){
		return lua_tonumber(L, index);
	}
	static void push(lua_State* L, double data){
		lua_pushnumber(L, data);
	}
};

template <>
struct LuaParsing<const char*>{
	static bool check(lua_State* L, int index){
		return lua_isstring(L, index);
	}
	static const char* get(lua_State* L, int index){
		return lua_tostring(L, index);
	}
	static void push(lua_State* L, const char* data){
		lua_pushstring(L, data);
	}
};

template <>
struct LuaParsing<std::string>{
	static bool check(lua_State* L, int index){
		return lua_isstring(L, index);
	}
	static std::string get(lua_State* L, int index){
		return std::string(lua_tostring(L, index));
	}
	static void push(lua_State* L, std::string data){
		lua_pushstring(L, data.c_str());
	}
};

template <>
struct LuaParsing<const std::string&>{
	static bool check(lua_State* L, int index){
		return lua_isstring(L, index);
	}
	static void push(lua_State* L, const std::string& data){
		lua_pushstring(L, data.c_str());
	}
};

#ifdef UG_FOR_LUA
template <>
struct LuaParsing<LuaFunctionHandle>{
	static bool check(lua_State* L, int index){
		return lua_isfunction(L, index);
	}
	static LuaFunctionHandle get(lua_State* L, int index){
		LuaFunctionHandle tmp;
		lua_pushvalue(L, index);
		tmp.ref = luaL_ref(L, LUA_REGISTRYINDEX);
		return tmp;
	}
	static void push(lua_State* L, LuaFunctionHandle data){
		UG_THROW("Return value of type LuaFunctionHandle not implemented.");
	}
};

template <>
struct LuaParsing<LuaTableHandle>{
	static bool check(lua_State* L, int index){
		return lua_istable(L, index);
	}
	static LuaTableHandle get(lua_State* L, int index){
		LuaTableHandle tmp(L, index);
		untested();
		return tmp;
	}
	static void push(lua_State* L, LuaTableHandle data){
		UG_THROW("Return value of type LuaTableHandle not implemented.");
	}
};
#endif

template <>
struct LuaParsing<void*>{
	static bool checkAndGet(std::pair<void*, const ClassNameNode*>& res,
	                        lua_State* L, int index, const char* baseClassName){
		if(!lua_isuserdata(L, index)) return false;

		UserDataWrapper* udata =
			reinterpret_cast<UserDataWrapper*>(lua_touserdata(L, index));

		if(udata->is_const()) return false;

		//	extract the pointer to the object.
		//	udata is either a RawUserData or a SmartUserDataWrapper
		void* obj = nullptr;
		if(udata->is_raw_ptr())
			obj = static_cast<RawUserDataWrapper*>(udata)->obj;
		else if(udata->is_smart_ptr() && IMLPICIT_SMART_PTR_TO_PTR_CONVERSION)
			obj = static_cast<SmartUserDataWrapper*>(udata)->smartPtr.get();
		else return false;

		//	get the object and its metatable. Make sure that obj can be cast to
		//	the type that is expected by the paramsTemplate.
		if(lua_getmetatable(L, index) == 0) return false;

		lua_pushstring(L, "class_name_node");
		lua_rawget(L, -2);
		const ClassNameNode* classNameNode = (const ClassNameNode*) lua_touserdata(L, -1);
		lua_pop(L, 2);

		if(!classNameNode) return false;
		if(classNameNode->empty()) return false;
		if(!ClassNameTreeContains(*classNameNode, baseClassName)) return false;

		res.first = obj;
		res.second = classNameNode;

		return true;
	}

	static void push(lua_State* L, void* data, const char* className){
		CreateNewUserData(L, data, className, NULL, false);
	}
};

template <>
struct LuaParsing<const void*>{
	static bool checkAndGet(std::pair<const void*, const ClassNameNode*>& res,
	                        lua_State* L, int index, const char* baseClassName){
		if(!lua_isuserdata(L, index)) return false;

		auto* udata = reinterpret_cast<UserDataWrapper*>(lua_touserdata(L, index));

		//	extract the pointer to the object.
		//	udata is either a RawUserData or a SmartUserDataWrapper
		const void* obj = nullptr;

		if(udata->is_raw_ptr())
			obj = static_cast<RawUserDataWrapper*>(udata)->obj;
		else if(udata->is_smart_ptr() && IMLPICIT_SMART_PTR_TO_PTR_CONVERSION){
			//	we have to distinguish between const and non-const smart pointers.
			if(udata->is_const())
				obj = static_cast<ConstSmartUserDataWrapper*>(udata)->smartPtr.get();
			else
				obj = static_cast<SmartUserDataWrapper*>(udata)->smartPtr.get();
		}
		else return false;

		//	get the object and its metatable. Make sure that obj can be cast to
		//	the type that is expected by the paramsTemplate.
		if(lua_getmetatable(L, index) == 0) return false;

		lua_pushstring(L, "class_name_node");
		lua_rawget(L, -2);
		auto* classNameNode = static_cast<const ClassNameNode *>(lua_touserdata(L, -1));
		lua_pop(L, 2);

		if(!classNameNode) return false;
		if(classNameNode->empty()) return false;
		if(!ClassNameTreeContains(*classNameNode, baseClassName)) return false;

		res.first = obj;
		res.second = classNameNode;

		return true;
	}

	static void push(lua_State* L, const void* data, const char* className){
	//	we're removing const with a cast. However, it was made sure that
	//	obj is treated as a const value.
		CreateNewUserData(L, const_cast<void *>(data), className, nullptr, true);
	}
};

template <>
struct LuaParsing<SmartPtr<void> >{
	static bool checkAndGet(std::pair<SmartPtr<void>, const ClassNameNode*>& res,
	                        lua_State* L, int index, const char* baseClassName){
		if(!lua_isuserdata(L, index)) return false;

		auto* udata = static_cast<UserDataWrapper*>(lua_touserdata(L, index));

		if(!udata->is_smart_ptr()) return false;
		if(udata->is_const()) return false;

		SmartPtr<void>& obj = static_cast<SmartUserDataWrapper *>(lua_touserdata(L, index))->smartPtr;

		if(lua_getmetatable(L, index) == 0) return false;
		lua_pushstring(L, "class_name_node");
		lua_rawget(L, -2);
		const auto* classNameNode = static_cast<const ClassNameNode *>(lua_touserdata(L, -1));
		lua_pop(L, 2);

		if(!classNameNode) return false;
		if(classNameNode->empty()) return false;
		if(!ClassNameTreeContains(*classNameNode, baseClassName)) return false;

		res.first = obj;
		res.second = classNameNode;

		return true;
	}

	static void push(lua_State* L, SmartPtr<void> data, const char* className){
		CreateNewUserData(L, data, className);
	}
};

template <>
struct LuaParsing<ConstSmartPtr<void> >{
	static bool checkAndGet(std::pair<ConstSmartPtr<void>, const ClassNameNode*>& res,
	                        lua_State* L, int index, const char* baseClassName){
		if(!lua_isuserdata(L, index)) return false;

		auto* udata = static_cast<UserDataWrapper*>(lua_touserdata(L, index));

		if(!udata->is_smart_ptr()) return false;

		ConstSmartPtr<void> obj;
		if(static_cast<UserDataWrapper *>(lua_touserdata(L, index))->is_const())
			obj = static_cast<ConstSmartUserDataWrapper *>(lua_touserdata(L, index))->smartPtr;
		else
			obj = static_cast<SmartUserDataWrapper *>(lua_touserdata(L, index))->smartPtr;

		if(lua_getmetatable(L, index) == 0) return false;
		lua_pushstring(L, "class_name_node");
		lua_rawget(L, -2);
		auto* classNameNode = static_cast<const ClassNameNode *>(lua_touserdata(L, -1));
		lua_pop(L, 2);

		if(!classNameNode) return false;
		if(classNameNode->empty()) return false;
		if(!ClassNameTreeContains(*classNameNode, baseClassName)) return false;

		res.first = obj;
		res.second = classNameNode;

		return true;
	}

	static void push(lua_State* L, ConstSmartPtr<void> data, const char* className){
		CreateNewUserData(L, data, className);
	}
};


}
}
}
#endif
