/*
 * Copyright (c) 2015:  G-CSC, Goethe University Frankfurt
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

#ifndef LUA_STACK_H_
#define LUA_STACK_H_

#include "registry/registry.h"
#include "lua_parsing.h"
// #include "bindings_lua.h"
#include "common/util/smart_pointer.h"

namespace ug {
namespace bridge {

///	copies parameter values from the lua-stack to a parameter-list.
/**	\returns	The index of the first bad parameter starting from 1.
 *				Returns 0 if everything went right.
 *				Returns -1 if the number of parameters did not match.
 */
int LuaStackToParams(ParameterStack& ps,
							const ParameterInfo& psInfo,
							lua_State* L,
							int offsetToFirstParam = 0);

///	Pushes the parameter-values to the Lua-Stack.
/**
 * \returns The number of items pushed to the stack.
 */
int ParamsToLuaStack(const ParameterStack& ps, lua_State* L);

template <typename T>
static void ParamStackEntryToLuaStack(const ParameterStack& ps, lua_State* L,
                                      int index, bool bIsVector)
{
	if(!bIsVector){
		lua::LuaParsing<T>::push(L, ps.to<T>(index));
	}
	else {
		const std::vector<T>& vec = ps.to<std::vector<T> >(index);
		lua_createtable(L, vec.size(), 0);
		int newTable = lua_gettop(L);
		for(int i=0; i < (int)vec.size(); i++) {
			lua::LuaParsing<T>::push(L, vec[i]);
			lua_rawseti(L, newTable, i + 1);
		}
	}
}

template <typename T>
static void ParamStackPointerEntryToLuaStack(const ParameterStack& ps, lua_State* L,
                                             int index, bool bIsVector)
{
	const char* className = ps.class_name(index);
	if(!bIsVector){
		lua::LuaParsing<T>::push(L, ps.to<T>(index), className);
	}
	else {
		SmartPtr<std::vector<std::pair<T, const ClassNameNode*> > > spVec = ps.to<SmartPtr<std::vector<std::pair<T, const ClassNameNode*> > > >(index);
		lua_createtable(L, spVec->size(), 0);
		int newTable = lua_gettop(L);
		for(int i=0; i < (int)spVec->size(); i++) {
			lua::LuaParsing<T>::push(L, (*spVec)[i].first, (*spVec)[i].second->name().c_str());
			lua_rawseti(L, newTable, i + 1);
		}
	}
}

template <typename T>
static bool PushLuaStackEntryToParamStack(ParameterStack& ps, lua_State* L,
                                          int index, bool bIsVector)
{
	if(!bIsVector){
		if(lua::LuaParsing<T>::check(L, index)){
			ps.push(lua::LuaParsing<T>::get(L, index));
		}
		else return false;
	}
	else {
		if (lua_istable(L, index)){
			SmartPtr<std::vector<T> > spVec
							= SmartPtr<std::vector<T> >(new std::vector<T>());
			lua_pushnil(L);
			while (lua_next(L, index) != 0) {
				if(!lua::LuaParsing<T>::check(L, -1)) {
					lua_pop(L, 1);
					while (lua_next(L, index) != 0) lua_pop(L, 1);
					return false;
				}
				spVec->push_back(lua::LuaParsing<T>::get(L, -1));
				lua_pop(L, 1);
		   }
		   ps.push(spVec);
		}
		else return false;
	}
	return true;
}

template <typename T>
static bool PushLuaStackPointerEntryToParamStack(ParameterStack& ps, lua_State* L,
                                                 int index, const char* baseClassName,
                                                 bool bIsVector)
{
	using result_type = std::pair<T, const ClassNameNode*>;

	result_type res;
	if(!bIsVector){
		if(lua::LuaParsing<T>::checkAndGet(res, L, index, baseClassName)){
			ps.push(res.first, res.second);
		}
		else return false;
	}
	else {
		if (lua_istable(L, index)){
			SmartPtr<std::vector<result_type> > spVec
				= SmartPtr<std::vector<result_type> >(new std::vector<result_type>());
			lua_pushnil(L);
			while (lua_next(L, index) != 0) {
				if(!lua::LuaParsing<T>::checkAndGet(res, L, -1, baseClassName)) {
					lua_pop(L, 1);
					while (lua_next(L, index) != 0) lua_pop(L, 1);
					return false;
				}
				spVec->push_back(res);
				lua_pop(L, 1);
		   }
		   ps.push(spVec);
		}
		else return false;
	}
	return true;
}

}
}
#endif