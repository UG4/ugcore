/*
 * lua_stack.h
 *
 *  Created on: 01.10.2014
 *      Author: mrupp
 */

#ifndef LUA_STACK_H_
#define LUA_STACK_H_

#include "registry/registry.h"
#include "lua_parsing.h"
#include "bindings_lua.h"
#include "common/util/smart_pointer.h"

namespace ug
{
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
	typedef std::pair<T, const ClassNameNode*> result_type;

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
#endif /* LUA_STACK_H_ */
