/*
 * lua_stack.cpp
 *
 *  Created on: 01.10.2014
 *      Author: mrupp
 */

#include "bindings_lua.h"
#include "registry/registry.h"
#include "registry/class_helper.h"
#include "common/common.h"
#include "info_commands.h"
#include "lua_util.h"
#include "lua_parsing.h"

#include "lua_stack.h"

namespace ug
{
namespace bridge {


int LuaStackToParams(ParameterStack& ps,
							const ParameterInfo& psInfo,
							lua_State* L,
							int offsetToFirstParam)
{
//	make sure that we have the right amount of parameters
//	if the sizes do not match, return -1.
	if((lua_gettop(L)) - offsetToFirstParam != (int)psInfo.size())
		return -1;

//	initialize temporary variables
	int badParam = 0;

//	iterate through parameter list and copy the value in the associated stack entry.
	for(int i = 0; i < psInfo.size(); ++i){
	//	get type and vectorMode (if bIsVector==true a lua table is expected)
		int type = psInfo.type(i);
		bool bIsVector = psInfo.is_vector(i);

	//	compute stack index, stack-indices start with 1
		int index = (int)i + offsetToFirstParam + 1;

	//	check for nil
		if(lua_type(L, index) == LUA_TNONE) return (int)i + 1;

	//	try to read in expected type
		switch(type){
			case Variant::VT_BOOL:{
				if(!PushLuaStackEntryToParamStack<bool>(ps, L, index, bIsVector))
					badParam = (int)i + 1;
			}break;
			case Variant::VT_INT:{
				if(!PushLuaStackEntryToParamStack<int>(ps, L, index, bIsVector))
					badParam = (int)i + 1;
			}break;
			case Variant::VT_SIZE_T:{
				if(!PushLuaStackEntryToParamStack<size_t>(ps, L, index, bIsVector))
					badParam = (int)i + 1;
			}break;
			case Variant::VT_FLOAT:{
				if(!PushLuaStackEntryToParamStack<float>(ps, L, index, bIsVector))
					badParam = (int)i + 1;
			}break;
			case Variant::VT_DOUBLE:{
				if(!PushLuaStackEntryToParamStack<double>(ps, L, index, bIsVector))
					badParam = (int)i + 1;
			}break;
			case Variant::VT_CSTRING:{
				if(!PushLuaStackEntryToParamStack<const char*>(ps, L, index, bIsVector))
					badParam = (int)i + 1;
			}break;
			case Variant::VT_STDSTRING:{
				if(!PushLuaStackEntryToParamStack<std::string>(ps, L, index, bIsVector))
					badParam = (int)i + 1;
			}break;
			case Variant::VT_POINTER:{
			//	NOTE that smart-ptrs are implicitly used as normal pointers here.
			//	This can cause severe problems in regard to reference-counting.
			//	This behaviour was introduced, since the registry does not
			//	allow by-value arguments. (Small temporary objects profit from
			//	this strategy).
				if(!PushLuaStackPointerEntryToParamStack<void*>
					(ps, L, index, psInfo.class_name(i), bIsVector))
					badParam = (int)i + 1;
			}break;
			case Variant::VT_CONST_POINTER:{
			//	NOTE that smart-ptrs are implicitly used as normal pointers here.
			//	This can cause severe problems in regard to reference-counting.
			//	This behaviour was introduced, since the registry does not
			//	allow by-value arguments. (Small temporary objects profit from
			//	this strategy).
				if(!PushLuaStackPointerEntryToParamStack<const void*>
					(ps, L, index, psInfo.class_name(i), bIsVector))
					badParam = (int)i + 1;
			}break;
			case Variant::VT_SMART_POINTER:{
				if(!PushLuaStackPointerEntryToParamStack<SmartPtr<void> >
					(ps, L, index, psInfo.class_name(i), bIsVector))
					badParam = (int)i + 1;
			}break;
			case Variant::VT_CONST_SMART_POINTER:{
				if(!PushLuaStackPointerEntryToParamStack<ConstSmartPtr<void> >
					(ps, L, index, psInfo.class_name(i), bIsVector))
					badParam = (int)i + 1;
			}break;

			default:{//	unknown type
				badParam = (int)i + 1;
			}break;
		}

	//	check if param has been read correctly
		if(badParam)
			return badParam;
	}

//	return result flag
	return badParam;
}




int ParamsToLuaStack(const ParameterStack& ps, lua_State* L)
{
//	push output parameters to the lua stack
	for(int i = 0; i < ps.size(); ++i){
		const int type = ps.type(i);
		const bool bIsVector = ps.is_vector(i);
		switch(type){
			case Variant::VT_BOOL:{
				ParamStackEntryToLuaStack<bool>(ps, L, i, bIsVector);
			}break;
			case Variant::VT_INT:{
				ParamStackEntryToLuaStack<int>(ps, L, i, bIsVector);
			}break;
			case Variant::VT_SIZE_T:{
				ParamStackEntryToLuaStack<size_t>(ps, L, i, bIsVector);
			}break;
			case Variant::VT_FLOAT:{
				ParamStackEntryToLuaStack<float>(ps, L, i, bIsVector);
			}break;
			case Variant::VT_DOUBLE:{
				ParamStackEntryToLuaStack<double>(ps, L, i, bIsVector);
			}break;
			case Variant::VT_CSTRING:{
				ParamStackEntryToLuaStack<const char*>(ps, L, i, bIsVector);
			}break;
			case Variant::VT_STDSTRING:{
				ParamStackEntryToLuaStack<std::string>(ps, L, i, bIsVector);
			}break;
			case Variant::VT_POINTER:{
				ParamStackPointerEntryToLuaStack<void*>(ps, L, i, bIsVector);
			}break;
			case Variant::VT_CONST_POINTER:{
				ParamStackPointerEntryToLuaStack<const void*>(ps, L, i, bIsVector);
			}break;
			case Variant::VT_SMART_POINTER:{
				ParamStackPointerEntryToLuaStack<SmartPtr<void> >(ps, L, i, bIsVector);
			}break;
			case Variant::VT_CONST_SMART_POINTER:{
				ParamStackPointerEntryToLuaStack<ConstSmartPtr<void> >(ps, L, i, bIsVector);
			}break;
			default:{
				UG_THROW("ParamsToLuaStack: invalid type in ParamStack.")
			}break;
		}
	}

	return (int)ps.size();
}

}
}
