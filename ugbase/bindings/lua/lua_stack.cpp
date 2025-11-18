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

#include "bindings_lua.h"
#include "registry/registry.h"
#include "registry/class_helper.h"
#include "common/common.h"
#include "info_commands.h"
#include "lua_util.h"
#include "lua_parsing.h"

#include "lua_stack.h"

#ifdef UG_FOR_LUA
#include "bindings/lua/lua_function_handle.h"
#include "bindings/lua/lua_table_handle.h"
#endif

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

#ifdef UG_FOR_LUA
			case Variant::VT_LUA_FUNCTION_HANDLE:{
				if(!PushLuaStackEntryToParamStack<LuaFunctionHandle>(ps, L, index, bIsVector))
					badParam = (int)i + 1;
			}break;
			case Variant::VT_LUA_TABLE_HANDLE:{
				if(!PushLuaStackEntryToParamStack<LuaTableHandle>(ps, L, index, bIsVector))
					badParam = (int)i + 1;
			}break;
#endif

			default:{//	unknown type
				badParam = (int)i + 1;
			}break;
		}

	//	check if param has been read correctly
		if(badParam){
			return badParam;
		}else{
		}
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
