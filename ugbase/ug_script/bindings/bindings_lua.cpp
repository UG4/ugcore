//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m09 d16

#include <sstream>
#include <cstring>
#include "bindings_lua.h"
#include "ug_script/ug_script.h"
#include "ug_interface/registry.h"
#include "common/common.h"

using namespace std;

namespace ug
{
namespace interface
{
namespace lua
{

///	copies parameter values from the lua-stack to a parameter-list.
/**	\returns	The index of the first bad parameter starting from 1.
 *				Returns 0 if everything went right.
 */
static int LuaStackToParams(ParameterStack& params,
							const ParameterStack& paramsTemplate,
							lua_State* L,
							 int offsetToFirstParam = 0)
{
	int badParam = 0;
	bool printDefaultParamErrorMsg = true;
	
//	iterate through the parameter list and copy the value in the associated
//	stack entry.
	for(int i = 0; i < paramsTemplate.size(); ++i){
		int type = paramsTemplate.get_type(i);
		int index = (int)i + offsetToFirstParam + 1;//stack-indices start with 1
		
		switch(type){
			case PT_INTEGER:{
				if(lua_isnumber(L, index))
					params.push_integer(lua_tonumber(L, index));
				else
					badParam = (int)i + 1;
			}break;
			case PT_NUMBER:{
				if(lua_isnumber(L, index)){
					params.push_number(lua_tonumber(L, index));
				}
				else
					badParam = (int)i + 1;
			}break;
			case PT_STRING:{
				if(lua_isstring(L, index))
					params.push_string(lua_tostring(L, index));
				else
					badParam = (int)i + 1;
			}break;/*
			case PTID_OBJECT:{
				if(lua_isuserdata(L, index)){
					IObject* obj = ((UserData_IObject*)lua_touserdata(L, index))->obj;
				//	check whether the given object supports the correct interface
					if(obj->type_check(params.get_interface_type(i))){
					//	types did match
						param->set_object(obj);
					}
					else{
					//	types did not match.
						UG_LOG("ERROR: type mismatch in argument " << i + 1);
						UG_LOG(": Expected type that supports " << params.get_interface_type(i));
						string strTypes;
						obj->collect_supported_types(strTypes);
						UG_LOG(", but given object supports " << strTypes << " only.\n");
						printDefaultParamErrorMsg = false;
						badParam = (int)i + 1;
					}
				}
				else
					badParam = (int)i + 1;
			}break;*/
			default:{//	unknown type
				badParam = (int)i + 1;
			}break;
		}
		
	//	if badParam has been set, we'll return immediatly
		if(badParam){
			if(printDefaultParamErrorMsg){
				UG_LOG("ERROR: Bad parameter " << badParam << ": ");
				//UG_LOG(params.param(badParam - 1)->get_type_string() << " expected.\n");
			}
			return badParam;
		}
	}
	
	return 0;
}

///	Pushes the parameter-values to the Lua-Stack.
/**
 * \returns The number of items pushed to the stack.
 */
static int ParamsToLuaStack(const ParameterStack& params, lua_State* L)
{
//	push output parameters to the stack	
	for(int i = 0; i < params.size(); ++i){
		int type = params.get_type(i);
		switch(type){
			case PT_INTEGER:{
				lua_pushnumber(L, params.to_integer(i));
			}break;
			case PT_NUMBER:{
				lua_pushnumber(L, params.to_number(i));
			}break;
			case PT_STRING:{
				lua_pushstring(L, params.to_string(i));
			}break;/*
			case PTID_OBJECT:{
			//	check whether the object has already been associated with a userdata
				IObject* obj = param->to_object();
			//	note that GetRegisteredUserData will leave the data on the stack, if
			//	it already existed
				if(!GetRegisteredUserData(L, obj)){
				//	it is not. We have to register a new one
				//	Note that the new userdata will be on the top of the stack after
				//	CreateAndRegisterUserData returns.
					CreateAndRegisterUserData(L, obj);
				}
			}break;*/
			default:{
				UG_LOG("ERROR in ParamsToLuaStack: Unknown parameter in ParameterList. ");
				UG_LOG("Return-values may be incomplete.\n");
				return (int)i;
			}break;
		}
	}
	
	return (int)params.size();
}

static int LuaProxyFunction(lua_State* L)
{
	const ExportedFunction* func = (const ExportedFunction*)lua_touserdata(L, lua_upvalueindex(1));

	ParameterStack paramsIn;;
	ParameterStack paramsOut;
	
	int badParam = LuaStackToParams(paramsIn, func->params_in(), L, 0);
	
//	check whether the parameter was correct
	if(badParam > 0){
		UG_LOG("ERROR occured during call to " << func->name() << endl);
		return 0;
	}

	func->execute(paramsIn, paramsOut);
	
	return ParamsToLuaStack(paramsOut, L);
}

static int LuaProxyConstructor(lua_State* L)
{
	const IExportedClass* inst = (const IExportedClass*)lua_touserdata(L, lua_upvalueindex(1));
	
	lua_pushlightuserdata(L, inst->create());
	
	return 1;
}

static int LuaProxyMethod(lua_State* L)
{
	ExportedMethod* m = (ExportedMethod*)lua_touserdata(L, lua_upvalueindex(1));

	if(!lua_isuserdata(L, 1)){
		UG_LOG("ERROR in call to LuaProxyMethod: No object specified.\n");
		return 0;
	}
	
	void* self = lua_touserdata(L, 1);
	
	ParameterStack paramsIn;;
	ParameterStack paramsOut;
	
	int badParam = LuaStackToParams(paramsIn, m->params_in(), L, 1);
	
//	check whether the parameter was correct
	if(badParam > 0){
		UG_LOG("ERROR occured during call to " << m->name() << endl);
		return 0;
	}

	m->execute(self, paramsIn, paramsOut);
	
	return ParamsToLuaStack(paramsOut, L);
}

bool CreateBindings_LUA(lua_State* L, InterfaceRegistry& reg)
{
//	registers a meta-object for each object found in the ObjectRegistry.
//	Global functions are registered for all GlobalFunction-objects in the registry.

//	iterate through all registered objects
	
/*
//	create the ug-table in which all our methods will go.
	lua_newtable(L);
//	set the name of the table
	lua_setglobal(L, "ug");
	*/
////////////////////////////////
//	create a method for each global function, which calls LuaProxyFunction with
//	a an index to the function.
	size_t numFuncs = reg.num_functions();
	
	for(size_t i = 0; i < numFuncs; ++i){
		ExportedFunction* func = &reg.get_function(i);

		lua_pushlightuserdata(L, func);
		lua_pushcclosure(L, LuaProxyFunction, 1);
		lua_setglobal(L, func->name().c_str());
	}


	size_t numClasses = reg.num_classes();
	for(size_t i = 0; i < numClasses; ++i){
		const IExportedClass* inst = &reg.get_class(i);
		lua_pushlightuserdata(L, (void*)inst);
		lua_pushcclosure(L, LuaProxyConstructor, 1);
		lua_setglobal(L, inst->name());
		
		lua_newtable(L);
	//	register methods
		for(size_t j = 0; j < inst->num_methods(); ++j){
			const ExportedMethod& m = inst->get_method(j);
			lua_pushstring(L, m.name().c_str());
			lua_pushlightuserdata(L, (void*)&m);
			lua_pushcclosure(L, LuaProxyMethod, 1);
			lua_settable(L, -3);
		}
	//	set the name of the table
		lua_setglobal(L, "tmp");
	}
	
	return true;
}

}//	end of namespace
}//	end of namespace
}//	end of namespace
