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

struct UserDataWrapper
{
	void*	obj;
};

///	creates a new UserData_IObject and associates it with ptr in luas registry
/**
 * Creates a new userdata in lua, which encapsulates the given pointer.
 * It then assigns the specified metatable to the userdata.
 * When the function is done, the userdata is left on luas stack.
 */
static UserDataWrapper* CreateNewUserData(lua_State* L, void* ptr,
										  const char* metatableName)
{
//	create the userdata
	UserDataWrapper* udata = (UserDataWrapper*)lua_newuserdata(L,
											sizeof(UserDataWrapper));
		
//	associate the object with the userdata.
	udata->obj = ptr;
		
//	associate the metatable (userdata is already on the stack)
	luaL_getmetatable(L, metatableName);
	lua_setmetatable(L, -2);
	
	return udata;
}

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
			}break;
			case PT_POINTER:{
				if(lua_isuserdata(L, index)){
				//	get the object and its metatable. Make sure that obj can be cast to
				//	the type that is expected by the paramsTemplate.
					void* obj = ((UserDataWrapper*)lua_touserdata(L, index))->obj;
					if(lua_getmetatable(L, index) != 0){

						lua_pushstring(L, "names");
						lua_rawget(L, -2);
						const std::vector<const char*>* names = (const std::vector<const char*>*) lua_touserdata(L, -1);
						lua_pop(L, 2);
						bool typeMatch = false;
						if(names){
							if(!names->empty()){
								if(ClassNameVecContains(*names, paramsTemplate.class_name(i)))
									typeMatch = true;
							}
						}

						if(typeMatch)
							params.push_pointer(obj, names);
						else{
							UG_LOG("ERROR: type mismatch in argument " << i + 1);
							UG_LOG(": Expected type that supports " << paramsTemplate.class_name(i));
							bool gotone = false;
							if(names){
								if(!names->empty()){
									gotone = true;
									UG_LOG(", but given object of type " << names->at(0) << ".\n");
								}
							}

							if(!gotone){
								UG_LOG(", but given object of unknown type.\n");
							}
							printDefaultParamErrorMsg = false;
							badParam = (int)i + 1;
						}
					}
					else{
						UG_LOG("METATABLE not found. ");
						badParam = (int)i + 1;
					}
				}
				else
					badParam = (int)i + 1;
			}break;
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
			}break;
			case PT_POINTER:{
				void* obj = params.to_pointer(i);
				CreateNewUserData(L, obj, params.class_name(i));
			}break;
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
	IExportedClass* c = (IExportedClass*)lua_touserdata(L, lua_upvalueindex(1));
	
	CreateNewUserData(L, c->create(), c->name());

	return 1;
}

static int LuaProxyMethod(lua_State* L)
{
	ExportedMethod* m = (ExportedMethod*)lua_touserdata(L, lua_upvalueindex(1));

	if(!lua_isuserdata(L, 1)){
		UG_LOG("ERROR in call to LuaProxyMethod: No object specified.\n");
		return 0;
	}
	
	UserDataWrapper* self = (UserDataWrapper*)lua_touserdata(L, 1);
	
	ParameterStack paramsIn;;
	ParameterStack paramsOut;
	
	int badParam = LuaStackToParams(paramsIn, m->params_in(), L, 1);
	
//	check whether the parameter was correct
	if(badParam > 0){
		UG_LOG("ERROR occured during call to " << m->name() << endl);
		return 0;
	}

	m->execute(self->obj, paramsIn, paramsOut);
	
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
		const IExportedClass* c = &reg.get_class(i);
		
	//	set the constructor-function
		lua_pushlightuserdata(L, (void*)c);
		lua_pushcclosure(L, LuaProxyConstructor, 1);
		lua_setglobal(L, c->name());
		
	//	create the meta-table for the object
	//	overwrite index and store the class-name
		luaL_newmetatable(L, c->name());
		lua_pushvalue(L, -1);
		lua_setfield(L, -2, "__index");
		lua_pushstring(L, "names");
		lua_pushlightuserdata(L, (void*)c->class_names());
		lua_settable(L, -3);
		
		
	//	register methods
		for(size_t j = 0; j < c->num_methods(); ++j){
			const ExportedMethod& m = c->get_method(j);
			lua_pushstring(L, m.name().c_str());
			lua_pushlightuserdata(L, (void*)&m);
			lua_pushcclosure(L, LuaProxyMethod, 1);
			lua_settable(L, -3);
		}
		
	//	pop the metatable from the stack.
		lua_pop(L, 1);
	}
	
	return true;
}

}//	end of namespace
}//	end of namespace
}//	end of namespace
