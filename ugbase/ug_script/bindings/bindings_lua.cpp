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

///	This class will be used as Userdata, holding a pointer to IObject*
class UserData_IObject
{
	public:
		IObject* obj;
};

///	Returns the user-data that is associated with a given instance of IObject in luas registry.
/**	leaves the user-data on top of the stack - if there was one.*/
static UserData_IObject* GetRegisteredUserData(lua_State* L, IObject* ptr)
{
//	check whether an entry associated with ptr exists in the lua-registry
	lua_pushlightuserdata(L, ptr);
	lua_gettable(L, LUA_REGISTRYINDEX);
	UserData_IObject* udata = (UserData_IObject*)lua_touserdata(L, -1);
	if(!udata)
		lua_pop(L, 1);
	return udata;
}

static void UnregisterUserData(lua_State* L, IObject* ptr)
{
//	the pointer is the key
	lua_pushlightuserdata(L, ptr);
//	write nil to the table
	lua_pushnil(L);
//	write the table entry
	lua_settable(L, LUA_REGISTRYINDEX);
}

///	creates a new UserData_IObject and associates it with ptr in luas registry
/**
 * The correct metatable is assigned on the fly.
 * leaves the newly created userdata on the stack.
 * Note that this makes lua garbage collection impossible for IObject*.*/
static UserData_IObject* CreateAndRegisterUserData(lua_State* L, IObject* ptr)
{
//	create the userdata
	UserData_IObject* udata = (UserData_IObject*)lua_newuserdata(L, sizeof(UserData_IObject));
	
//	associate the object with the userdata.
	udata->obj = ptr;
	
//	register the udata in the table.
//	the pointer is the key
	lua_pushlightuserdata(L, ptr);
//	copy the userdata, so that one copy resides on the stack
	lua_pushvalue(L, -2);
//	create the table entry
	lua_settable(L, LUA_REGISTRYINDEX);
	
//	associate the metatable (userdata is already on the stack)
	string metaName("UGMETA_");
	metaName.append(ptr->type_name());
	luaL_getmetatable(L, metaName.c_str());
	lua_setmetatable(L, -2);
	
	return udata;
}

///	copies parameter values from the lua-stack to a parameter-list.
/**	\returns	The index of the first bad parameter starting from 1.
 *				Returns 0 if everything went right.
 */
static int LuaStackToParams(const ParameterList& params, lua_State* L,
							 int offsetToFirstParam = 0)
{
	int badParam = 0;
	bool printDefaultParamErrorMsg = true;
	
//	iterate through the parameter list and copy the value in the associated
//	stack entry.
	for(size_t i = 0; i < params.size(); ++i){
		IParameter* param = params.param(i);
		int type = param->get_type_id();
		int index = (int)i + offsetToFirstParam + 1;//stack-indices start with 1
		
		switch(type){
			case PTID_INT:{
				if(lua_isnumber(L, index))
					param->set_int(lua_tonumber(L, index));
				else
					badParam = (int)i + 1;
			}break;
			case PTID_DOUBLE:{
				if(lua_isnumber(L, index)){
					param->set_double(lua_tonumber(L, index));
				}
				else
					badParam = (int)i + 1;
			}break;
			case PTID_STRING:{
				if(lua_isstring(L, index))
					param->set_string(lua_tostring(L, index));
				else
					badParam = (int)i + 1;
			}break;
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
			}break;
			default:{//	unknown type
				badParam = (int)i + 1;
			}break;
		}
		
	//	if badParam has been set, we'll return immediatly
		if(badParam){
			if(printDefaultParamErrorMsg){
				UG_LOG("ERROR: Bad parameter " << badParam << ": ");
				UG_LOG(params.param(badParam - 1)->get_type_string() << " expected.\n");
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
static int ParamsToLuaStack(const ParameterList& params, lua_State* L)
{
//	push output parameters to the stack	
	for(size_t i = 0; i < params.size(); ++i){
		IParameter* param = params.param(i);
		int type = param->get_type_id();
		switch(type){
			case PTID_INT:{
				lua_pushnumber(L, param->to_int());
			}break;
			case PTID_DOUBLE:{
				lua_pushnumber(L, param->to_double());
			}break;
			case PTID_STRING:{
				lua_pushstring(L, param->to_string());
			}break;
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

///	called by lua whenever a global function is called
/**	The first parameter has to be the index of the IGlobalFunction object
 *	that shall be executed, followed by the methods parameters.*/
static int GlobalFunction(lua_State* L)
{
	//UG_LOG("global function called...\n");
	if(!lua_isnumber(L, 1)){
		UG_LOG("ERROR in call to GlobalFunction: No method index given.\n");
		return 0;
	}
	
	int funcIndex = lua_tointeger(L, 1);
	IGlobalFunction* func = Registry::inst().get_function(funcIndex);
	const MethodDesc& md = func->get_method();
	
//	read the parameters from the lua-stack
	const ParameterList& in = md.params_in();
	int badParam = LuaStackToParams(in, L, 1);
	
//	check whether the parameter was correct
	if(badParam > 0){
		UG_LOG("ERROR occured during call to " << md.name() << endl);
		return 0;
	}
	
//	execute the function
	func->execute();
	
//	push output parameters to Lua-Stack and return the number of parameters.
	return ParamsToLuaStack(md.params_out(), L);
}

///	called by lua whenever an objects member function is called
static int ObjectFunction(lua_State* L)
{
	if(!lua_isuserdata(L, 1)){
		UG_LOG("ERROR in call to ObjectFunction: No object specified.\n");
		return 0;
	}
	if(!lua_isstring(L, 2)){
		UG_LOG("ERROR in call to ObjectFunction: Can't extract name of metatable.\n");
		return 0;		
	}
	if(!lua_isnumber(L, 3)){
		UG_LOG("ERROR in call to ObjectFunction: No member-method index given.\n");
		return 0;
	}
	
	void* self = lua_touserdata(L, 1);
	const char* metaName = lua_tostring(L, 2);
	int methodIndex = lua_tointeger(L, 3);
	
//	make sure that the correct metatable is associated with self.
//	if it is, we can safely cast to IObject.
	if(!luaL_checkudata(L, 1, metaName)){
		UG_LOG("ERROR in call to ObjectFunction: Object passed as self has incompatible type.\n");
		return 0;
	}
	
//	everythings fine - perform the cast
	IObject* obj = reinterpret_cast<UserData_IObject*>(self)->obj;
	
	const MethodDesc& md = obj->get_method(methodIndex);
	
//	read the parameters from the lua-stack
	const ParameterList& in = md.params_in();
	int badParam = LuaStackToParams(in, L, 3);

//	check whether the parameter was correct
	if(badParam > 0){
		UG_LOG("ERROR occured during call to " << obj->type_name() << endl;);
		return 0;
	}

	obj->execute(methodIndex);
	
//	push output parameters to Lua-Stack and return the number of parameters.
	return ParamsToLuaStack(md.params_out(), L);
}

///	creates the object that is specified by the index
static int ObjectCreateFunction(lua_State* L)
{
	Registry& reg = Registry::inst();
	
	int index = lua_tointeger(L, 1);
		
//	create the object and associate it with the userdata.
	IObject* srcObj = reg.get_object(index);
	IObject* obj = srcObj->clone();
	CreateAndRegisterUserData(L, obj);

	return 1;
}

///	called by lua whenever an objects member function is called
static int ObjectDeleteFunction(lua_State* L)
{
	if(!lua_isuserdata(L, 1)){
		UG_LOG("ERROR in call to ObjectFunction: No object specified.\n");
		return 0;
	}
	
	UserData_IObject* udata = reinterpret_cast<UserData_IObject*>(lua_touserdata(L, 1));
	IObject* obj = udata->obj;
	
//	unregister the pointer and the associated userdata
	UnregisterUserData(L, obj);
	delete obj;
	udata->obj = NULL;

	return 0;
}


bool CreateBindings_LUA(lua_State* L)
{
//todo:	Instead of creating lua-methods and compiling them, one
//		should investigate lua_pushcclosure closer.

//	registers a meta-object for each object found in the ObjectRegistry.
//	Global functions are registered for all GlobalFunction-objects in the registry.

//	iterate through all registered objects
	Registry& reg = Registry::inst();
	size_t numFuncs = reg.num_functions();

//	set up the functions callbacks in lua
	lua_pushcfunction(L, GlobalFunction);
	lua_setglobal(L, "UG_GLOBAL_FUNCTION");

//	create the ug-table in which all our methods will go.
	lua_newtable(L);
//	add the delete method
	lua_pushstring(L, "delete");
	lua_pushcfunction(L, ObjectDeleteFunction);
	lua_settable(L, -3);
//	set the name of the table
	lua_setglobal(L, "ug");
	
	
////////////////////////////////
//	create a method for each global function, which calls GlobalFunction with
//	a pointer to the function object.
	for(size_t i = 0; i < numFuncs; ++i){
		const MethodDesc& md = reg.get_function(i)->get_method();
		stringstream ss;
	//	write the function head and parameters
		ss << "function ug." << md.name() << "(";

		const ParameterList& plist = md.params_in();
		for(size_t j = 0; j < plist.size(); ++j){
			ss << "p" << j;
			if(j < plist.size() - 1)
				ss << ", ";
		}
		ss << ")";
		
	//	write the function body
	//	call GlobalFunctionCalled with the function pointer and all parameters.
		ss << "return UG_GLOBAL_FUNCTION(" << i;
		for(size_t j = 0; j < plist.size(); ++j)
			ss << ", p" << j;
		ss << ") end";
	//	parse the buffer
		script::ParseBuffer(ss.str().c_str());
	}
	
////////////////////////////////
//	register classes at lua
//	first we have to setup the metatable for each object and 
//	create a creation-method for each object
	lua_pushcfunction(L, ObjectCreateFunction);
	lua_setglobal(L, "UG_CREATE_OBJECT");

	lua_pushcfunction(L, ObjectFunction);
	lua_setglobal(L, "UG_OBJECT_FUNCTION");
		
	size_t numObjs = reg.num_objects();
	for(size_t i = 0; i < numObjs; ++i){
		IObject* obj = reg.get_object(i);
		
	//	create the metatable
	//	todo: the name of the metatable should be stored, so that it can be
	//		  accessed given the index of the object.
		string metaName("UGMETA_");
		metaName.append(obj->type_name());
		luaL_newmetatable(L, metaName.c_str());
		lua_pushvalue(L, -1);
		lua_setfield(L, -2, "__index");
		lua_setglobal(L, metaName.c_str());
		
		for(size_t j = 0; j < obj->num_methods(); ++j)
		{
			const MethodDesc& md = obj->get_method(j);
			stringstream ss;
		
		//	create a function in the form of (Class and funcname are example names)
		//	function UGMETA_Class.funcname(self, p0, p1) return UG_OBJECT_FUNCTION(self, "UGMETA_Class", 2, p0, p1) end
		//	where self contains a pointer to the instance of the class
		//	on which the method is executed, "UGMETA_Class" specifies the name
		//	of the metatable that should be associated with self,
		//	and 2 is a sample index for the called method in the class.
			ss << "function " << metaName << "." << md.name() << "(self";

		//	write parameters
			const ParameterList& params = md.params_in();
			for(size_t k = 0; k < params.size(); ++k)
				ss << ", p" << k;

			ss << ") return UG_OBJECT_FUNCTION(self, " << "\"" << metaName << "\""  << ", " << j;

		//	forward parameters
			for(size_t k = 0; k < params.size(); ++k)
				ss << ", p" << k;
			
			ss << ") end";
			
		//	the function is complete and can now be parsed.
			script::ParseBuffer(ss.str().c_str());
		}
		
	//	build and parse the create-method
		lua_newtable(L);
		lua_setglobal(L, obj->type_name());
		{
			stringstream ss;
			ss << "function " << obj->type_name() << ".new() return UG_CREATE_OBJECT(" << i << ") end";
		//	parse the buffer
			script::ParseBuffer(ss.str().c_str());
		}
	}
	
	return true;
}

}//	end of namespace
}//	end of namespace
}//	end of namespace
