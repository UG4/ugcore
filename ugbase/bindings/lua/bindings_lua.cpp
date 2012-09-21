//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m09 d16

#include <sstream>
#include <cstring>
#include <cstdlib>
#include <queue>
#include "bindings_lua.h"
#include "registry/registry.h"
#include "registry/class_helper.h"
#include "common/common.h"
#include "info_commands.h"
#include "lua_util.h"

//#define __UG__BINDINGS_LUA__CATCH_UNKNOWN_EXCEPTIONS__

using namespace std;
using namespace ug::script;

///	throw mechanism for lua related errors.
#define UG_LUA_THROW_EMPTY(luaState)	luaL_error(luaState, "")
#define UG_LUA_THROW(luaState, msg)		luaL_error(luaState, "\n%s", msg)

namespace ug{
namespace bridge{
namespace lua{

//	set this variable to true if smart-ptr arguments shall be automatically
//	converted to raw-ptrs where required.
const bool IMLPICIT_SMART_PTR_TO_PTR_CONVERSION = true;


//	a symbol preceding error messages
const char* errSymb = " % ";


///	creates a new UserDataWrapper and associates it with ptr in luas registry
/**
 * Creates a new userdata in lua, which encapsulates the given pointer / smart-pointer.
 * It then assigns the specified metatable to the userdata.
 * When the function is done, the userdata is left on luas stack.
 * \{
 */
static SmartUserDataWrapper* CreateNewUserData(lua_State* L, const SmartPtr<void>& ptr,
											  const char* metatableName)
{
//	create the userdata
	SmartUserDataWrapper* udata = (SmartUserDataWrapper*)lua_newuserdata(L,
											sizeof(SmartUserDataWrapper));
	new(udata) SmartUserDataWrapper;

//	associate the object with the userdata.
	udata->smartPtr = ptr;
	udata->type = SMART_POINTER;

//	associate the metatable (userdata is already on the stack)
	luaL_getmetatable(L, metatableName);
	lua_setmetatable(L, -2);

	return udata;
}

static ConstSmartUserDataWrapper* CreateNewUserData(lua_State* L, const ConstSmartPtr<void>& ptr,
											  const char* metatableName)
{
//	create the userdata
	ConstSmartUserDataWrapper* udata = (ConstSmartUserDataWrapper*)lua_newuserdata(L,
											sizeof(ConstSmartUserDataWrapper));

	new(udata) ConstSmartUserDataWrapper;

//	associate the object with the userdata.
	udata->smartPtr = ptr;
	udata->type = SMART_POINTER | IS_CONST;

//	associate the metatable (userdata is already on the stack)
	luaL_getmetatable(L, metatableName);
	lua_setmetatable(L, -2);

	return udata;
}

static RawUserDataWrapper* CreateNewUserData(lua_State* L, void* ptr,
										  const char* metatableName,
										  void (*deleteFunc)(const void*),
										  bool is_const)
{
//	create the userdata
	RawUserDataWrapper* udata = (RawUserDataWrapper*)lua_newuserdata(L,
											sizeof(RawUserDataWrapper));

	new(udata) RawUserDataWrapper;

//	associate the object with the userdata.
	udata->obj = ptr;
	udata->deleteFunc = deleteFunc;
	udata->type = RAW_POINTER;
	if(is_const)
		udata->type |= IS_CONST;
	
//	associate the metatable (userdata is already on the stack)
	luaL_getmetatable(L, metatableName);
	lua_setmetatable(L, -2);
	
	return udata;
}
/** \} */

/**
 *
 * \brief returns a String describing the parameters on the lua stack
 * ex. "GlobalMultiGridRefiner*, LuaUserNumber2d*, number, string"
 */
// returns a String describing the parameters on the lua stack
static string GetLuaParametersString(lua_State* L, int offsetToFirstParam = 0)
{
	string str;
	bool bFirst=true;
	int index = offsetToFirstParam + 1; // stack-indices start with 1
	for(; lua_type(L, index) != LUA_TNONE; index++)
	{
		if(!bFirst) str.append(", "); else bFirst = false;
		str.append(GetLuaTypeString(L, index));
	}
	return str;
}

/**
 *
 * \returns	String describing the reason why LuaStackToParams failed.
 * \param paramsTempalte
 * \param L
 * \param offsetToFirstParam
 * \param badParamOneBased : return value as in LuaStackParams
 * \sa LuaStackToParams
 */
static string GetTypeMismatchString(const ParameterStack& paramsTemplate,
									lua_State* L, int offsetToFirstParam,
									int badParamOneBased)
{
	std::stringstream ss;

	if(badParamOneBased == -1)
		ss << "number of parameters did not match (got "
			<< lua_gettop(L) - offsetToFirstParam
			<< ", but needs " << paramsTemplate.size() << ").";
	else
	{
		int i = badParamOneBased-1; // i is zero-based.
		int index = (int)i + offsetToFirstParam + 1;
		ss << "type mismatch in argument " << badParamOneBased
			<< ": expected " << ParameterToString(paramsTemplate, i)
			<< ", but given " << GetLuaTypeString(L, index);
	}
	return ss.str();
}

///	copies parameter values from the lua-stack to a parameter-list.
/**	\returns	The index of the first bad parameter starting from 1.
 *				Returns 0 if everything went right.
 *				Returns -1 if the number of parameters did not match.
 */
static int LuaStackToParams(ParameterStack& params,
							const ParameterStack& paramsTemplate,
							lua_State* L,
							int offsetToFirstParam = 0)
{
//	I disabled all output, since we might have overloads (sreiter).

//	make sure that we have the right amount of parameters
//	if the sizes do not match, return -1.
	if((lua_gettop(L)) - offsetToFirstParam != (int)paramsTemplate.size())
		return -1;

//	initialize temporary variables
	int badParam = 0;
	// MR: commented out, unused?
	//bool printDefaultParamErrorMsg = true;	
	
//	iterate through the parameter list and copy the value in the associated
//	stack entry.
	for(int i = 0; i < paramsTemplate.size(); ++i){
		int type = paramsTemplate.get_type(i);
		int index = (int)i + offsetToFirstParam + 1;//stack-indices start with 1
		
		if(lua_type(L, index) == LUA_TNONE)
			return (int)i + 1;

		switch(type){
			case PT_BOOL:{
				if(lua_isboolean(L, index))
					params.push_bool(lua_toboolean(L, index));
				else badParam = (int)i + 1;

			}break;
			case PT_INTEGER:{
				if(lua_isnumber(L, index))
					params.push_integer(lua_tointeger(L, index));
				else{
					badParam = (int)i + 1;
				}
			}break;
			case PT_NUMBER:{
				if(lua_isnumber(L, index)){
					params.push_number(lua_tonumber(L, index));
				}
				else{
					badParam = (int)i + 1;
				}
			}break;
			case PT_CSTRING:{
				if(lua_isstring(L, index))
					params.push_cstring(lua_tostring(L, index));
				else{
					badParam = (int)i + 1;
				}
			}break;
			case PT_STD_STRING:{
				if(lua_isstring(L, index))
					params.push_std_string(lua_tostring(L, index));
				else{
					badParam = (int)i + 1;
				}
			}break;
			case PT_POINTER:{
			//	NOTE that smart-ptrs are implicitly used as normal pointers here.
			//	This can cause severe problems in regard to reference-counting.
			//	This behaviour was introduced, since the registry does not
			//	allow by-value arguments. (Small temporary objects profit from
			//	this strategy).
				if(lua_isuserdata(L, index)){
					void* obj = NULL;
					UserDataWrapper* udata =
						reinterpret_cast<UserDataWrapper*>(lua_touserdata(L, index));

					if(udata->is_const()){
						badParam = (int)i + 1;
						break;
					}

				//	extract the pointer to the object.
				//	udata is either a RawUserData or a SmartUserDataWrapper
					if(udata->is_raw_ptr())
						obj = static_cast<RawUserDataWrapper*>(udata)->obj;
					else if(udata->is_smart_ptr() && IMLPICIT_SMART_PTR_TO_PTR_CONVERSION)
						obj = static_cast<SmartUserDataWrapper*>(udata)->smartPtr.get();
					else{
						badParam = (int)i + 1;
						break;
					}
					
				//	get the object and its metatable. Make sure that obj can be cast to
				//	the type that is expected by the paramsTemplate.
					if(lua_getmetatable(L, index) != 0){

						lua_pushstring(L, "class_name_node");
						lua_rawget(L, -2);
						const ClassNameNode* classNameNode = (const ClassNameNode*) lua_touserdata(L, -1);
						lua_pop(L, 2);
						bool typeMatch = false;
						if(classNameNode){
							if(!classNameNode->empty()){
								if(ClassNameTreeContains(*classNameNode, paramsTemplate.class_name(i)))
									typeMatch = true;
							}
						}

						if(typeMatch)
							params.push_pointer(obj, classNameNode);
						else{
							//printDefaultParamErrorMsg = false;
							badParam = (int)i + 1;
						}
					}
					else{
						badParam = (int)i + 1;
					}
				}
				else
				{
					badParam = (int)i + 1;
				}
			}break;
			case PT_CONST_POINTER:{
			//	NOTE that smart-ptrs are implicitly used as normal pointers here.
			//	This can cause severe problems in regard to reference-counting.
			//	This behaviour was introduced, since the registry does not
			//	allow by-value arguments. (Small temporary objects profit from
			//	this strategy).
				if(lua_isuserdata(L, index)){
				//	extract the object-pointer from the user-data
					const void* obj = NULL;
					UserDataWrapper* udata =
						reinterpret_cast<UserDataWrapper*>(lua_touserdata(L, index));

					if(udata->is_raw_ptr())
						obj = static_cast<RawUserDataWrapper*>(udata)->obj;
					else if(udata->is_smart_ptr() && IMLPICIT_SMART_PTR_TO_PTR_CONVERSION){
					//	we have to distinguish between const and non-const smart pointers.
						if(udata->is_const())
							obj = static_cast<ConstSmartUserDataWrapper*>(udata)->smartPtr.get();
						else
							obj = static_cast<SmartUserDataWrapper*>(udata)->smartPtr.get();
					}
					else{
						badParam = (int)i + 1;
						break;
					}

				//	get the objects metatable. Make sure that obj can be cast to
				//	the type that is expected by the paramsTemplate.
					if(lua_getmetatable(L, index) != 0){

						lua_pushstring(L, "class_name_node");
						lua_rawget(L, -2);
						const ClassNameNode* classNameNode = (const ClassNameNode*) lua_touserdata(L, -1);
						lua_pop(L, 2);
						bool typeMatch = false;
						if(classNameNode){
							if(!classNameNode->empty()){
								if(ClassNameTreeContains(*classNameNode, paramsTemplate.class_name(i)))
									typeMatch = true;
							}
						}

						if(typeMatch)
							params.push_const_pointer(obj, classNameNode);
						else{
							//printDefaultParamErrorMsg = false;
							badParam = (int)i + 1;
						}
					}
					else{
						badParam = (int)i + 1;
					}
				}
				else
				{
					badParam = (int)i + 1;
				}
			}break;
			case PT_SMART_POINTER:{
				if(lua_isuserdata(L, index)){
				//	Check whether this is really a smart pointer
					if(!((UserDataWrapper*)lua_touserdata(L, index))->is_smart_ptr()){
						badParam = (int)i + 1;
						break;
					}

				//	if the object is a const object, we can't use it here.
					if(((UserDataWrapper*)lua_touserdata(L, index))->is_const()){
						badParam = (int)i + 1;
						break;
					}

				//	get the object and its metatable. Make sure that obj can be cast to
				//	the type that is expected by the paramsTemplate.
					SmartPtr<void>& obj = ((SmartUserDataWrapper*)lua_touserdata(L, index))->smartPtr;
					if(lua_getmetatable(L, index) != 0){

						lua_pushstring(L, "class_name_node");
						lua_rawget(L, -2);
						const ClassNameNode* classNameNode = (const ClassNameNode*) lua_touserdata(L, -1);
						lua_pop(L, 2);
						bool typeMatch = false;
						if(classNameNode){
							if(!classNameNode->empty()){
								if(ClassNameTreeContains(*classNameNode, paramsTemplate.class_name(i)))
									typeMatch = true;
							}
						}

						if(typeMatch)
							params.push_smart_pointer(obj, classNameNode);
						else{
							/*bool gotone = false;
							if(classNameNode){
								if(!classNameNode->empty()){
									gotone = true;
								}
							}
							printDefaultParamErrorMsg = false;*/
							badParam = (int)i + 1;
						}
					}
					else{
						badParam = (int)i + 1;
					}
				}
				else{
					badParam = (int)i + 1;
				}
			}break;
			case PT_CONST_SMART_POINTER:{
				if(lua_isuserdata(L, index)){
				//	Check whether this is really a smart pointer
					if(!((UserDataWrapper*)lua_touserdata(L, index))->is_smart_ptr()){
						badParam = (int)i + 1;
						break;
					}

				//	get the object and its metatable. Make sure that obj can be cast to
				//	the type that is expected by the paramsTemplate.
					ConstSmartPtr<void> obj;
					if(((UserDataWrapper*)lua_touserdata(L, index))->is_const())
						obj = ((ConstSmartUserDataWrapper*)lua_touserdata(L, index))->smartPtr;
					else
						obj = ((SmartUserDataWrapper*)lua_touserdata(L, index))->smartPtr;

					if(lua_getmetatable(L, index) != 0){

						lua_pushstring(L, "class_name_node");
						lua_rawget(L, -2);
						const ClassNameNode* classNameNode = (const ClassNameNode*) lua_touserdata(L, -1);
						lua_pop(L, 2);
						bool typeMatch = false;
						if(classNameNode){
							if(!classNameNode->empty()){
								if(ClassNameTreeContains(*classNameNode, paramsTemplate.class_name(i)))
									typeMatch = true;
							}
						}

						if(typeMatch)
							params.push_const_smart_pointer(obj, classNameNode);
						else{
							/*bool gotone = false;
							if(classNameNode){
								if(!classNameNode->empty()){
									gotone = true;
								}
							}
							printDefaultParamErrorMsg = false;*/
							badParam = (int)i + 1;
						}
					}
					else{
						badParam = (int)i + 1;
					}
				}
				else{
					badParam = (int)i + 1;
				}
			}break;

			default:{//	unknown type
				badParam = (int)i + 1;
			}break;
		}
		if(badParam)
			return badParam;
	}
		
	return badParam;
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
			case PT_BOOL:{
				lua_pushboolean(L, (params.to_bool(i)) ? 1 : 0);
			}break;
			case PT_INTEGER:{
				lua_pushnumber(L, params.to_integer(i));
			}break;
			case PT_NUMBER:{
				lua_pushnumber(L, params.to_number(i));
			}break;
			case PT_CSTRING:{
				lua_pushstring(L, params.to_cstring(i));
			}break;
			case PT_STD_STRING:{
				lua_pushstring(L, params.to_std_string(i).c_str());
			}break;
			case PT_POINTER:{
				void* obj = params.to_pointer(i);
				CreateNewUserData(L, obj, params.class_name(i), NULL, false);
			}break;
			case PT_CONST_POINTER:{
			//	we're removing const with a cast. However, it was made sure that
			//	obj is treated as a const value.
				void* obj = (void*)params.to_const_pointer(i);
				CreateNewUserData(L, obj, params.class_name(i), NULL, true);
			}break;
			case PT_SMART_POINTER:{
				CreateNewUserData(L, params.to_smart_pointer(i), params.class_name(i));
			}break;
			case PT_CONST_SMART_POINTER:{
				CreateNewUserData(L, params.to_const_smart_pointer(i), params.class_name(i));
			}break;
			default:{
				return (int)i;
			}break;
		}
	}
	
	return (int)params.size();
}


static void PrintUGErrorTraceback(UGError &err)
{
//	header
	UG_LOG(errSymb<<"  Error traceback (innermost first): \n");

//	padding to insert
	std::string pad(errSymb); pad.append("     ");

//	print each message
	for(size_t i=0;i<err.num_msg();++i)
	{
	//	get copy of original sting
		std::string msg = err.get_msg(i);

	//	add paddings
		std::string::size_type pos = 0;
		while (1) {
		    pos = msg.find('\n', pos);
		    if (pos == std::string::npos) break;
		    pos++;
		    msg.insert(pos, pad);
		}

	//	write message
		UG_LOG(errSymb<<std::setw(3)<<i<<": "<<msg<<endl);

	//	write file and line
		UG_LOG(pad << "[at "<<err.get_file(i)<<", line "<<err.get_line(i)<<"]\n");
	}
}

//	global functions are handled here
//	Note that not the best matching, but the first matchin overload is chosen!
static int LuaProxyFunction(lua_State* L)
{
	const ExportedFunctionGroup* funcGrp = (const ExportedFunctionGroup*)
											lua_touserdata(L, lua_upvalueindex(1));
//	we have to try each overload!
	int badParam = -2;
	for(size_t i = 0; i < funcGrp->num_overloads(); ++i){
		const ExportedFunction* func = funcGrp->get_overload(i);

		ParameterStack paramsIn;
		ParameterStack paramsOut;

		badParam = LuaStackToParams(paramsIn, func->params_in(), L, 0);

	//	check whether the parameter was correct
		if(badParam != 0){
		//	parameters didn't match. Try the next overload.
			continue;
		}

		try{
			func->execute(paramsIn, paramsOut);
		}
		catch(LuaError& err){
			UG_LUA_THROW(L, err.get_msg().c_str());
		}
		catch(UGError& err){
			UG_LOG(errSymb<<"Error at " << GetLuaFileAndLine(L) << ":\n");
			UG_LOG(errSymb<<"UGError thrown in call to function '");
			PrintFunctionInfo(*func); UG_LOG("'.\n");
			PrintUGErrorTraceback(err);

			//UG_LOG(errSymb<<"Call stack:\n");	LuaStackTrace(L);
			UG_LUA_THROW_EMPTY(L);
		}
		catch(bad_alloc& ba)
		{
			UG_LOG(errSymb<<"Error at " << GetLuaFileAndLine(L) << ":\n");
			UG_LOG(errSymb<<"std::bad_alloc thrown in call to function '");
			PrintFunctionInfo(*func); UG_LOG("'.\n");
			UG_LOG(errSymb<<"bad_alloc description: " << ba.what() << endl);

			//UG_LOG(errSymb<<"Call stack:\n");	LuaStackTrace(L);
			UG_LUA_THROW_EMPTY(L);
		}
#ifdef __UG__BINDINGS_LUA__CATCH_UNKNOWN_EXCEPTIONS__
		catch(...)
		{
			UG_LOG(errSymb<<"Error at " << GetLuaFileAndLine(L) << ":\n");
			UG_LOG(errSymb<<"Unknown Exception thrown in call to function '");
			PrintFunctionInfo(*func); UG_LOG("'.\n");
			UG_LUA_THROW_EMPTY(L);
		}
#endif

	//	if we reach this point, then the method was successfully executed.
		return ParamsToLuaStack(paramsOut, L);
	}

	if(badParam != 0)
	{
		UG_LOG(errSymb<<"Error at "<<GetLuaFileAndLine(L) << ":\n");
		UG_LOG(errSymb<<"ERROR occured when trying to call '"
			   << funcGrp->name() << "(" << GetLuaParametersString(L, 0) << "):'\n");
		UG_LOG(errSymb<<"No matching overload found! Candidates are:\n");
		for(size_t i = 0; i < funcGrp->num_overloads(); ++i)
		{
			const ExportedFunction* func = funcGrp->get_overload(i);
			ParameterStack paramsIn;
			badParam = LuaStackToParams(paramsIn, func->params_in(), L, 0);
			UG_LOG(errSymb<<" - ");
			PrintFunctionInfo(*func);
			UG_LOG(": " << GetTypeMismatchString(func->params_in(), L, 0, badParam) << "\n");
		}
		//UG_LOG(errSymb<<"Call stack:\n");	LuaStackTrace(L);

		UG_LUA_THROW_EMPTY(L);
	}

//	this point shouldn't be reached
	UG_LUA_THROW(L, "Unknown internal error!");
	return 0;
}


static int LuaConstructor(lua_State* L, IExportedClass* c, const char *groupname=NULL)
{
//	try each constructor overlaod
	int badParam = -2;
	for(size_t i = 0; i < c->num_constructors(); ++i)
	{
	//	get overload
		const ExportedConstructor& constr = c->get_constructor(i);

		ParameterStack paramsIn;
		badParam = LuaStackToParams(paramsIn, constr.params_in(), L, 0);

	//	check whether the parameter was correct
		if(badParam != 0)
		{
		//	parameters didn't match. Try the next overload.
			continue;
		}

		try{
			if(c->construct_as_smart_pointer()){
				CreateNewUserData(L,
					SmartPtr<void>(constr.create(paramsIn), c->get_delete_function()),
					c->name().c_str());
			}
			else{
				CreateNewUserData(L, constr.create(paramsIn), c->name().c_str(),
								  c->get_delete_function(), false);
			}
		}
		catch(LuaError& err){
			UG_LUA_THROW(L, err.get_msg().c_str());
		}
		catch(UGError& err)
		{
			UG_LOG(errSymb<<"Error in " << GetLuaFileAndLine(L) <<":\n");
			UG_LOG(errSymb<<"UGError thrown while creating class '"<<c->name());
			UG_LOG("'.\n");
			PrintUGErrorTraceback(err);

			//UG_LOG(errSymb<<"Call stack:\n"); LuaStackTrace(L);
			UG_LUA_THROW_EMPTY(L);
		}

	//	object created
		return 1;
	}

//	no matching overload found
	if(badParam < 0)
	{
		UG_LOG(errSymb<<"Error at "<<GetLuaFileAndLine(L) << ":\n");
		UG_LOG(errSymb<<"ERROR occured when trying to create object of ");
		if(groupname)
		{ 	UG_LOG("group " << groupname << " (default class " << c->name() << ") with constructor '" << c->name());	}
		else
		{ 	UG_LOG("class " << c->name() << " with constructor '" << c->name());	}
		UG_LOG("(" << GetLuaParametersString(L, 0) << ")':\n");
		UG_LOG(errSymb<<"No matching overload found! Candidates are:\n");
		for(size_t i = 0; i < c->num_constructors(); ++i)
		{
			const ExportedConstructor& constr = c->get_constructor(i);
			ParameterStack paramsIn;
			badParam = LuaStackToParams(paramsIn, constr.params_in(), L, 0);
			UG_LOG(errSymb<<" - ");
			PrintConstructorInfo(constr, c->name().c_str());
			UG_LOG(": " << GetTypeMismatchString(constr.params_in(), L, 0, badParam) << "\n");
		}
		//UG_LOG(errSymb<<"Call stack:\n");	LuaStackTrace(L);
		UG_LUA_THROW_EMPTY(L);
	}

	//UG_LOG(errSymb<<"Call stack:\n"); LuaStackTrace(L);
	UG_LUA_THROW_EMPTY(L);

	return 0;
}

static int LuaProxyConstructor(lua_State* L)
{
//	get class
	IExportedClass* c = (IExportedClass*)lua_touserdata(L, lua_upvalueindex(1));
	return LuaConstructor(L, c);
}


//	creates the class which is set as default class for the specified group.
//	we assume that the first upvalue is a ClassGroupDesc*
static int LuaProxyGroupConstructor(lua_State* L)
{
//	get the group and make sure that it contains data
	const ClassGroupDesc* group = (ClassGroupDesc*)lua_touserdata(L, lua_upvalueindex(1));

	if(group->empty()){
		UG_LOG(errSymb<<"Error at "<<GetLuaFileAndLine(L) << ":\n")
		UG_LOG(errSymb<<"Can't create default instance of group '" << group->name());
		UG_LOG("': Group is empty!\n");
		lua_pushnil(L);
		return 1;
	}

//	get the associated default class
	IExportedClass* c = group->get_default_class();

	if(!c){
		UG_LOG(errSymb<<"Error at "<<GetLuaFileAndLine(L) << ":\n")
		UG_LOG(errSymb<<"Can't create default instance of group '" << group->name());
		UG_LOG("': No default class set!\n");
		lua_pushnil(L);
		return 1;
	}

	return LuaConstructor(L, c, group->name().c_str());
}

//	This method is not called by lua, but a helper to LuaProxyMethod.
//	It recursivly calls itself until a matching overload was found.
static int ExecuteMethod(lua_State* L, const ExportedMethodGroup* methodGrp,
						UserDataWrapper* self, const ClassNameNode* classNameNode,
						bool errorOutput)
{
	ParameterStack paramsIn;;
	ParameterStack paramsOut;

	//int badParam = LuaStackToParams(paramsIn, m->params_in(), L, 1);

//	we have to try each overload!
	int badParam = -2;
	for(size_t i = 0; i < methodGrp->num_overloads(); ++i){
		const ExportedMethod* m = methodGrp->get_overload(i);

		ParameterStack paramsIn;
		ParameterStack paramsOut;

		badParam = LuaStackToParams(paramsIn, m->params_in(), L, 1);

	//	check whether the parameter was correct
		if(badParam != 0){
		//	parameters didn't match. Try the next overload.
			continue;
		}

		try
		{
		//	raw pointer
			if(self->is_raw_ptr())
			{
			//	cast to the needed base class
				void* objPtr = ClassCastProvider::cast_to_base_class(
											((RawUserDataWrapper*)self)->obj,
											classNameNode, m->class_name().c_str());

				m->execute(objPtr, paramsIn, paramsOut);
			}
		//	smart pointer
			else if(self->is_smart_ptr())
			{
				if(self->is_const())
				{
				//	cast to the needed base class
					void* objPtr = ClassCastProvider::cast_to_base_class(
												(void*)((ConstSmartUserDataWrapper*)self)->smartPtr.get(),
												classNameNode, m->class_name().c_str());

					m->execute(objPtr, paramsIn, paramsOut);
				}
				else
				{
				//	cast to the needed base class
					void* objPtr = ClassCastProvider::cast_to_base_class(
												((SmartUserDataWrapper*)self)->smartPtr.get(),
												classNameNode, m->class_name().c_str());

					m->execute(objPtr, paramsIn, paramsOut);
				}
			}
		}
		catch(LuaError& err){
			UG_LUA_THROW(L, err.get_msg().c_str());
		}
		catch(UGError& err)
		{
			UG_LOG(errSymb << GetLuaFileAndLine(L) << ":\n");
			UG_LOG(errSymb << "UGError thrown in call to method '");
			PrintLuaClassMethodInfo(L, 1, *m); UG_LOG("'.\n");
			PrintUGErrorTraceback(err);
			UG_LUA_THROW_EMPTY(L);
		}
		catch(std::bad_alloc& ba)
		{
			UG_LOG(errSymb << GetLuaFileAndLine(L) << ":\n");
			UG_LOG(errSymb << "std::bad_alloc thrown in call to '");
			UG_LOG(errSymb<<"bad_alloc description: " << ba.what() << endl);
			PrintLuaClassMethodInfo(L, 1, *m); UG_LOG(".'\n");
			UG_LUA_THROW_EMPTY(L);
		}
#ifdef 	__UG__BINDINGS_LUA__CATCH_UNKNOWN_EXCEPTIONS__
		catch(...)
		{
			UG_LOG(errSymb << GetLuaFileAndLine(L) << ":\n");
			UG_LOG(errSymb << "Unknown Exception thrown in call to '");
			PrintLuaClassMethodInfo(L, 1, *m); UG_LOG("'.\n");
			UG_LUA_THROW_EMPTY(L);
		}
#endif

	//	if we reach this point, then the method was successfully executed.
		return ParamsToLuaStack(paramsOut, L);
	}

//	check whether the parameters were correct
	if(badParam != 0)
	{
	//	they were not. If the class has a base class, then we can try to
	//	to find a method-group in one of the base classes and recursively
	//	call this method.
		if(classNameNode != NULL){
		//	check whether a base-class contains overloads of this method-group
		//	push all base classes to this queue of class name nodes
			std::queue<const ClassNameNode*> qClassNameNodes;
			for(size_t i = 0; i < classNameNode->num_base_classes(); ++i)
				qClassNameNodes.push(&classNameNode->base_class(i));

		//	now visit the whole base-class hierarchy to find overloads of
		//	this method. Stop if one was successfully executed.
			while(!qClassNameNodes.empty()){
				const ClassNameNode* curClassName = qClassNameNodes.front();
				qClassNameNodes.pop();

			//	get the metatable of this class
				luaL_getmetatable(L, curClassName->name().c_str());

			//	check whether the metatable contains a method-group with
			//	the given name
				const ExportedMethodGroup* newMethodGrp = NULL;
				if(!self->is_const()){
				//	access the table which stores method-groups
					lua_pushstring(L, "__method_grps");
					lua_rawget(L, -2);

					if(lua_istable(L, -1)){
						lua_pushstring(L, methodGrp->name().c_str());
						lua_rawget(L, -2);

					//	if we retrieved something != nil, we've found one.
						if(!lua_isnil(L, -1)){
							newMethodGrp = (const ExportedMethodGroup*)
											lua_touserdata(L, -1);
						}

					//	pop the result
						lua_pop(L, 1);
					}
				//	pop the table
					lua_pop(L, 1);
				}

			//	if the object is const or if no non-const member was found,
			//	we'll check the const methods
				if(!newMethodGrp){
				//	the method is const
					lua_pushstring(L, "__const_method_grps");
					lua_rawget(L, -2);

					if(lua_istable(L, -1)){
					//	check whether the entry is contained in the table
						lua_pushstring(L, methodGrp->name().c_str());
						lua_rawget(L, -2);

					//	if we retrieved something != nil, we're done.
						if(!lua_isnil(L, -1)){
							newMethodGrp = (const ExportedMethodGroup*)
											lua_touserdata(L, -1);
						}

					//	remove result from stack
						lua_pop(L, 1);
					}
				//	remove __const table from stack
					lua_pop(L, 1);
				}

			//	remove metatable from stack
				lua_pop(L, 1);

			//	if we found a base-implementation, call it now.
			//	if not, add all base-classes to the queue again.
			//	NOTE: If a base class contains the implementation, we
			//	don't have to add it to the queue, since the method
			//	is recursive.
				if(newMethodGrp){
					int retVal = ExecuteMethod(L, newMethodGrp, self,
												curClassName, errorOutput);
					if(retVal >= 0)
						return retVal;
				}
				else{
					for(size_t i = 0; i < curClassName->num_base_classes(); ++i)
						qClassNameNodes.push(&curClassName->base_class(i));
				}
			}

		}

	//	neither the given class nor one of its base classes contains a matching
	//	overload of the given method. We thus have to output errors.
	//	Here we only print the overload-infos. The rest is done in LuaProxyMethod.
		if(errorOutput){
			for(size_t i = 0; i < methodGrp->num_overloads(); ++i)
			{
				const ExportedMethod* func = methodGrp->get_overload(i);
				ParameterStack paramsIn;
				badParam = LuaStackToParams(paramsIn, func->params_in(), L, 1);
				UG_LOG("- ");
				PrintFunctionInfo(*func);
				UG_LOG(": " << GetTypeMismatchString(func->params_in(), L, 1, badParam) << "\n");
			}
		}
	}
	
	return -1;
}

//	member methods of classes are handled here
static int LuaProxyMethod(lua_State* L)
{
	const ExportedMethodGroup* methodGrp = (const ExportedMethodGroup*)
											lua_touserdata(L, lua_upvalueindex(1));

	if(!lua_isuserdata(L, 1))
	{
		UG_LOG(errSymb<<"Error at "<<GetLuaFileAndLine(L) << ":\n")
		UG_LOG(errSymb<<"Error in call to LuaProxyMethod: No object specified in call to '");
		PrintLuaClassMethodInfo(L, 1, *methodGrp->get_overload(0));
		UG_LOG("'.\n");
		return 0;
	}

	UserDataWrapper* self = (UserDataWrapper*)lua_touserdata(L, 1);

//	get metatable of object and extract the class name node
	lua_getmetatable(L, 1);
	lua_pushstring(L, "class_name_node");
	lua_rawget(L, -2);
	const ClassNameNode* classNameNode
		= (const ClassNameNode*) lua_touserdata(L, -1);
	lua_pop(L, 2);

	int retVal = ExecuteMethod(L, methodGrp, self, classNameNode, false);
	if(retVal >= 0)
		return retVal;

//	The call failed. We have to output errors
	const char *classname = "(unknown class)";
	if(classNameNode != NULL)
		classname = classNameNode->name().c_str();

	UG_LOG(errSymb<<"Error at "<<GetLuaFileAndLine(L) << ":\n");
	UG_LOG(errSymb<<"There is no member function '"
		   << classname << ":" << methodGrp->name() << "(" <<
		   GetLuaParametersString(L, 1) << ")':\n");

	UG_LOG(errSymb<<"No matching overload found! Candidates in class "
			<< classname << " are:\n");

//	since no matching method was found, we'll run ExecuteMethod again, this
//	time outputting the errors
	ExecuteMethod(L, methodGrp, self, classNameNode, true);

	//UG_LOG(errSymb<<"Call stack:\n"); LuaStackTrace(L);
	UG_LUA_THROW_EMPTY(L);

	return 0;
}

//TODO:	The metatable indexer could probably be integrated into ExecuteMethod,
//		(called by LuaProxyMethod) and thus be avoided completely.
//		Note that a recursion over the base classes is performed in
//		ExecuteMethod, too.
static int MetatableIndexer(lua_State*L)
{
//	the stack contains the object and the requested key (the method name).
//	we have to make sure to only call const methods on const objects
	bool is_const = ((UserDataWrapper*)lua_touserdata(L, 1))->is_const();
	
//	first we push the objects metatable onto the stack.
	lua_getmetatable(L, 1);

//	queue of class name nodes
	std::queue<const ClassNameNode*> qClassNameNodes;

//	now we have to traverse the class hierarchy
	while(1)
	{
	//	if the object is not const, check whether the key is in the current metatable
		if(!is_const){
			lua_pushvalue(L, 2);
			lua_rawget(L, -2);
		
		//	if we retrieved something != nil, we're done.
			if(!lua_isnil(L, -1))
				return 1;
	
			lua_pop(L, 1);
		}
		
	//	the next thing to check are the const methods
		lua_pushstring(L, "__const");
		lua_rawget(L, -2);
		if(!lua_isnil(L, -1)){
		//	check whether the entry is contained in the table
			lua_pushvalue(L, 2);
			lua_rawget(L, -2);
		
		//	if we retrieved something != nil, we're done.
			if(!lua_isnil(L, -1)){
				return 1;
			}
				
		//	remove nil from stack
			lua_pop(L, 1);
		}
		lua_pop(L, 1);
		
	//	the requested key has not been in the metatable.
	//	we thus have to check the parents metatable.

		lua_pushstring(L, "class_name_node");
		lua_rawget(L, -2);
		const ClassNameNode* pClassNameNode = (const ClassNameNode*) lua_touserdata(L, -1);
		lua_pop(L, 2);

	//	push all base classes to queue
		for(size_t i = 0; i < pClassNameNode->num_base_classes(); ++i)
		{
			qClassNameNodes.push(&pClassNameNode->base_class(i));
		}

	//	check if some base class left; if not, return nil
		if(qClassNameNodes.empty())
		{
			lua_pushnil(L);
			return 1;
		}

	//	get metatable of first base class
		luaL_getmetatable(L, qClassNameNodes.front()->name().c_str());
		qClassNameNodes.pop();
	}

//	this point should never be reached
	return 0;
}

static int LuaProxyRelease(lua_State* L)
{
	void* ptr = lua_touserdata(L, 1);
//	we only proceed if the userdata encapsulates a smart pointer
	if(((UserDataWrapper*)ptr)->is_smart_ptr()){
	//	invalidate the associated smart-pointer
		if(((UserDataWrapper*)ptr)->is_const())
			((ConstSmartUserDataWrapper*)ptr)->smartPtr.invalidate();
		else
			((SmartUserDataWrapper*)ptr)->smartPtr.invalidate();
	}
	return 0;
}

static int LuaProxyDelete(lua_State* L)
{
	void* ptr = lua_touserdata(L, 1);
//	we perform delete if the user-data is a raw pointer
	if(((UserDataWrapper*)ptr)->is_raw_ptr()){
		RawUserDataWrapper* udata = (RawUserDataWrapper*)ptr;
		if(udata->deleteFunc){
			udata->deleteFunc(udata->obj);
			((RawUserDataWrapper*)ptr)->obj = NULL;
		}
		else{
			UG_LOG("WARNING in LuaProxyDelete: Can't delete object, since"
					"object was not created from script.\n");
		}
	}
	else if(((UserDataWrapper*)ptr)->is_smart_ptr()){
	//	invalidate the associated smart-pointer
		if(((UserDataWrapper*)ptr)->is_const())
			((ConstSmartUserDataWrapper*)ptr)->smartPtr.invalidate();
		else
			((SmartUserDataWrapper*)ptr)->smartPtr.invalidate();
	}
	return 0;
}

bool CreateBindings_LUA(lua_State* L, Registry& reg)
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
//	create some common methods
	lua_getglobal(L, "delete");
	if(!lua_isnil(L, -1)){
	//	the method already exists. Don't recreate it.
		lua_pop(L, 1);
	}
	else{
	//	the method doesn't exist yet. Create it.
		lua_pop(L, 1);
		lua_pushcfunction(L, LuaProxyDelete);
		lua_setglobal(L, "delete");
	}

////////////////////////////////
//	create a method for each global function, which calls LuaProxyFunction with
//	a an index to the function.
	size_t numFuncs = reg.num_functions();
	
	for(size_t i = 0; i < numFuncs; ++i){
		ExportedFunctionGroup* funcGrp = &reg.get_function_group(i);

	//	check whether the function already exists
		lua_getglobal(L, funcGrp->name().c_str());
		if(!lua_isnil(L, -1)){
		//	the method already exists. Don't recreate it.
			lua_pop(L, 1);
			continue;
		}
		lua_pop(L, 1);
		
	//	the function is new. Register it.
		lua_pushlightuserdata(L, funcGrp);
		lua_pushcclosure(L, LuaProxyFunction, 1);
		lua_setglobal(L, funcGrp->name().c_str());
	}


	size_t numClasses = reg.num_classes();
	for(size_t i = 0; i < numClasses; ++i){
		const IExportedClass* c = &reg.get_class(i);
		
	//	check whether the class already exists
		lua_getglobal(L, c->name().c_str());
		if(!lua_isnil(L, -1)){
		//	the class already exists. Don't recreate it.
			lua_pop(L, 1);
			continue;
		}
		lua_pop(L, 1);
		
	//	The class is new. Register it.
		if(c->is_instantiable())
		{
		//	set the constructor-function
			lua_pushlightuserdata(L, (void*)c);
			lua_pushcclosure(L, LuaProxyConstructor, 1);
			lua_setglobal(L, c->name().c_str());
		}

	//	create the meta-table for the object
	//	overwrite index and store the class-name
		luaL_newmetatable(L, c->name().c_str());
	//	we use our custom indexing method to allow method-derivation
		lua_pushcfunction(L, MetatableIndexer);
		lua_setfield(L, -2, "__index");
	//	in order to support smart-pointers we have to overload the garbage-collection
		lua_pushcfunction(L, LuaProxyRelease);
		lua_setfield(L, -2, "__gc");
	//	we have to store the class-names of the class hierarchy
		lua_pushstring(L, "class_name_node");
		lua_pushlightuserdata(L, (void*)&(c->class_name_node()));
		lua_settable(L, -3);

	//	add class name hierarchy
	//	\todo: REMOVE, only needed by info commands.
		lua_pushstring(L, "names");
		lua_pushlightuserdata(L, (void*)c->class_names());
		lua_settable(L, -3);

	//	register methods
	//	NOTE: A C-Closure is registered (a function-pointer with default argument)
		for(size_t j = 0; j < c->num_methods(); ++j){
			const ExportedMethodGroup& m = c->get_method_group(j);
			lua_pushstring(L, m.name().c_str());
			lua_pushlightuserdata(L, (void*)&m);
			lua_pushcclosure(L, LuaProxyMethod, 1);
			lua_settable(L, -3);
		}

	//	register const-methods
	//	NOTE: A C-Closure is registered (a function-pointer with default argument)
		if(c->num_const_methods() > 0)
		{
		//	create a new table in the meta-table and store it in the entry __const
			lua_newtable(L);
			for(size_t j = 0; j < c->num_const_methods(); ++j){
				const ExportedMethodGroup& m = c->get_const_method_group(j);
				lua_pushstring(L, m.name().c_str());
				lua_pushlightuserdata(L, (void*)&m);
				lua_pushcclosure(L, LuaProxyMethod, 1);
				lua_settable(L, -3);
			}
			lua_setfield(L, -2, "__const");
		}

	//	register method-groups
	//	NOTE: Pointers to ExportedMethodGroups are registered (lua userdatas)
		if(c->num_methods() > 0){
			lua_newtable(L);
			for(size_t j = 0; j < c->num_methods(); ++j){
				const ExportedMethodGroup& m = c->get_method_group(j);
				lua_pushstring(L, m.name().c_str());
				lua_pushlightuserdata(L, (void*)&m);
				lua_settable(L, -3);
			}
			lua_setfield(L, -2, "__method_grps");
		}

	//	register const-methods-groups
	//	NOTE: Pointers to ExportedMethodGroups are registered (lua userdatas)
		if(c->num_const_methods() > 0)
		{
		//	create a new table in the meta-table and store it in the entry __const
			lua_newtable(L);
			for(size_t j = 0; j < c->num_const_methods(); ++j){
				const ExportedMethodGroup& m = c->get_const_method_group(j);
				lua_pushstring(L, m.name().c_str());
				lua_pushlightuserdata(L, (void*)&m);
				lua_settable(L, -3);
			}
			lua_setfield(L, -2, "__const_method_grps");
		}
	//	pop the metatable from the stack.
		lua_pop(L, 1);
	}
	
//	now register constructors for class-groups
	size_t numClassGroups = reg.num_class_groups();
	for(size_t i = 0; i < numClassGroups; ++i){
		const ClassGroupDesc* cg = reg.get_class_group(i);

	//	check whether the class-group already exists
		lua_getglobal(L, cg->name().c_str());
		if(!lua_isnil(L, -1)){
		//	the class-group already exists. Don't recreate it.
			lua_pop(L, 1);
			continue;
		}
		lua_pop(L, 1);

	//	The class-group is new. Register the proxy-constructor.
		lua_pushlightuserdata(L, (void*)cg);
		lua_pushcclosure(L, LuaProxyGroupConstructor, 1);
		lua_setglobal(L, cg->name().c_str());
	}

	return true;
}

}//	end of namespace
}//	end of namespace
}//	end of namespace
