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
#include "lua_parsing.h"
#ifdef UG_PARALLEL
#include "pcl/pcl_base.h"
#endif
#include "bridge/bridge.h"

#include "common/util/stringify.h"
#include "lua_stack.h"

//#define __UG__BINDINGS_LUA__CATCH_UNKNOWN_EXCEPTIONS__

using namespace std;
using namespace ug::script;
using namespace ug::bridge;

//	a symbol preceding error messages
static const char* errSymb = " % ";

///	throw mechanism for lua related errors.
#define UG_LUA_THROW_EMPTY(luaState)	luaL_error(luaState, "%s", "__UG__LUA__EMPTY__MSG__")
#define UG_LUA_THROW(luaState, msg)		luaL_error(luaState, "\n%s", msg)

#define UG_LUA_BINDINGS_THROW(L) \
			ug_throw_error(); \
			UG_LUA_THROW_EMPTY(L);

#define UG_LUA_BINDINGS_CATCH(msg, args)\
		catch(LuaError& err){ \
			UG_LUA_THROW(L, err.get_msg().c_str()); \
		} \
		catch(SoftAbort& err){ \
			throw err; \
		} \
		catch(UGError& err){ \
			ug::LogAssistant& la = ug::GetLogAssistant();\
			la.set_output_process(la.get_process_rank());\
			la.flush();\
			UG_LOG((Stringify() << "\n" << errSymb << "EXCEPTION on Proc "<<la.get_process_rank()<<":" \
				<< " UGError thrown\n" \
				<< UGErrorTraceback(err) \
				<< errSymb << msg << "\n" \
				<< errSymb <<"   ARGS: ("<< args << ")\n" \
				<< errSymb<<"In FILE: " << GetLuaFileAndLineNumber(L) << "\n" \
				<< errSymb<<"At LINE: '" << GetLuaLine(L) << "'\n").str()); \
			/*PrintLUACallStack();*/\
			la.flush();\
			UG_LUA_BINDINGS_THROW(L)\
		}\
		catch(std::exception& ex)\
		{\
			ug::LogAssistant& la = ug::GetLogAssistant();\
			la.set_output_process(la.get_process_rank());\
			la.flush();\
			UG_LOG((Stringify() << "\n" << errSymb << "EXCEPTION on Proc "<<la.get_process_rank()<<":\n" \
				<< ErrorStringFromStdException(&ex) << "\n" \
				<< errSymb << msg << "\n" \
				<< errSymb <<"   ARGS: ("<< args << ")\n" \
				<< errSymb<<"In FILE: " << GetLuaFileAndLineNumber(L) << "\n" \
				<< errSymb<<"At LINE: '" << GetLuaLine(L) << "'\n").str()); \
			/*PrintLUACallStack();*/\
			la.flush();\
			UG_LUA_BINDINGS_THROW(L);\
		}\
		catch(...)\
		{\
			ug::LogAssistant& la = ug::GetLogAssistant();\
			la.set_output_process(la.get_process_rank());\
			la.flush();\
			UG_LOG((Stringify() << "\n" << errSymb << "EXCEPTION on Proc "<<la.get_process_rank()<<":" \
				<< " Unknown Exception thrown\n" \
				<< errSymb << msg << "\n" \
				<< errSymb <<"With ARGUMENTS: ("<< args << ")\n" \
				<< errSymb<<"In FILE: " << GetLuaFileAndLineNumber(L) << "\n" \
				<< errSymb<<"At LINE: '" << GetLuaLine(L) << "'\n").str()); \
				/*PrintLUACallStack();*/\
				la.flush();\
			UG_LUA_BINDINGS_THROW(L);\
		}


namespace ug{
namespace bridge{
namespace lua{

//	set this variable to true if smart-ptr arguments shall be automatically
//	converted to raw-ptrs where required.
const bool IMLPICIT_SMART_PTR_TO_PTR_CONVERSION = true;


std::string ParameterStackString(ParameterStack &s)
{
	std::stringstream ss;
	for(int i=0; i<s.size(); i++)
	{
		if(i != 0) ss << ", ";
		if(s.is_vector(i)){
			// todo: more specific output for vectors
			ss << "std::vector";
		} else {
			switch(s.type(i))
			{
				case Variant::VT_INVALID: ss << "invalid"; break;
				case Variant::VT_BOOL: ss << "bool " << s.to<bool>(i); break;
				case Variant::VT_INT: ss << "int " << s.to<int>(i); break;
				case Variant::VT_SIZE_T: ss << "size_t " << s.to<size_t>(i); break;
				case Variant::VT_FLOAT: ss << "float " << s.to<float>(i); break;
				case Variant::VT_DOUBLE: ss << "double " << s.to<double>(i); break;
				case Variant::VT_CSTRING: ss << "cstring \"" << s.to<const char*>(i) << "\""; break;
				case Variant::VT_STDSTRING: ss << "std::string \"" << s.to<std::string>(i) << "\""; break;
				case Variant::VT_POINTER: ss << "pointer " << s.to<void*>(i); break;
				case Variant::VT_CONST_POINTER: ss << "ConstPointer " << s.class_name(i) << " " << s.to<const void*>(i); break;
				case Variant::VT_SMART_POINTER: ss << "SmartPtr<" << s.class_name(i) << "> " << s.to<SmartPtr<void> >(i).get(); break;
				case Variant::VT_CONST_SMART_POINTER: ss << "ConstSmartPtr<" << s.class_name(i) << "> " << s.to<ConstSmartPtr<void> >(i).get(); break;
			#ifdef UG_FOR_LUA
				case Variant::VT_LUA_FUNCTION_HANDLE: ss << "LuaFunctionHandle " << s.to<LuaFunctionHandle>(i).ref; break;
			#endif
				default: ss << "unknown"; break;
			}
		}
	}
	return ss.str();
}

///	creates a new UserDataWrapper and associates it with ptr in luas registry
/**
 * Creates a new userdata in lua, which encapsulates the given pointer / smart-pointer.
 * It then assigns the specified metatable to the userdata.
 * When the function is done, the userdata is left on luas stack.
 * \{
 */
SmartUserDataWrapper* CreateNewUserData(lua_State* L, const SmartPtr<void>& ptr,
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

ConstSmartUserDataWrapper* CreateNewUserData(lua_State* L, const ConstSmartPtr<void>& ptr,
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

RawUserDataWrapper* CreateNewUserData(lua_State* L, void* ptr,
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
 * if one of the parameters is nil, this returns an error
 * to warn the user about not initialised variables
 * @param L						the lua state
 * @param offsetToFirstParam	offset to first param
 * @return error string or ""
 */
string GetNilWarning(lua_State* L, int offsetToFirstParam)
{
	int index = offsetToFirstParam + 1; // stack-indices start with 1
	std::vector<int> indices;
	for(; lua_type(L, index) != LUA_TNONE; index++)
		if(lua_isnil(L, index))
			indices.push_back(index-offsetToFirstParam);

	if(indices.size() > 0)
	{
		std::stringstream ss;

		if(indices.size() == 1)
			ss << errSymb << "WARNING: Argument " << indices[0] << " is nil, it is ";
		else
		{
			ss << errSymb << "WARNING: Arguments ";
			for(size_t i=0; i<indices.size()-1; i++)
			{
				if(i>0) ss << ", ";
				ss << indices[i];
			}
			ss << " and " << indices[indices.size()-1] << " are nil, they are ";
		}
		ss << "possibly not initialized or initialized with an error.\n";
		return ss.str();
	}
	else return "";
}


/**
 * users sometimes use the . operator (like in myObject.the_method(parameter) ), which sometimes print very weird error messages.
 * this function tries to catch some of those errors.
 * NOTE: this function can NOT catch errors like approxSpace.print_statistic() ( note the wrong . )
 * @param name the name of the function. will be used to search "."+name in the code.
 * @return error string or ""
 */
std::string GetColonWarning(std::string name)
{
	std::string line = GetLuaLine(GetDefaultLuaState());
	if(line.find(std::string(".") + name) != std::string::npos)
		return std::string(errSymb) + "WARNING: You seem to have used the point operator. Remember that for calling method b of an object A, you must use A:b(c), NOT A.b(c)\n";
	else
		return "";
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
static string GetTypeMismatchString(const ParameterInfo& paramsTemplate,
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




/**
 * @param err
 * @return traceback of errors in UGError err
 */
string UGErrorTraceback(UGError &err)
{
	std::stringstream ss;
//	header
	ss << errSymb<<"  Error traceback (innermost first): \n";

//	padding to insert
	std::string pad(errSymb); pad.append("     ");

//	print each message
	for(size_t i=0;i<err.num_msg();++i)
	{
	//	get copy of original string
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
		ss << errSymb<<std::setw(3)<<i<<": "<<msg<<endl;

	//	write file and line
		ss << pad << "[at "<<SnipStringFront(err.get_file(i), 62, 3);
		ss<<":"<<err.get_line(i)<<"]"<<endl;
	}
	return ss.str();
}

/**
 * prints the LUA Call Stack, i.e.
 * the call stack of functions calls in LUA up to this point,
 * files and filelines.
 */
void PrintLUACallStack()
{
	std::string s = LuaStackTraceString();
	if(s.length())
	{
		UG_LOG("Call Stack:\n" << s << "\n");
	}
}

/**
 * LuaProxyFunction handling calls to global functions.
 * Note that not the best matching, but the first matching overload is chosen!
 * @param L
 * @return The number of items pushed to the stack
 */
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
		UG_LUA_BINDINGS_CATCH("In CALL to function '" << FunctionInfo(*func) << "'", ParameterStackString(paramsIn));

	//	if we reach this point, then the method was successfully executed.
		return ParamsToLuaStack(paramsOut, L);
	}

	if(badParam != 0)
	{
		UG_LOG(errSymb<<"Error at "<<GetLuaFileAndLine(L) << ":\n");
		UG_LOG(errSymb<<"ERROR occured when trying to call '"
			   << funcGrp->name() << "(" << GetLuaParametersString(L, 0) << "):'\n");
		UG_LOG(GetNilWarning(L, 0));

		UG_LOG(errSymb<<"No matching overload found! Candidates are:\n");
		for(size_t i = 0; i < funcGrp->num_overloads(); ++i)
		{
			const ExportedFunction* func = funcGrp->get_overload(i);
			ParameterStack paramsIn;
			badParam = LuaStackToParams(paramsIn, func->params_in(), L, 0);
			UG_LOG(errSymb<<" - ");
			UG_LOG(FunctionInfo(*func));
			UG_LOG(": " << GetTypeMismatchString(func->params_in(), L, 0, badParam) << "\n");
		}

		PrintLUACallStack();

		UG_LUA_THROW_EMPTY(L);
	}

//	this point shouldn't be reached
	UG_LUA_THROW(L, "Unknown internal error!");
	return 0;
}

/**
 * Helper function of LuaProxyConstructor and LuaProxyGroupConstructor
 * @param L
 * @param c the class to create and object of
 * @param groupname if not nil, c is the default class of this group
 * @return The number of items pushed to the stack (should be one = 1 object).
 */
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
		UG_LUA_BINDINGS_CATCH("When CREATING class '"<< c->name()  << "'", ParameterStackString(paramsIn));

	//	object created
		return 1;
	}

//	no matching overload found
	if(badParam != 0)
	{
		UG_LOG(errSymb<<"Error at "<<GetLuaFileAndLine(L) << ":\n");
		UG_LOG(errSymb<<"ERROR occured when trying to create object of ");
		if(groupname)
		{ 	UG_LOG("group " << groupname << " (default class " << c->name() << ") with constructor '" << c->name());	}
		else
		{ 	UG_LOG("class " << c->name() << " with constructor '" << c->name());	}
		UG_LOG("(" << GetLuaParametersString(L, 0) << ")':\n");
		UG_LOG(GetNilWarning(L, 0));
		UG_LOG(errSymb<<"No matching overload found! Candidates are:\n");
		for(size_t i = 0; i < c->num_constructors(); ++i)
		{
			const ExportedConstructor& constr = c->get_constructor(i);
			ParameterStack paramsIn;
			badParam = LuaStackToParams(paramsIn, constr.params_in(), L, 0);
			UG_LOG(errSymb << " - ");
			UG_LOG(ConstructorInfo(constr, c->name().c_str()));
			UG_LOG(": " << GetTypeMismatchString(constr.params_in(), L, 0, badParam) << "\n");
		}
		PrintLUACallStack();
		//UG_LOG(errSymb<<"Call stack:\n");	LuaStackTrace(L);
		UG_LUA_THROW_EMPTY(L);
	}

	//UG_LOG(errSymb<<"Call stack:\n"); LuaStackTrace(L);
	UG_LUA_THROW_EMPTY(L);

	return 0;
}

/**
 * creates a object of a class
 * @param L
 * @return The number of items pushed to the stack (should be one = 1 object).
 */
static int LuaProxyConstructor(lua_State* L)
{
//	get class
	IExportedClass* c = (IExportedClass*)lua_touserdata(L, lua_upvalueindex(1));
	return LuaConstructor(L, c);
}


/**
 * creates the class which is set as default class for the specified group.
 * we assume that the first upvalue is a ClassGroupDesc*
 * @param L
 * @return The number of items pushed to the stack (should be one = 1 object).
 */
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
		UG_LOG("': No default class set! (Have you called InitUG?)\n");
		lua_pushnil(L);
		return 1;
	}

	return LuaConstructor(L, c, group->name().c_str());
}

/**
 * This method is not called by lua, but a helper to LuaProxyMethod.
 * It recursively calls itself until a matching overload was found.
 * @param L
 * @param methodGrp
 * @param self
 * @param classNameNode
 * @param errorOutput
 * @return The number of items pushed to the stack
 */
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
		UG_LUA_BINDINGS_CATCH("In CALL to method '" << LuaClassMethodInfo(L, 1, *m)  << "'", ParameterStackString(paramsIn));

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
				UG_LOG("- " << FunctionInfo(*func) << ": " << GetTypeMismatchString(func->params_in(), L, 1, badParam) << "\n");
			}
		}
	}
	
	return -1;
}
/**
 * a default __tostring method which shows classname: \<adress\>
 * __tostring is used in all print(object) and tostring(object) calls in LUA
 * @param L
 * @return nr of parameters (its one string)
 */
static int LuaToStringDefault(lua_State *L)
{
	IExportedClass* c = (IExportedClass*)lua_touserdata(L, lua_upvalueindex(1));
	ParameterStack out;
	char buf[255];
	try{
		sprintf(buf, "[ %s: %p", c->name().c_str(), c);

	}
	catch(...)
	{
		sprintf(buf, "%p", c);

	}
	string b = buf;
	out.push(b);
	return ParamsToLuaStack(out, L);
}



/**
 * member methods of classes are handled here
 * @param L the lua State
 * @return number of parameters returned
 */
static int LuaProxyMethod(lua_State* L)
{
	const ExportedMethodGroup* methodGrp = (const ExportedMethodGroup*)
											lua_touserdata(L, lua_upvalueindex(1));

	if(!lua_isuserdata(L, 1))
	{
		UG_LOG(GetColonWarning(methodGrp->name()));
		UG_LOG(errSymb<<"Error at "<<GetLuaFileAndLine(L) << ":\n")
		UG_LOG(errSymb<<"Error in call to LuaProxyMethod: No object specified in call to '");
		UG_LOG(LuaClassMethodInfo(L, 1, *methodGrp->get_overload(0)));
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

	UG_LOG(GetNilWarning(L, 1));
	UG_LOG(GetColonWarning(methodGrp->name()));
	UG_LOG(errSymb<<"No matching overload found! Candidates in class "
			<< classname << " are:\n");

//	since no matching method was found, we'll run ExecuteMethod again, this
//	time outputting the errors
	ExecuteMethod(L, methodGrp, self, classNameNode, true);

	PrintLUACallStack();

	//UG_LOG(errSymb<<"Call stack:\n"); LuaStackTrace(L);
	UG_LUA_THROW_EMPTY(L);

	return 0;
}

/**
 * This function is used if the user enters a member function wrong and returns
 * the best guess what he means.
 * Also searches in base classes.
 * @param c			(in) the class
 * @param name		(in) the name of the member function we are searching
 * @param minname	(out) the name of the closest match to name (in Levenshtein-Distance)
 * @param mind		(out) number of edits needed from name -> minname
 */
void GetBestMatchingMember(const IExportedClass *c, const char *name, std::string &minname, int &mind)
{
	const char *funcname;
	for(size_t i=0; i<c->num_const_methods(); i++)
	{
		funcname = c->get_const_method(i).name().c_str();
		int d = LevenshteinDistance(funcname, name );
		if(d < mind)
		{
			minname = funcname;
			mind = d;
		}
	}
	for(size_t i=0; i<c->num_methods(); i++)
	{
		funcname = c->get_method(i).name().c_str();

		int d = LevenshteinDistance(funcname, name );
		if(d < mind)
		{
			minname = funcname;
			mind = d;
		}
	}
	// search in the classes derived from
	const vector<const char *> *pNames = c->class_names();
	if(pNames)
	{
		bridge::Registry &reg = ug::bridge::GetUGRegistry();
		for(size_t i=1; i<pNames->size(); i++) // skip the class itself
			GetBestMatchingMember(reg.get_class(pNames->at(i)), name, minname, mind);
	}
}

/**
 * @sa GetBestMatchingMember
 */
void GetBestMatchingMember(const char *classname, const char *name, std::string &minname, int &mind)
{
	bridge::Registry &reg = ug::bridge::GetUGRegistry();
	const IExportedClass *c = reg.get_class(classname);
	minname = "";
	mind = 10;
	GetBestMatchingMember(c, name, minname, mind);
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
			lua_getmetatable(L, 1);
			lua_pushstring(L, "class_name_node");
			lua_rawget(L, -2);
			const ClassNameNode* pClassNameNode = (const ClassNameNode*) lua_touserdata(L, -1);
			lua_pop(L, 2);
			const char *classname = pClassNameNode->name().c_str();
			const char *funcname = lua_tostring(L, -1);
			UG_LOG("LUA ERROR! Could not find member function \"" << funcname << "\" of class " << classname << ". ");
			lua_pushnil(L);

			std::string minname;
			int mind =0;
			GetBestMatchingMember(classname, funcname, minname, mind);
			if(mind < (int)strlen(funcname)/2)
			{ UG_LOG("Did you mean " << minname << " ?\n"); }
			else { UG_LOG("\n"); }
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



/**
 * Registers a meta-object for each class and class group found in the ug registry reg.
 * Global functions are registered for all GlobalFunction-objects in the registry reg.
 * @param L the Lua State
 * @param reg the ug registry
 * @return true on success
 */
bool CreateBindings_LUA(lua_State* L, Registry& reg)
{
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
		lua_pushstring(L, "__names");
		lua_pushlightuserdata(L, (void*)c->class_names());
		lua_settable(L, -3);

		bool bToStringFound =false;

	//	register methods
	//	NOTE: A C-Closure is registered (a function-pointer with default argument)
		for(size_t j = 0; j < c->num_methods(); ++j){
			const ExportedMethodGroup& m = c->get_method_group(j);
			lua_pushstring(L, m.name().c_str());
			lua_pushlightuserdata(L, (void*)&m);
			lua_pushcclosure(L, LuaProxyMethod, 1);
			lua_settable(L, -3);
			if(m.name().compare("__tostring") == 0) bToStringFound = true;
		}

		if(bToStringFound==false)
		{
			// we did not find a __tostring function for lua
			// for print(myObject) to use myObject:__tostring() it is necessary
			// that __tostring is in the direct metatable of myObject, and not in __const
			// (somehow lua does not use the metatable indexer provided here)
			// so if there is a __tostring() const function, but none non-const,
			// we also put the __tostring() const function in the non-const part,
			// which should be no harm, since a const function could also
			// be a non-const function (the reverse is not true)
			for(size_t j = 0; j < c->num_const_methods(); ++j)
			{
				const ExportedMethodGroup& m = c->get_const_method_group(j);
				if(m.name().compare("__tostring") == 0)
				{
					lua_pushstring(L, m.name().c_str());
					lua_pushlightuserdata(L, (void*)&m);
					lua_pushcclosure(L, LuaProxyMethod, 1);
					lua_settable(L, -3);
					bToStringFound = true;
					break;
				}
			}
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
		
		// if we did not find a __tostring function, we
		// add a default __tostring method which shows 
		// classname: <adress>
		if(bToStringFound == false)
		{
			lua_pushstring(L, "__tostring");
			lua_pushlightuserdata(L, (void*)c);
			lua_pushcclosure(L, LuaToStringDefault, 1);
			lua_settable(L, -3);
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
