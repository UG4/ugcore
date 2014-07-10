/**
 * \file lua_serialization.cpp
 *
 * \author Martin Rupp
 *
 * \date 23.06.2014
 *
 * Goethe-Center for Scientific Computing 2014.
 *
 * Serialization stuff for lua.
 * not fully developed.
 * the idea is to call Serialize(obj, serializationObj) for all UserData-Objects, if this Serialize is registered in the registry.
 */

#include <iomanip>
#include "bindings/lua/lua_util.h"
#include "bindings/lua/bindings_lua.h"
#include "bridge/bridge.h"
#include "bridge/util.h"
#include "registry/class_helper.h"
#include "common/util/sort_util.h"
#include "common/util/string_util.h"
#include "lua_stack_check.h"
#include "registry/class_name_provider.h"
#include "common/util/string_util.h"

#include "common/util/stringify.h"
//#include "bridge/misc_bridges/serialization.h"
#include "bindings_lua.h"
//#include "lua_parsing.h"

#include "info_commands.h"

using namespace std;

namespace ug{
size_t SerializerGetLastID();

namespace bridge{

string LUAStringEscape(string s)
{
	s = ReplaceAll(s, "\\", "\\\\");
	s = ReplaceAll(s, "\"", "\\\"");
	s = ReplaceAll(s, "\n", "\\n");
	return s;
}


#if 0
namespace lua{
struct MyLuaParsing
{
	static bool checkAndGet(std::pair<ConstSmartPtr<void>, const ClassNameNode*>& res,
	                        lua_State* L, int index, const char* baseClassName){
		if(!lua_isuserdata(L, index)) return false;

		UserDataWrapper* udata =
			reinterpret_cast<UserDataWrapper*>(lua_touserdata(L, index));

		if(!udata->is_smart_ptr()) return false;

		ConstSmartPtr<void> obj;
		if(((UserDataWrapper*)lua_touserdata(L, index))->is_const())
			obj = ((ConstSmartUserDataWrapper*)lua_touserdata(L, index))->smartPtr;
		else
			obj = ((SmartUserDataWrapper*)lua_touserdata(L, index))->smartPtr;

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

	static void push(lua_State* L, ConstSmartPtr<void> data, const char* className){
		CreateNewUserData(L, data, className);
	}
};
}


template <typename T>
static bool PushLuaStackPointerEntryToParamStack(ParameterStack& ps, lua_State* L,
                                                 int index, const char* baseClassName,
                                                 bool bIsVector)
{
	typedef std::pair<T, const ClassNameNode*> result_type;

	result_type res;
	if(!bIsVector){

		UG_LOG(GetLuaTypeString(L, -1) << " " << baseClassName << "\n");
		if(lua::MyLuaParsing::checkAndGet(res, L, index, baseClassName)){
			ps.push(res.first, res.second);
		}
		else return false;
	}
	else {UG_THROW("ll");
	}
	return true;
}

const void *getPtr(lua::UserDataWrapper *self)
{
	//	raw pointer
	if(self->is_raw_ptr())
	{
	//	cast to the needed base class
		return ((lua::RawUserDataWrapper*)self)->obj;
	}
//	smart pointer
	else if(self->is_smart_ptr())
	{
		if(self->is_const())
		//	cast to the needed base class
			return (void*)((lua::ConstSmartUserDataWrapper*)self)->smartPtr.get();
		else
			return ((lua::SmartUserDataWrapper*)self)->smartPtr.get();
	}
}

TheSerializer theSerializer;
#endif

void WriteLUAObject2(lua_State* L, std::string name, std::ostream &s)
{
	//UG_LOG(" " << lua_tostring(L, -1) << "\n");
	if(lua_isnil(L, -1))
		s << name << " = nil\n";
	else if(lua_isboolean(L, -1))
		s << name << " = " << (lua_toboolean(L, -1)==true ? "true" : "false") << "\n";
	else if(lua_isnumber(L, -1))
		s << name << " = " << lua_tostring(L, -1) << "\n";
	else if(lua_isstring(L, -1))
		s << name << " = \"" << LUAStringEscape(lua_tostring(L, -1)) << "\"\n";
	else if(lua_istable(L, -1))
	{
		lua_pushvalue(L, -1);
		/* table is in the stack at index 't' */
		s << name << " = " << name << " or {}\n";
		//UG_LOG(name << " = " << name << " or {}\n");
		lua_pushnil(L);
		while (lua_next(L, -2) != 0)
		{
			std::string nextname;
			// copy the key so that lua_tostring does not modify the original
			lua_pushvalue(L, -2);
			if(lua_type(L, -1) == LUA_TNUMBER)
				nextname = Stringify() << "\t" << name << "[" << lua_tostring(L, -1) << "]";
			else
				nextname = Stringify() << "\t" << name << "[\"" << lua_tostring(L, -1) << "\"]";
			lua_pop(L, 1);

			// also copy value
			lua_pushvalue(L, -1);
			WriteLUAObject2(L, nextname, s);
		    lua_pop(L, 2);
		 }
		lua_pop(L, 1);

	}
	else if(lua_isuserdata(L, -1))
	{
#if 0
		const ClassNameNode&  cnn = SerializeObject::class_name_node();
		lua::UserDataWrapper* self = (lua::UserDataWrapper*)lua_touserdata(L, 1);
		size_t id = theSerializer.is_serialized(self);
		if(id == 0)
		{
			const ExportedFunctionGroup* funcGrp = GetUGRegistry().get_exported_function_group("Serialize");
			if(funcGrp)
			{
			//	we have to try each overload!
				for(size_t i = 0; i < funcGrp->num_overloads(); ++i){
					const ExportedFunction* func = funcGrp->get_overload(i);

					ParameterStack paramsIn;
					ParameterStack paramsOut;
					const ParameterInfo &psInfo = func->params_in();

					if(PushLuaStackPointerEntryToParamStack<ConstSmartPtr<void> >
							(paramsIn, L, -1, psInfo.class_name(0), false))
					{
						size_t id = theSerializer.create_id(self);
						SerializeObject so(theSerializer, id);


						paramsIn.push(&so, &cnn);
						func->execute(paramsIn, paramsOut);

						theSerializer.serialize(so.buffer(), id, psInfo.class_name(0));
						break;
					}

				}
			}

		}


		if(id == 0)
			s << name << " = nil -- no Serialize(" << GetLuaTypeString(L, -1) << ")\n";
		else
			s << name << " = SerializerGet(" << id << ")\n";
#else
		s << name << " = nil -- Userdata serialize not yet supported\n";
#endif

	}
	else
	{
		const char *p =lua_tostring(L, -1);
		s << "-- " << name << " = " << GetLuaTypeString(L, -1) << " ";
		if(p) s << p;
		s << "\n";
	}

}
void WriteLUAObject(lua_State* L, std::string name, std::ostream &s)
{
	lua_getglobal(L, name.c_str());
	WriteLUAObject2(L, name, s);
	lua_pop(L, 1);
}

/*void LuaList_writeObjects(const char*filename)
{
	fstream file(filename, ios::out);
	lua_State* L = script::GetDefaultLuaState();
	std::vector<std::string> luaObjects;
	GetLuaGlobal_luaObjects(luaObjects);

	if(luaObjects.empty()) return;
	int maxLength = (*max_element(luaObjects.begin(), luaObjects.end(), IsLonger)).size();
	for(size_t i=0; i<luaObjects.size(); i++)
	{
		if(luaObjects[i].compare("_G") == 0) continue;
		if(luaObjects[i].compare("package") == 0) continue;
		WriteLUAObject(L, luaObjects[i].c_str(), std::cout);
	}
}*/

void LuaWrite(const char *filename, const char* obj)
{
	fstream file(filename, ios::out);
	file << "function LoadTheCheckpoint()\n";
	lua_State* L = script::GetDefaultLuaState();
	WriteLUAObject(L, obj, file);
	file << "return " << obj << "\n";
	file << "end\n";
}

void LuaWriteCout(const char* obj)
{
	lua_State* L = script::GetDefaultLuaState();
	WriteLUAObject(L, obj, std::cout);
}


bool RegisterSerializationCommands(Registry &reg, const char* parentGroup)
{
	stringstream grpSS; grpSS << parentGroup << "/Serialization";
	std::string grp = grpSS.str();

	try
	{
//		reg.add_function("LuaList_writeObjects", &LuaList_writeObjects, grp.c_str());
		reg.add_function("LuaWrite", &LuaWrite, grp.c_str());
		reg.add_function("LuaWriteCout", &LuaWriteCout, grp.c_str());
	}
	UG_REGISTRY_CATCH_THROW(grp);

	return true;
}

}
}
