/*
 * lua_parsing.h
 *
 */

#ifndef LUA_PARSING_H_
#define LUA_PARSING_H_
#include "externals/lua/lua.h"



namespace ug{
namespace bridge{
namespace lua{
extern const bool IMLPICIT_SMART_PTR_TO_PTR_CONVERSION;
static SmartUserDataWrapper* CreateNewUserData(lua_State* L, const SmartPtr<void>& ptr,
											  const char* metatableName);

static ConstSmartUserDataWrapper* CreateNewUserData(lua_State* L, const ConstSmartPtr<void>& ptr,
											  const char* metatableName);

static RawUserDataWrapper* CreateNewUserData(lua_State* L, void* ptr,
										  const char* metatableName,
										  void (*deleteFunc)(const void*),
										  bool is_const);
template <typename T>
struct LuaParsing;

template <>
struct LuaParsing<bool>{
	static bool check(lua_State* L, int index){
		return lua_isboolean(L, index);
	}
	static bool get(lua_State* L, int index){
		return lua_toboolean(L, index);
	}
	static void push(lua_State* L, bool data){
		lua_pushboolean(L, (data ? 1 : 0));
	}
};

template <>
struct LuaParsing<int>{
	static bool check(lua_State* L, int index){
		return lua_isnumber(L, index);
	}
	static int get(lua_State* L, int index){
		return lua_tointeger(L, index);
	}
	static void push(lua_State* L, int data){
		lua_pushnumber(L, data);
	}
};

template <>
struct LuaParsing<size_t>{
	static bool check(lua_State* L, int index){
		return lua_isnumber(L, index);
	}
	static size_t get(lua_State* L, int index){
		return lua_tointeger(L, index);
	}
	static void push(lua_State* L, size_t data){
		lua_pushnumber(L, data);
	}
};

template <>
struct LuaParsing<float>{
	static bool check(lua_State* L, int index){
		return lua_isnumber(L, index);
	}
	static float get(lua_State* L, int index){
		return lua_tonumber(L, index);
	}
	static void push(lua_State* L, float data){
		lua_pushnumber(L, data);
	}
};

template <>
struct LuaParsing<double>{
	static bool check(lua_State* L, int index){
		return lua_isnumber(L, index);
	}
	static double get(lua_State* L, int index){
		return lua_tonumber(L, index);
	}
	static void push(lua_State* L, double data){
		lua_pushnumber(L, data);
	}
};

template <>
struct LuaParsing<const char*>{
	static bool check(lua_State* L, int index){
		return lua_isstring(L, index);
	}
	static const char* get(lua_State* L, int index){
		return lua_tostring(L, index);
	}
	static void push(lua_State* L, const char* data){
		lua_pushstring(L, data);
	}
};

template <>
struct LuaParsing<std::string>{
	static bool check(lua_State* L, int index){
		return lua_isstring(L, index);
	}
	static std::string get(lua_State* L, int index){
		return std::string(lua_tostring(L, index));
	}
	static void push(lua_State* L, std::string data){
		lua_pushstring(L, data.c_str());
	}
};

template <>
struct LuaParsing<const std::string&>{
	static bool check(lua_State* L, int index){
		return lua_isstring(L, index);
	}
	static void push(lua_State* L, const std::string& data){
		lua_pushstring(L, data.c_str());
	}
};

template <>
struct LuaParsing<void*>{
	static bool checkAndGet(std::pair<void*, const ClassNameNode*>& res,
	                        lua_State* L, int index, const char* baseClassName){
		if(!lua_isuserdata(L, index)) return false;

		UserDataWrapper* udata =
			reinterpret_cast<UserDataWrapper*>(lua_touserdata(L, index));

		if(udata->is_const()) return false;

		//	extract the pointer to the object.
		//	udata is either a RawUserData or a SmartUserDataWrapper
		void* obj = NULL;
		if(udata->is_raw_ptr())
			obj = static_cast<RawUserDataWrapper*>(udata)->obj;
		else if(udata->is_smart_ptr() && IMLPICIT_SMART_PTR_TO_PTR_CONVERSION)
			obj = static_cast<SmartUserDataWrapper*>(udata)->smartPtr.get();
		else return false;

		//	get the object and its metatable. Make sure that obj can be cast to
		//	the type that is expected by the paramsTemplate.
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

	static void push(lua_State* L, void* data, const char* className){
		CreateNewUserData(L, data, className, NULL, false);
	}
};

template <>
struct LuaParsing<const void*>{
	static bool checkAndGet(std::pair<const void*, const ClassNameNode*>& res,
	                        lua_State* L, int index, const char* baseClassName){
		if(!lua_isuserdata(L, index)) return false;

		UserDataWrapper* udata =
			reinterpret_cast<UserDataWrapper*>(lua_touserdata(L, index));

		//	extract the pointer to the object.
		//	udata is either a RawUserData or a SmartUserDataWrapper
		const void* obj = NULL;

		if(udata->is_raw_ptr())
			obj = static_cast<RawUserDataWrapper*>(udata)->obj;
		else if(udata->is_smart_ptr() && IMLPICIT_SMART_PTR_TO_PTR_CONVERSION){
			//	we have to distinguish between const and non-const smart pointers.
			if(udata->is_const())
				obj = static_cast<ConstSmartUserDataWrapper*>(udata)->smartPtr.get();
			else
				obj = static_cast<SmartUserDataWrapper*>(udata)->smartPtr.get();
		}
		else return false;

		//	get the object and its metatable. Make sure that obj can be cast to
		//	the type that is expected by the paramsTemplate.
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

	static void push(lua_State* L, const void* data, const char* className){
	//	we're removing const with a cast. However, it was made sure that
	//	obj is treated as a const value.
		CreateNewUserData(L, (void*)data, className, NULL, true);
	}
};

template <>
struct LuaParsing<SmartPtr<void> >{
	static bool checkAndGet(std::pair<SmartPtr<void>, const ClassNameNode*>& res,
	                        lua_State* L, int index, const char* baseClassName){
		if(!lua_isuserdata(L, index)) return false;

		UserDataWrapper* udata =
			reinterpret_cast<UserDataWrapper*>(lua_touserdata(L, index));

		if(!udata->is_smart_ptr()) return false;
		if(udata->is_const()) return false;

		SmartPtr<void>& obj = ((SmartUserDataWrapper*)lua_touserdata(L, index))->smartPtr;

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

	static void push(lua_State* L, SmartPtr<void> data, const char* className){
		CreateNewUserData(L, data, className);
	}
};

template <>
struct LuaParsing<ConstSmartPtr<void> >{
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
}
}
#endif /* LUA_PARSING_H_ */

