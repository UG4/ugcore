//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m09 d16

#ifndef __H__UG__INTERFACE__BINDINGS_LUA__
#define __H__UG__INTERFACE__BINDINGS_LUA__

#include <vector>
#include <string>

extern "C" {
#include "externals/lua/lua.h"
#include "externals/lua/lauxlib.h"
#include "externals/lua/lualib.h"
}


#include "common/common.h"
#include "registry/registry.h"

namespace ug
{
namespace bridge
{
namespace lua
{

extern const bool IMLPICIT_SMART_PTR_TO_PTR_CONVERSION;

enum UserDataWrapperTypes{
	RAW_POINTER = 1,
	SMART_POINTER = 1 << 1,
	IS_CONST = 1 << 2
};

struct UserDataWrapper{
	byte type;

	bool is_const()		{return (type & IS_CONST) == IS_CONST;}
	bool is_raw_ptr()	{return (type & RAW_POINTER) == RAW_POINTER;}
	bool is_smart_ptr()	{return (type & SMART_POINTER) == SMART_POINTER;}
};

struct SmartUserDataWrapper : public UserDataWrapper
{
	SmartPtr<void>	smartPtr;
};

struct ConstSmartUserDataWrapper : public UserDataWrapper
{
	ConstSmartPtr<void>	smartPtr;
};

struct RawUserDataWrapper : public UserDataWrapper
{
	void*	obj;
	void (*deleteFunc)(const void*);
};


///	creates bindings for ug_interface and a given lua-state.
/**	If you use ug::script, this method will be invoked automatically.*/
bool CreateBindings_LUA(lua_State* L, Registry& reg);

SmartUserDataWrapper* CreateNewUserData(lua_State* L, const SmartPtr<void>& ptr,
										const char* metatableName);

ConstSmartUserDataWrapper* CreateNewUserData(lua_State* L, const ConstSmartPtr<void>& ptr,
											 const char* metatableName);

RawUserDataWrapper* CreateNewUserData(lua_State* L, void* ptr,
									  const char* metatableName,
									  void (*deleteFunc)(const void*),
									  bool is_const);


}//	end of namespace
}//	end of namespace
}//	end of namespace

#endif
