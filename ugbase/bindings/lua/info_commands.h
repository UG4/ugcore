
#include "registry/registry.h"
#include <string>

#ifndef __H__UG_SCRIPT__INFO_COMMANDS__
#define __H__UG_SCRIPT__INFO_COMMANDS__

#include "common/ug_config.h"

namespace ug
{
namespace bridge
{

///	registers info commands TypeInfo, ClassUsage and others
UG_API bool RegisterInfoCommands(bridge::Registry &reg,
						  const char* parentGroup = "/ug4");

UG_API const ClassNameNode* GetClassNameNode(lua_State *L, int index);

UG_API const std::vector<const char*> *GetClassNames(lua_State *L, int index);
UG_API const std::vector<const char*> *GetClassNames(lua_State* L, const char *name);

UG_API void PrintLuaClassMethodInfo(lua_State *L, int index, const ExportedMethod &thefunc);

UG_API std::string GetFileLine(const char *filename, size_t line);

UG_API int PrintFunctionInfo(lua_State *L, bool bComplete);

UG_API int UGTypeInfo(const char *p);

UG_API bool GetLuaNamespace(lua_State* L, std::string name);

UG_API double LuaGetNumber(lua_State *L, const char *name, double notAvailable);
UG_API std::string LuaGetString(lua_State *L, const char *name, const char *notAvailable);
UG_API bool LuaGetBoolean(lua_State *L, const char *name, bool notAvailable);

UG_API void lua_printCurrentLine(lua_State* L);
UG_API void lua_getLastLine(lua_State* L, lua_Debug entry);


/// returns a String describing the content of the lua stack at a given index
std::string GetLuaTypeString(lua_State* L, int index);

/// prints information about lua's call stack (file:line source).
void lua_stacktrace(lua_State* L);

/// returns the current file and line ( \sa lua_stacktrace ).
std::string GetLuaFileAndLine(lua_State* L);

}
}
						  
#endif // __H__UG_SCRIPT__INFO_COMMANDS__
