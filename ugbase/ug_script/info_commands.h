
#include "ug_bridge/registry.h"
#include <string>

#ifndef __H__UG_SCRIPT__INFO_COMMANDS__
#define __H__UG_SCRIPT__INFO_COMMANDS__

namespace ug
{
namespace bridge
{

///	registers info commands TypeInfo, ClassUsage and others
bool RegisterInfoCommands(bridge::Registry &reg,
						  const char* parentGroup = "/ug4");

const ClassNameNode* GetClassNameNode(lua_State *L, int index);

const std::vector<const char*> *GetClassNames(lua_State *L, int index);
const std::vector<const char*> *GetClassNames(lua_State* L, const char *name);

void PrintLuaClassMethodInfo(lua_State *L, int index, const ExportedMethod &thefunc);

std::string GetFileLine(const char *filename, size_t line);

int PrintFunctionInfo(lua_State *L, bool bComplete);

int UGTypeInfo(const char *p);

bool GetLuaNamespace(lua_State* L, std::string name);

int LuaGetNumber(lua_State *L, const char *name, int notAvailable);
std::string LuaGetString(lua_State *L, const char *name, const char *notAvailable);
bool LuaGetBoolean(lua_State *L, const char *name, bool notAvailable);


}
}
						  
#endif // __H__UG_SCRIPT__INFO_COMMANDS__
