
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


const std::vector<const char*> *GetClassNames(lua_State *L, int index);
void PrintLuaClassMethodInfo(lua_State *L, int index, const ExportedMethod &thefunc);

std::string GetFileLine(const char *filename, size_t line);
}
}
						  
#endif // __H__UG_SCRIPT__INFO_COMMANDS__
