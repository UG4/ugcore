
#include "ug_bridge/registry.h"

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
void PrintFunctionInfo(const bridge::ExportedFunctionBase &thefunc, bool isConst=false, const char *classname=NULL, const char *highlightclassname=NULL);
void PrintLuaClassMethodInfo(lua_State *L, int index, const ExportedMethod &thefunc);
}
}
						  
#endif // __H__UG_SCRIPT__INFO_COMMANDS__
