
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
const std::vector<const char*> *GetClassNames(lua_State* L, const char *name);

void PrintLuaClassMethodInfo(lua_State *L, int index, const ExportedMethod &thefunc);

std::string GetFileLine(const char *filename, size_t line);

int PrintFunctionInfo(lua_State *L, bool bComplete);

int UGTypeInfo(const char *p);

bool GetLuaNamespace(lua_State* L, std::string &name);


class LuaStackCheck
{
public:
	LuaStackCheck(lua_State *L, const char *msg, int growth=0)
	{
		m_L = L;
		m_topBefore = lua_gettop(L);
		m_growth = growth;
		m_msg = msg;
	}
	~LuaStackCheck()
	{
		int m_topAfter = lua_gettop(m_L);
		UG_ASSERT(m_topBefore+m_growth == m_topAfter, m_msg << ": " << m_topBefore << " + " << m_growth << " != " << m_topAfter << "\n");
	}
private:
	lua_State *m_L;
	int m_topBefore;
	int m_growth;
	const char *m_msg;
};

#define LUA_STACK_CHECK_STRINGIFY(x) #x
#define LUA_STACK_CHECK_TOSTRING(x) LUA_STACK_CHECK_STRINGIFY(x)

#define LUA_STACK_CHECK(L, growth) LuaStackCheck check##__LINE__(L, __FILE__ ":" LUA_STACK_CHECK_TOSTRING(__LINE__), growth)


}
}
						  
#endif // __H__UG_SCRIPT__INFO_COMMANDS__
