/*
 * lua_stack_check.h
 *
 *  Created on: 29.03.2011
 *      Author: mrupp
 */

#ifndef LUA_STACK_CHECK_H_
#define LUA_STACK_CHECK_H_

namespace ug
{
namespace bridge
{

#ifndef NDEBUG
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

#define LUA_STACK_CHECK(L, growth) ::ug::bridge::LuaStackCheck check##__LINE__(L, __FILE__ ":" LUA_STACK_CHECK_TOSTRING(__LINE__), growth)

#else

#define LUA_STACK_CHECK(L, growth)

#endif


}
}
#endif /* LUA_STACK_CHECK_H_ */
