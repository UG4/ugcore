/*
 * Copyright (c) 2011:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef LUA_STACK_CHECK_H_
#define LUA_STACK_CHECK_H_

namespace ug {
namespace bridge {

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
#endif