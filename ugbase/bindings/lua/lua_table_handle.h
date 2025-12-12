/*
 * Copyright (c) 2021  G-CSC, Goethe University Frankfurt
 * Author: Felix Salfelder
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

#ifndef H__UG__LUA_TABLE_HANDLE__
#define H__UG__LUA_TABLE_HANDLE__

#include <string>
#include <cstddef>

using std::size_t;

typedef struct lua_State lua_State;

namespace ug {

class Variant;

namespace impl {
struct LuaTableHandle_;
}

/// Handle for a lua reference
class LuaTableHandle /* public SuitableBaseClass */ {
public:
	LuaTableHandle() = delete;
	LuaTableHandle(const LuaTableHandle&);
	LuaTableHandle(LuaTableHandle&&) noexcept;
	explicit LuaTableHandle(lua_State* ref, int idx);
	~LuaTableHandle();

public:
	[[nodiscard]] size_t size() const;
	[[nodiscard]] Variant get(std::string const& key) const;
	[[nodiscard]] Variant get(int const& key) const;

	LuaTableHandle& operator = (const LuaTableHandle&);
	LuaTableHandle& operator = (LuaTableHandle&&);

private:
	impl::LuaTableHandle_* _data;
};

} // ug

#endif
