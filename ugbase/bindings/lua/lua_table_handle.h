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

namespace ug{
class Variant;
namespace impl{
struct LuaTableHandle_;
}

/// Handle for a lua reference
class LuaTableHandle /* public SuitableBaseClass */ {
public:
	class const_iter /* erase in SuitableBaseClass */ {
		typedef std::pair<std::string, ug::Variant> value_type;
	public:
		const_iter();
		~const_iter();
		const_iter(impl::LuaTableHandle_ /*const*/ * _data);
		const_iter(const_iter const&);

	public:
		value_type operator*() const;
		const_iter& operator++();
		const_iter& operator=(const_iter const& o);
		bool operator==(const_iter const& o) const;
		bool operator!=(const_iter const& o) const{
			return !operator==(o);
		}

	private:
		impl::LuaTableHandle_ /*const*/ * _data;
		int _keyref;
		int _valref;
		int _cnt;
	}; // const_iter
public:
	LuaTableHandle() : _data(nullptr) {}
	LuaTableHandle(LuaTableHandle const&);
	LuaTableHandle(LuaTableHandle&&);
	explicit LuaTableHandle(lua_State* ref, int idx);
	~LuaTableHandle();

public:
	bool operator==(LuaTableHandle const& o) const;
	size_t size() const;
	size_t count() const;
	ug::Variant get(std::string const& key) const;
	ug::Variant get(int const& key) const;

	const_iter begin();
	const_iter end();

	LuaTableHandle& operator=(LuaTableHandle const&);
	LuaTableHandle& operator=(LuaTableHandle&&);

private:
	impl::LuaTableHandle_* _data;
};

} // ug

#endif // guard
