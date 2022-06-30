// LICENSE HEADER

#include "lua_table_handle.h"
#include <iostream>
#include <assert.h>
#include <sstream>
#include <iomanip>
#include "common/util/variant.h"
#include "common/util/trace.h"

#define TRACE_UNTESTED
#include "common/util/trace.h"

extern "C" {
#include "externals/lua/lua.h"
#include "externals/lua/lauxlib.h"
#include "externals/lua/lualib.h"
#include "externals/lua/ldo.h"
}

namespace ug {
namespace impl {

static ug::Variant pop2var(lua_State* _L)
{
	int t = lua_type(_L, -1);
	ug::Variant ret;

	if(t == LUA_TTABLE){
		LuaTableHandle h(_L, -1);
		ret = ug::Variant(h);
	}else if(t == LUA_TNUMBER){
		number n = lua_tonumber(_L, -1);
		ret = ug::Variant(n);
	}else if(t == LUA_TBOOLEAN){
		bool b = lua_toboolean(_L, -1);
		ret = ug::Variant(b);
	}else if(lua_isstring(_L, -1)){
		std::string s(lua_tostring(_L, -1));
		ret = ug::Variant(s);
	}else if(t == LUA_TNIL){
	}else{
		std::cerr << "type not handled " << t << "\n";
	}
	lua_pop(_L, 1); // pop value

	return ret;
}

struct LuaTableHandle_{
public:
	LuaTableHandle_(lua_State*, int);
	~LuaTableHandle_();
	lua_State* _L{nullptr};
	int _ref{0};
	int _index{0};

	bool operator==(LuaTableHandle_ const& o) const{
		return _ref == o._ref;
	}
	bool operator!=(LuaTableHandle_ const& o) const{
		return !operator==(o);
	}

	static void attach(LuaTableHandle_* c, LuaTableHandle_** to) {
		assert(to);
		if (c == *to) {
		}else if (!c) { untested();
			detach(to);
		}else if (!*to) {
			++(c->_attach_count);
			*to = c;
		}else if (*c != **to) {
			// They are different, usually by edit.
			detach(to);
			++(c->_attach_count);
			*to = c;
		}else if (c->_attach_count == 0) { untested();
			// The new and old are identical.
			// Use the old one.
			// The new one is not used anywhere, so throw it away.
			delete c;
		}else{ untested();
			// The new and old are identical.
			// Use the old one.
			// The new one is also used somewhere else, so keep it.
		}
	}
	static void detach(LuaTableHandle_** from) {
		assert(from);
		if (*from) {
			assert((**from)._attach_count > 0);
			--((**from)._attach_count);
			if ((**from)._attach_count == 0) { untested();
				delete *from;
			}else{
			}
			*from = NULL;
		}else{
		}
	}

	size_t size() const{
		lua_rawgeti(_L, LUA_REGISTRYINDEX, _ref);
		size_t ret = lua_objlen(_L, -1);
		lua_pop(_L, 1);
		return ret;
	}

	size_t count() const{ untested();
		lua_rawgeti(_L, LUA_REGISTRYINDEX, _ref);
		size_t n = 0;
		lua_pushnil(_L);
		while (lua_next(_L, -2) != 0) { untested();
			lua_pop(_L, 1);
			++n;
		}
		lua_pop(_L, 1);
		return n;
	}

	ug::Variant get(int const& key) const{
		lua_rawgeti(_L, LUA_REGISTRYINDEX, _ref);
		lua_rawgeti(_L, -1, key+1); // lua starts at 1.

		ug::Variant ret = pop2var(_L);

		lua_pop(_L, 1); // pop ref

		return ret;
	}
	ug::Variant get(std::string const& key) const{
		// std::cerr << "LTH::get " << key << " ref " << _ref << "idx" << _index << "\n";
		lua_rawgeti(_L, LUA_REGISTRYINDEX, _ref);
		lua_getfield(_L, -1, key.c_str());

		ug::Variant ret = pop2var(_L);

		lua_pop(_L, 1); // pop ref
		return ret;
	}

//	int get_int() const { }

private:
	int _attach_count{0};
};

LuaTableHandle_::LuaTableHandle_(lua_State* L, int index)
{
	lua_pushvalue(L, index); // copy table to top of stack.
	int ref = luaL_ref(L, LUA_REGISTRYINDEX); // pops copy from stack

	_ref = ref;
	_L = L;
	_index = index;

	_attach_count = 0;
}

LuaTableHandle_::~LuaTableHandle_()
{ untested();
	assert(!_attach_count);

	// drop reference to table from index.
	luaL_unref(_L, LUA_REGISTRYINDEX, _ref);
}

} // impl

LuaTableHandle::LuaTableHandle(lua_State* L, int index)
	: _data(nullptr)
{
	auto d = new impl::LuaTableHandle_(L, index);
	impl::LuaTableHandle_::attach(d, &_data);
	//lua_settable(L, 0); // lapi.c?
	//lua_gettop(L); // lapi.c?
}

LuaTableHandle::LuaTableHandle(LuaTableHandle const& p)
    : _data(nullptr)
{
	impl::LuaTableHandle_::attach(p._data, &_data);
}

LuaTableHandle::LuaTableHandle(LuaTableHandle&& p)
    : _data(nullptr)
{
	impl::LuaTableHandle_::attach(p._data, &_data);
	impl::LuaTableHandle_::detach(&p._data);
}

LuaTableHandle::~LuaTableHandle()
{
	impl::LuaTableHandle_::detach(&_data);
}

bool LuaTableHandle::operator==(LuaTableHandle const& o) const
{
	if(_data == o._data){
		return true;
	}else if(_data && o._data){
		return *_data == *o._data;
	}else{
		return false;
	}
}

LuaTableHandle& LuaTableHandle::operator=(LuaTableHandle const& p)
{
	impl::LuaTableHandle_::attach(p._data, &_data);
	return *this;
}

LuaTableHandle& LuaTableHandle::operator=(LuaTableHandle&& p)
{
	impl::LuaTableHandle_::attach(p._data, &_data);
	impl::LuaTableHandle_::detach(&p._data);
	return *this;
}

size_t LuaTableHandle::count() const
{ untested();
	assert(_data);
	return _data->count();
}

size_t LuaTableHandle::size() const
{
	assert(_data);
	return _data->size();
}

ug::Variant LuaTableHandle::get(std::string const& key) const
{
	assert(_data);
	return _data->get(key);
}

ug::Variant LuaTableHandle::get(int const& key) const
{
	assert(_data);
	return _data->get(key);
}

LuaTableHandle::const_iter LuaTableHandle::begin()
{
	return const_iter(_data);
}

LuaTableHandle::const_iter LuaTableHandle::end()
{
	return const_iter();
}

LuaTableHandle::const_iter::const_iter()
    : _data(nullptr), _keyref(-1), _valref(-1), _cnt(-1)
{
}

LuaTableHandle::const_iter::const_iter(LuaTableHandle::const_iter const& o)
    : _data(nullptr), _keyref(o._keyref), _valref(o._valref), _cnt(o._cnt)
{ untested();
	impl::LuaTableHandle_::attach(o._data, &_data);

	if(!_data){ untested();
	}else if(_cnt!=-1){ untested();
		auto L = _data->_L;
		lua_rawgeti(L, LUA_REGISTRYINDEX, _keyref);
		_keyref = luaL_ref(L, LUA_REGISTRYINDEX);
		lua_rawgeti(L, LUA_REGISTRYINDEX, _valref);
		_valref = luaL_ref(L, LUA_REGISTRYINDEX);
	}else{ untested();
	}
}

LuaTableHandle::const_iter::~const_iter()
{
	if(_cnt==-1){
	}else if(_data){ untested();
		luaL_unref(_data->_L, LUA_REGISTRYINDEX, _keyref);
		luaL_unref(_data->_L, LUA_REGISTRYINDEX, _valref);
		impl::LuaTableHandle_::attach(nullptr, &_data);
	}else{
	}
}

LuaTableHandle::const_iter::const_iter(impl::LuaTableHandle_ /* const */ * d)
    : _data(nullptr), _keyref(-1), _valref(-1), _cnt(-1)
{
	impl::LuaTableHandle_::attach(d, &_data);
	assert(_data);
	auto L = _data->_L;
	lua_rawgeti(L, LUA_REGISTRYINDEX, _data->_ref); // get table.

	lua_pushnil(L);

	if(lua_next(L, -2)) {
		_cnt = 0;
		_valref = luaL_ref(L, LUA_REGISTRYINDEX); // pop value
		_keyref = luaL_ref(L, LUA_REGISTRYINDEX); // pop key
	}else{
	}
}

LuaTableHandle::const_iter&
LuaTableHandle::const_iter::operator=(LuaTableHandle::const_iter const& o)
{ untested();
	_keyref = o._keyref;
	_valref = o._valref;
	_cnt = o._cnt;
	impl::LuaTableHandle_::attach(o._data, &_data);

	if(!_data){ untested();
	}else if(_cnt!=-1){ untested();
		auto L = _data->_L;
		lua_rawgeti(L, LUA_REGISTRYINDEX, _keyref);
		_keyref = luaL_ref(L, LUA_REGISTRYINDEX);
		lua_rawgeti(L, LUA_REGISTRYINDEX, _valref);
		_valref = luaL_ref(L, LUA_REGISTRYINDEX);
	}else{ untested();
	}

	return *this;
}

bool LuaTableHandle::const_iter::operator==(LuaTableHandle::const_iter const& o) const
{
	return _cnt == o._cnt;
}

static std::string lua_like_string(double number)
{
	std::ostringstream oss;
	oss << std::fixed << std::setprecision(2) << number;
	std::string s = oss.str();
	if(s.find('.') != std::string::npos) {
		s = s.substr(0, s.find_last_not_of('0')+1);
		if(s.find('.') == s.size()-1) {
			s = s.substr(0, s.size()-1);
		}else{
		}
	}else{
	}
	return s;
}

LuaTableHandle::const_iter::value_type
LuaTableHandle::const_iter::operator*() const
{
	assert(_data);
	auto L = _data->_L;

	lua_rawgeti(L, LUA_REGISTRYINDEX, _keyref); // get key
	ug::Variant f = ug::impl::pop2var(L); // string?

	// std::cerr << "deref, key " << f << "\n";
	std::string ff;
  	try{
		ff = f.to_std_string();
	}catch(...){
		ff = lua_like_string(f.to_float());
	}

	lua_rawgeti(L, LUA_REGISTRYINDEX, _valref); // get value
	ug::Variant s = ug::impl::pop2var(L);

	return std::make_pair(ff, s);
}

LuaTableHandle::const_iter& LuaTableHandle::const_iter::operator++()
{
	assert(_data);
	auto L = _data->_L;
	// std::cerr << "op++ " << _data->_ref << " " << _keyref << "\n";
	lua_rawgeti(L, LUA_REGISTRYINDEX, _data->_ref);
	lua_rawgeti(L, LUA_REGISTRYINDEX, _keyref); // get key.
	if(lua_next(L, -2)) {
		++_cnt;
		lua_rawseti(L, LUA_REGISTRYINDEX, _valref); // pops from stack
		lua_rawseti(L, LUA_REGISTRYINDEX, _keyref); // pops from stack
	}else{
		_cnt = -1;
		_keyref = -1;
		_valref = -1;
		lua_pop(L, 1);
	}
	lua_pop(L, 1); // pop table
	return *this;
}

} // ug
