

#ifndef __H__UG__LIB_DISC__TIME_DISC__TIME_INTEGRATOR_OBSERVERS__LUA_CALLBACK_OBSERVER
#define __H__UG__LIB_DISC__TIME_DISC__TIME_INTEGRATOR_OBSERVERS__LUA_CALLBACK_OBSERVER



#include "time_integrator_observer_interface.h"

namespace ug {


template<typename TDomain, typename TAlgebra>
class LuaCallbackObserver
: public ITimeIntegratorObserver<TDomain, TAlgebra>
{
public:
	using base_type = ITimeIntegratorObserver<TDomain, TAlgebra>;
	using grid_function_type = GridFunction<TDomain, TAlgebra>;
	using lua_function_type = LuaCallbackFunction<number>;

	LuaCallbackObserver()
	: m_lua_callback(nullptr), m_lua_id(0) {}

	explicit LuaCallbackObserver(int lua_id)
	: m_lua_callback(nullptr), m_lua_id(lua_id) {}

	~LuaCallbackObserver() override = default;


	bool step_process(SmartPtr<grid_function_type> uNew, int step, number time, number dt) override
	{
		if (!m_lua_callback.valid())
			return true;

		number lua_return_value;
		m_u = uNew;
		(*m_lua_callback)(lua_return_value, 4, step, time, dt, m_lua_id);

		return lua_return_value == 1;
	}

	void set_callback(const char* luaCallback)
	{
		m_lua_callback = make_sp(new lua_function_type());
		m_lua_callback->set_lua_callback(luaCallback, 4);
	}

	SmartPtr<grid_function_type> get_current_solution()
	{ return m_u; }

protected:
	SmartPtr<lua_function_type> m_lua_callback;
	SmartPtr<grid_function_type> m_u;
	const int m_lua_id;
};

}

#endif
