

#ifndef __H__UG__LIB_DISC__TIME_DISC__TIME_INTEGRATOR_OBSERVERS__LUA_CALLBACK_OBSERVER
#define __H__UG__LIB_DISC__TIME_DISC__TIME_INTEGRATOR_OBSERVERS__LUA_CALLBACK_OBSERVER

#ifdef UG_FOR_LUA

#include "time_integrator_observer_interface.h"

namespace ug {


template<class TDomain, class TAlgebra>
class LuaCallbackObserver
: public ITimeIntegratorObserver<TDomain, TAlgebra>
{
public:
	typedef ITimeIntegratorObserver<TDomain, TAlgebra> base_type;
	typedef GridFunction<TDomain, TAlgebra> grid_function_type;
	typedef LuaFunction<number, number> lua_function_type;

	LuaCallbackObserver()
	: m_lua_callback(SPNULL), m_lua_id(0) {}

	LuaCallbackObserver(int lua_id)
	: m_lua_callback(SPNULL), m_lua_id(lua_id) {}

	virtual ~LuaCallbackObserver()
	{}

	virtual bool step_process(SmartPtr<grid_function_type> uNew, int step, number time, number dt)
	{
		if (!m_lua_callback.valid())
			return true;

		number lua_return_value;
		m_u = uNew;
		(*m_lua_callback)(lua_return_value, 4, (number) step, time, dt, (number) m_lua_id);

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


#endif
