

#ifndef __H__UG__LIB_DISC__TIME_DISC__TIME_INTEGRATOR_OBSERVERS__PYTHON_CALLBACK_OBSERVER
#define __H__UG__LIB_DISC__TIME_DISC__TIME_INTEGRATOR_OBSERVERS__PYTHON_CALLBACK_OBSERVER



#ifdef UG_USE_PYBIND11
#include <Python.h>
#include <pybind11/pybind11.h>
#include <pybind11/functional.h>

namespace py = pybind11;

#include "time_integrator_observer_interface.h"

namespace ug {


template<class TDomain, class TAlgebra>
class PythonCallbackObserver
: public ITimeIntegratorObserver<TDomain, TAlgebra>
{
public:
	typedef ITimeIntegratorObserver<TDomain, TAlgebra> base_type;
	typedef GridFunction<TDomain, TAlgebra> grid_function_type;
	typedef py::object TFunctionHandle;

	PythonCallbackObserver() : m_callback(), m_u(SPNULL)
	{}

	PythonCallbackObserver(TFunctionHandle handle) : m_callback(handle), m_u(SPNULL)
	{}


	virtual ~PythonCallbackObserver()
	{}

	virtual bool step_process(SmartPtr<grid_function_type> uNew, int step, number time, number dt)
	{
		// store value.
		m_u = uNew;

		//	call python function
		// auto pyArgs = py::tuple(, py::cast(step), py::cast(time), py::cast(dt));
		py::object result_py = m_callback(py::cast(uNew), step, time, dt);
		return result_py.cast<bool>();
	}

	void set_callback(TFunctionHandle callback)
	{ m_callback = callback; }

	SmartPtr<grid_function_type> get_current_solution()
	{ return m_u; }

protected:
	TFunctionHandle m_callback;
	SmartPtr<grid_function_type> m_u;
};

}

#endif // UG_USE_PYBIND11

#endif // __H__UG__LIB_DISC__TIME_DISC__TIME_INTEGRATOR_OBSERVERS__PYTHON_CALLBACK_OBSERVER
