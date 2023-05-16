

#ifndef __H__UG__LIB_DISC__TIME_DISC__TIME_INTEGRATOR_OBSERVERS__PYTHON_CALLBACK_OBSERVER
#define __H__UG__LIB_DISC__TIME_DISC__TIME_INTEGRATOR_OBSERVERS__PYTHON_CALLBACK_OBSERVER



#ifdef UG_USE_PYBIND11
#include <Python.h>
#include <pybind11/pybind11.h>
#include <pybind11/functional.h>

namespace py = pybind11;

#include "time_integrator_observer_interface.h"

namespace ug {


//! This class implements a python-based observer for time integration.
/*! (Note: Could be implemented as a pure python class as well ...) */
template<class TDomain, class TAlgebra>
class PythonCallbackObserver
: public ITimeIntegratorObserver<TDomain, TAlgebra>
{
public:
	typedef ITimeIntegratorObserver<TDomain, TAlgebra> base_type;
	typedef GridFunction<TDomain, TAlgebra> grid_function_type;
	typedef py::object TFunctionHandle;

	PythonCallbackObserver() : m_callback()
	{}

	PythonCallbackObserver(TFunctionHandle handle) : m_callback(handle)
	{}


	virtual ~PythonCallbackObserver()
	{}

	virtual bool step_process(SmartPtr<grid_function_type> uNew, int step, number time, number dt)
	{
		//	Call python function
		py::object result_py = m_callback(py::cast(uNew), step, time, dt);
		return result_py.cast<bool>();
	}

	void set_callback(TFunctionHandle callback)
	{ m_callback = callback; }


protected:
	TFunctionHandle m_callback;

};

}

#endif // UG_USE_PYBIND11

#endif // __H__UG__LIB_DISC__TIME_DISC__TIME_INTEGRATOR_OBSERVERS__PYTHON_CALLBACK_OBSERVER
