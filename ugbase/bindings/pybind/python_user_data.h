/*
 * Copyright (c) 2022:  G-CSC, Goethe University Frankfurt
 * Author: Sam Gimbel
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

#include <iostream>

#include "lib_disc/spatial_disc/user_data/std_glob_pos_data.h"
#include "lib_disc/spatial_disc/user_data/user_function.h"
#include "lib_disc/spatial_disc/user_data/linker/linker.h"

#ifdef UG_USE_PYBIND11
#include <Python.h>
#include <pybind11/pybind11.h>
#include <pybind11/functional.h>

namespace py = pybind11;
namespace ug
{
////////////////////////////////////////////////////////////////////////////////
// Executing Python Code from Ug4
// Code has been passed as String to UG4
// More Info here: https://docs.python.org/3/extending/embedding.html#
// and also here: https://docs.python.org/3/c-api/ 
// also make sure to have python3-dev installed to use the python header:
// sudo apt install python3-dev
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// PythonUserFunction
// Please Note: This is just a quick fix
// the PythonUserFunction class should be extended like the lua equivalent class
////////////////////////////////////////////////////////////////////////////////


template<int dim>
struct PyUserDataTraits
{
	typedef typename std::function<double(double,int)> TCppFunction;
	template <typename TFunction>
	static py::object call(TFunction f, MathVector<dim>, double t, int si)
	{return f(t,si);}
};

template<>
struct PyUserDataTraits<1>
{
	typedef typename std::function<double(double, double, int)> TCppFunction;

	template <typename TFunction>
	static py::object  call(TFunction f, const MathVector<1> &x, double t, int si)
	{return f(x[0],t,si);}
};

template<>
struct PyUserDataTraits<2>
{
	typedef typename std::function<double(double, double, double, int)> TCppFunction;
	template <typename TFunction>
	static py::object call(TFunction f, const MathVector<2> &x, double t, int si)
	{return f(x[0],x[1],t,si);}
};

template<>
struct PyUserDataTraits<3>
{
	typedef typename std::function<double(double, double, double, double, int)> TCppFunction;

	template <typename TFunction>
	static py::object call(TFunction f, const MathVector<3> &x, double t, int si)
	{ return f(x[0],x[1],x[2],t,si); }
};

template <typename TData, typename TRet>
struct PyUserDataTypeTraits
{
	static TRet return_value(py::object result_py)
	{ return py::make_tuple(result_py)[0].cast<TRet>();}

	static TData data_value(py::object result_py)
	{ return py::make_tuple(result_py)[1].cast<TData>();}

};

template <typename TData>
struct PyUserDataTypeTraits<TData, void>
{
	static void return_value(py::object result_py)
	{return void();}

	static TData data_value(py::object result_py)
	{ return result_py.cast<TData>();}

};




//! This object allows to evaluate a function with arguments (x, t, si)
template <typename TData, int dim, typename TRet = void>
class PythonUserData
: public StdGlobPosData<PythonUserData<TData, dim, TRet>, TData, dim, TRet>
{

	public:

	typedef py::object TFunction;

	//!	Constructor
	/*! Creates a PythonUserData that uses a Python function to evaluate some data. */
	PythonUserData(TFunction f) : func(f)
	{
		UG_ASSERT(check_callback_returns(f), "Huhh: Function has invalid signature.")
	}

protected:

	static inline TRet evaluate_func(TFunction f, TData& value, const MathVector<dim>& x, number time, int si)
	{
		py::object result_py = PyUserDataTraits<dim>().call(f, x, time, si);
		value = PyUserDataTypeTraits<TData, TRet>::data_value(result_py);
		return PyUserDataTypeTraits<TData, TRet>::return_value(result_py);
	}

public:
	//!	Evaluates the data at a given point and time.
	inline TRet evaluate(TData& value, const MathVector<dim>& x, number time, int si) const
	{ return evaluate_func(func, value, x, time, si); };

protected:
	//! This function checks, if the Python function has the correct signature.
	/*! Performs dummy call and evaluates output. */
	static bool check_callback_returns(TFunction func)
	{
		TData data; MathVector<dim> x = 0.0; number time = 0.0; int si = 0;

		try {
			evaluate_func(func, data, x, time, si);
		} catch (std::exception &e) {
			return false;
		}
		return true;

	};

	//! Callback function
	TFunction func;

};

////////////////////////////////////////////////////////////////////////////////
// PythonUserFunction
////////////////////////////////////////////////////////////////////////////////

/// maps several data values to an output data value using a lua callback
/**
 * This class provides the evaluation of a user function, that is specified
 * in the script. Several data (of the same c++-type) can be used as input,
 * a data (of the same type) is returned.
 */
template <typename TData, int dim, typename TDataIn>
class PythonUserFunction : public StdDataLinker<PythonUserFunction<TData, dim, TDataIn>, TData, dim>
{
	public:
	//	type of base class
		typedef StdDataLinker<PythonUserFunction<TData, dim, TDataIn>, TData, dim> base_type;
		using base_type::set_input;

		typedef py::object TFunctionHandle;
	/**
	 * \brief constructor
	 * \param pyCallback name of the Lua function to use as callback
	 * \param numArgs number of arguments of the Lua callback
	 * \{
	 */
		//PythonUserFunction(const char* pyCallback, size_t numArgs);
		//PythonUserFunction(const char* pyCallback, size_t numArgs, bool bPosTimeNeed);
		PythonUserFunction(TFunctionHandle handle, size_t numArgs);
		// PythonUserFunction(LuaFunctionHandle handle, size_t numArgs, bool bPosTimeNeed);
	/// \}

	///	destructor frees the reference
		virtual ~PythonUserFunction();

	/**
	 * \brief sets the Lua function used to compute the derivative
	 * \param arg this is the derivative with respect to the parameter index \c arg
	 * \param pyCallback name of the Lua function to use as callback
	 */
	protected:
		int debug() const
		{
			UG_ASSERT(true, "This is a debug function");
			return 0;
		}
		/**
		 * \brief set input value for paramter \c i
		 * \param i parameter index this input is bind to
		 * \{
		 */
		void set_input(size_t i, SmartPtr<CplUserData<TDataIn, dim> > data);
		void set_deriv(size_t arg, TFunctionHandle handle);

		///	set number of needed inputs
		void set_num_input(size_t num);


	public:

		void set_input_and_deriv(size_t i, SmartPtr<CplUserData<TDataIn, dim> > input, TFunctionHandle deriv)
		{ set_input(i, input); set_deriv(i, deriv); }

	///	evaluates the data
		void operator() (TData& out, const std::vector<TDataIn> &) const;
		virtual void operator() (TData& out, int numArgs, ...) const;

		inline void evaluate (TData& value,
		                      const MathVector<dim>& globIP,
		                      number time, int si) const;

		template <int refDim>
		inline void evaluate(TData vValue[],
		                     const MathVector<dim> vGlobIP[],
		                     number time, int si,
		                     GridObject* elem,
		                     const MathVector<dim> vCornerCoords[],
		                     const MathVector<refDim> vLocIP[],
		                     const size_t nip,
		                     LocalVector* u,
		                     const MathMatrix<refDim, dim>* vJT = NULL) const;

		template <int refDim>
		void eval_and_deriv(TData vValue[],
		                    const MathVector<dim> vGlobIP[],
		                    number time, int si,
		                    GridObject* elem,
		                    const MathVector<dim> vCornerCoords[],
		                    const MathVector<refDim> vLocIP[],
		                    const size_t nip,
		                    LocalVector* u,
		                    bool bDeriv,
		                    int s,
		                    std::vector<std::vector<TData> > vvvDeriv[],
		                    const MathMatrix<refDim, dim>* vJT = NULL);

	protected:

	///	frees the callback-reference, if a callback was set.
	//	void free_callback_ref();

	///	frees callback-references for derivate callbacks
	//	void free_deriv_callback_ref(size_t arg);

	///	evaluates the data at a given point and time
		void eval_value(TData& out, const std::vector<TDataIn>& dataIn,
						const MathVector<dim>& x, number time, int si) const;

	///	evaluates the data at a given point and time
		void eval_deriv(TData& out, const std::vector<TDataIn>& dataIn,
						const MathVector<dim>& x, number time, int si, size_t arg) const;

	protected:
		///	number of arguments to use
		size_t m_numArgs;

		///	reference to function callback.
		TFunctionHandle m_cbValueRef;

		/// Reference to dervivative callbacks.
		std::vector<TFunctionHandle> m_cbDerivRef;

		///	data input
		std::vector<SmartPtr<CplUserData<TDataIn, dim> > > m_vpUserData;

		///	data input casted to dependent data
		std::vector<SmartPtr<DependentUserData<TDataIn, dim> > > m_vpDepentData;

};





namespace pybind{

	void RegisterPythonUserData(Registry& reg, std::string grp);

} // end namepace pybind


} // end of namespace ug

#include "python_user_data_impl.h"

#endif // UG_USE_PYBIND11
