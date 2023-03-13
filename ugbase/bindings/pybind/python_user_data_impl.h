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






template <typename TData, int dim, typename TRet = void>
class PythonUserData
	: public StdGlobPosData<PythonUserData<TData, dim, TRet>, TData, dim, TRet>
{

	public:

	///	Constructor
	/**
	 * Creates a PythonUserData that uses a Python function to evaluate some data.
	 */
	typedef py::object TFunction;

	PythonUserData(TFunction f) : func(f)
	{}

	///	evaluates the data at a given point and time
	inline TRet evaluate(TData& value, const MathVector<dim>& x, number time, int si) const
	{
		py::object result_py = PyUserDataTraits<dim>().call(func, x, time, si);
		value = result_py.cast<TData>();
		// std::cerr << "RESULT: " << value << std::endl;
		return;
	};

	//bool requires_grid_fct() const override
	//{ return false; }

	//protected:
	///	sets that PythonUserData is created by DataFactory (not implemented for python)
	//	void set_created_from_factory(bool bFromFactory=false) {m_bFromFactory = bFromFactory;}

	protected:
		TFunction func;


	///	flag, indicating if created from factory
	//	bool m_bFromFactory;

};


namespace pybind{

	void RegisterPythonUserData(Registry& reg, std::string grp);

} // end namepace pybind


} // end of namespace ug

#endif // UG_USE_PYBIND11
