/*
 * Copyright (c) 2025:  Goethe University Frankfurt
 * Author: Arne Naegel
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
#include <pybind11/numpy.h> // array_t


namespace py = pybind11;
namespace ug {
////////////////////////////////////////////////////////////////////////////////
// Executing Python Code from Ug4
// Code has been passed as String to UG4
// More Info here: https://docs.python.org/3/extending/embedding.html#
// and also here: https://docs.python.org/3/c-api/ 
// also make sure to have python3-dev installed to use the python header:
// sudo apt install python3-dev
////////////////////////////////////////////////////////////////////////////////



namespace pybind{

namespace PythonViews {


//! The NumpyVectorView provides an object view for Python.
/*! The ownership of the object remains in UG4. Creation is based on a SmarPtr.
 *  After the object is destroyed, the reference is deleted (and the UG4 object can be freed)
 * */
template <typename TVector>
class NumpyVectorView
{
public:
	typedef SmartPtr<TVector> TSmartPointer;

	//! Create view. (Note: Storing reference prevents deletion of underlying object)
	NumpyVectorView(TSmartPointer vec)
	: sp(vec),
	  view(py::memoryview::from_buffer(
			  &(*sp)[0],
			  sizeof(typename TVector::value_type),
			  py::format_descriptor<typename TVector::value_type>::format().c_str(),
			  { sp->size()},
			  {sizeof(typename TVector::value_type)})) {};


	//! Destroy view.
	~NumpyVectorView()
	{ sp = SPNULL; }

	//! Get a python view.
	py::memoryview& memory()
	{ return view; }

protected:
	//! Auxiliary
	TVector &reference()
	{ return (*sp)[0]; }

	// Data structures.
	TSmartPointer sp; 		// Reference to UG4 object.
	py::memoryview view;	// Python memory view.

};




} // end namepace PythonViews

void RegisterPythonViews(Registry& reg, std::string grp);

} // end namepace pybind


} // end of namespace ug


#endif // UG_USE_PYBIND11
