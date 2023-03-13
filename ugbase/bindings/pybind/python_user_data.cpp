/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
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
#include <sstream>
#include <string>

// include bridge
#include "bridge/bridge.h"
#include "bridge/util.h"

#include "python_user_data_impl.h"

using namespace std;

namespace ug
{
	


namespace pybind{

namespace PythonUserData{


template <typename TData, int dim>
void RegisterPythonUserDataType(Registry& reg, string type, string grp)
{
	string suffix = ug::bridge::GetDimensionSuffix<dim>();
	string tag = ug::bridge::GetDimensionTag<dim>();

//	PythonUser"Type"
	{
		typedef ug::PythonUserData<TData, dim> T;
		typedef CplUserData<TData, dim> TBase;
		typedef UserData<TData, dim> TBase2;

		string name = string("PythonUser").append(type).append(suffix);
		reg.add_class_<T, TBase, TBase2>(name, grp)
			.template add_constructor<void (*)(py::object f)>("Callback")
			.add_method("evaluate", &T::evaluate)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, string("PythonUser").append(type), tag);
	}

}

/**
 * Class exporting the functionality. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */
struct Functionality
{

/**
 * Function called for the registration of Dimension dependent parts.
 * All Functions and Classes depending on the Dimension
 * are to be placed here when registering. The method is called for all
 * available Dimension types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <int dim>
static void Dimension(Registry& reg, string grp)
{
	string suffix = ug::bridge::GetDimensionSuffix<dim>();
	string tag = ug::bridge::GetDimensionTag<dim>();

	RegisterPythonUserDataType<number, dim>(reg, "Number", grp);
	//RegisterPythonUserDataType<MathVector<dim>, dim>(reg, "Vector", grp);
	//RegisterPythonUserDataType<MathMatrix<dim,dim>, dim>(reg, "Matrix", grp);

}

/**
 * Function called for the registration of Domain and Algebra independent parts.
 * All Functions and Classes not depending on Domain and Algebra
 * are to be placed here when registering.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
static void Common(Registry& reg, string grp)
{

} // Common

}; // end Functionality
}// end PythonUserData

void RegisterPythonUserData(Registry& reg, string grp)
{
	typedef PythonUserData::Functionality Functionality;

	try{
		// RegisterCommon<Functionality>(reg,grp);
		RegisterDimensionDependent<Functionality>(reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

}// namespace bridge


}// namespace ug
