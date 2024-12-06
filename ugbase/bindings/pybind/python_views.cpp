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
#include <sstream>
#include <string>

// Include bridge.
#include "bridge/bridge.h"
#include "bridge/util.h"
#include "bridge/suffix_tag.h"
#include "bridge/util_domain_algebra_dependent.h"


#include "lib_disc/time_disc/theta_time_step.h"
#include "lib_disc/time_disc/time_integrator_observers/python_callback_observer.hpp"

// Include own header.
#include "python_views.h"


using namespace std;

namespace ug {

namespace pybind{

namespace PythonViews {
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

}

template <typename TDomain, typename TAlgebra, typename TRegistry>
static void DomainAlgebra(TRegistry& reg, string parentGroup)
{
	//	typedefs for Vector and Matrix
	typedef typename TAlgebra::vector_type vector_type;

	string suffix = ug::bridge::GetDomainAlgebraSuffix<TDomain,TAlgebra>();
	string tag = ug::bridge::GetDomainAlgebraTag<TDomain,TAlgebra>();

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

	{
		typedef NumpyVectorView<CPUAlgebra::vector_type> T;
		reg.template add_class_<T>("NumpyVectorView", grp)
		  .template add_constructor<void (*)(typename T::TSmartPointer) >("internal id=1")
		  .add_method("get_view", &T::get_view)
		  .add_method("set", &T::set);

	}
} // Common

}; // end Functionality


}// end PythonViews

void RegisterPythonViews(Registry& reg, string grp)
{
	typedef PythonViews::Functionality Functionality;

	try{
		RegisterCommon<Functionality>(reg,grp);
		// RegisterDimensionDependent<Functionality>(reg,grp);
		//ug::bridge::RegisterDomainAlgebraDependent<Functionality, Registry>(reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

}// namespace bridge


}// namespace ug
