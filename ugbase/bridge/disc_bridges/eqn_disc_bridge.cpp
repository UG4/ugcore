/*
 * Copyright (c) 2022:  G-CSC, Goethe University Frankfurt
 * Author: Felix Salfelder
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

#include <string>

#include "bridge/bridge.h"
#include "bridge/util.h"
#include "bridge/util_domain_algebra_dependent.h"

#include "lib_disc/eqn_disc.h"

namespace {

struct Functionality
{

template <class TDomain, class TAlgebra>
static void DomainAlgebra(ug::bridge::Registry& reg, std::string grp)
{
	std::string suffix = ug::bridge::GetDomainAlgebraSuffix<TDomain,TAlgebra>();
	std::string tag = ug::bridge::GetAlgebraTag<TAlgebra>();

	{
		typedef ug::IEquationDisc<TDomain, TAlgebra> T;
		std::string name = std::string("IEquationDisc");
		reg.add_class_<T>(name+suffix, grp)
			.add_method("get_dt", &T::get_dt)
			.add_method("set_dt", &T::set_dt)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name+suffix, name, tag);
	}
}

}; // Functionality

}

namespace ug{
namespace bridge{

void RegisterBridge_EquationDisc(Registry& reg, std::string grp)
{
	grp.append("/Discretization");

	try{
//		RegisterCommon<Functionality>(reg, grp);
//		RegisterDimensionDependent<Functionality>(reg, grp);
//		RegisterDomainDependent<Functionality>(reg, grp);
//		RegisterAlgebraDependent<Functionality>(reg, grp);
		RegisterDomainAlgebraDependent<Functionality>(reg, grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

} // bridge
} // ug
