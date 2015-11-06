/*
 * Copyright (c) 2014-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Stepniewski
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

// extern headers
#include <iostream>
#include <sstream>
#include <string>

// include bridge
#include "bridge/bridge.h"
#include "bridge/util.h"
#include "bridge/util_algebra_dependent.h"

// lib_disc includes
#include "lib_disc/domain.h"
#include "lib_disc/spatial_disc/manifold_assemble_util.h"


using namespace std;


namespace ug {
namespace bridge {

/**
 * Class exporting the functionality. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */
////////////////////////////////////////////////////////////////////////////////
//	Functionality
struct Functionality
{
	/**
	 * Function called for the registration of Domain dependent parts
	 * of the plugin. All Functions and Classes depending on the Domain
	 * are to be placed here when registering. The method is called for all
	 * available Domain types, based on the current build options.
	 *
	 * @param reg				registry
	 * @param parentGroup		group for sorting of functionality
	 */
	template <typename TAlgebra>
	static void Algebra(Registry& reg, string grp)
	{
		string suffix = GetAlgebraSuffix<TAlgebra>();
		string tag = GetAlgebraTag<TAlgebra>();

	//	registry of H slave exclusion functionality
		reg.add_function("MarkAllElemsForAssemblyButHSlaves",
				static_cast<void(*)(SmartPtr<IAssemble<TAlgebra> >, Grid&)>
		(&ug::MarkAllElemsForAssemblyButHSlaves<TAlgebra>), grp.c_str());
	}
};






////////////////////////////////////////////////////////////////////////////////
//	Register
void RegisterBridge_ManifoldUtil(Registry& reg, string grp)
{
	typedef Functionality Functionality;

	#if defined(UG_DIM_3)
		try{
			RegisterAlgebraDependent<Functionality>(reg,grp);
		}
		UG_REGISTRY_CATCH_THROW(grp);
	#endif
}

}//	end of namespace bridge
}//	end of namespace ug
