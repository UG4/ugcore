/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Scherer
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

#include "registry/registry.h"
#include "bridge/bridge.h"
#include "bridge/util.h"
#include "bridge/util_domain_dependent.h"
#include "lib_grid/tools/periodic_boundary_manager.h"

using namespace std;

namespace ug {
namespace bridge {
namespace periodicBoundary {

/**
 * \defgroup periodic_bridge Periodic Bounadry Bridge
 * \ingroup domain_bridge
 * \{
 */

void print_all_identifications(Grid& g) {
	if(!g.has_periodic_boundaries())
		return;
	PeriodicBoundaryManager& pi = *g.periodic_boundary_manager();
	pi.print_identification<Face>();
	pi.print_identification<Edge>();
	pi.print_identification<Vertex>();
}

/**
 * Class exporting the functionality. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */
struct Functionality {

	static void Common(Registry& reg, string grp) {
		reg.add_class_<PeriodicBoundaryManager>("PeriodicBoundaryManager", grp)
				.add_constructor();
		reg.add_function("PrintIdentification", &print_all_identifications, grp);
	}

	/**
	 * Function called for the registration of Domain dependent parts.
	 * All Functions and Classes depending on the Domain
	 * are to be placed here when registering. The method is called for all
	 * available Domain types, based on the current build options.
	 *
	 * @param reg				registry
	 * @param grp				group for sorting of functionality
	 */
	template<typename TDomain>
	static void Domain(Registry& reg, string grp) {
		reg.add_function("IdentifySubsets",
				static_cast<void(*)(TDomain&, int, int)>(&IdentifySubsets<TDomain>), grp)
		   .add_function("IdentifySubsets",
				static_cast<void(*)(TDomain&, const char*, const char*)>(&IdentifySubsets<TDomain>), grp);
	}
}; // end Functionality

// end group periodic_bridge
/// \}

}  // end periodicBoundary

/// \addtogroup periodic_bridge
void RegisterBridge_PeriodicBoundary(Registry& reg, string grp) {
	grp.append("/Periodic");

	try {
		RegisterCommon<periodicBoundary::Functionality>(reg, grp);
		RegisterDomainDependent<periodicBoundary::Functionality>(reg, grp);
	} UG_REGISTRY_CATCH_THROW(grp);
}

} // end of namespace bridge
} // end of namespace ug
