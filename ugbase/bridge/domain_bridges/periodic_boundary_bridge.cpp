/*
 * periodic_boundary_bridge.cpp
 *
 *  Created on: 26.11.2012
 *      Author: marscher
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
	pi.print_identification<EdgeBase>();
	pi.print_identification<VertexBase>();
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
	 * @param parentGroup		group for sorting of functionality
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
