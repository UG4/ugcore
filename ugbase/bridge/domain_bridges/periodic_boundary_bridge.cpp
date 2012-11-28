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
#include "lib_grid/tools/periodic_boundary_identifier.h"

using namespace std;

namespace ug {
namespace bridge {
namespace periodicBoundary {
/**
 * Class exporting the functionality. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */
struct Functionality {

	static void Common(Registry& reg, string grp) {
		reg.add_class_<PeriodicBoundaryIdentifier>("PeriodicBoundaryIdentifier", grp)
				.add_constructor();
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
		reg.add_function("IdentifySubsets", &IdentifySubsets<TDomain>, grp);
	}
}; // end Functionality

}  // end periodicBoundary

void RegisterBridge_PeriodicBoundary(Registry& reg, string grp) {
	grp.append("/Periodic");

	typedef boost::mpl::list<
#ifdef UG_DIM_2
	Domain2d
#endif
#if defined UG_DIM_2 && defined UG_DIM_3
	,
#endif
#ifdef UG_DIM_3
	Domain3d
#endif
	> CompileDomain2d3dList; // no periodic boundaries for 1d

	try {
		RegisterCommon<periodicBoundary::Functionality>(reg, grp);

		RegisterDomainDependent<periodicBoundary::Functionality,
				CompileDomain2d3dList>(reg, grp);

	} UG_REGISTRY_CATCH_THROW(grp);
}

} // end of namespace bridge
} // end of namespace ug
