// created by Sebastian Reiter
// s.b.reiter@gmail.com
// Mar 6, 2013 (d,m,y)

#include <iostream>
#include <sstream>
#include <vector>
#include <string>

// include bridge
#include "bridge/bridge.h"
#include "bridge/util.h"
#include "bridge/util_domain_dependent.h"

#ifdef UG_PARALLEL
	#include "lib_grid/parallelization/load_balancer.h"
	#include "lib_grid/parallelization/partitioner_bisection.h"
	#ifdef UG_PARMETIS
		#include "lib_grid/parallelization/partitioner_parmetis.h"
	#endif
#endif

using namespace std;

namespace ug{

/**
 * \defgroup loadbalance_bridge Load Balancing Bridge
 * \ingroup domain_bridge
 * \{
 */

static bool MetisIsAvailable()
{
	#ifdef UG_METIS
		return true;
	#endif
	return false;
}

static bool ParmetisIsAvailable()
{
	#ifdef UG_PARMETIS
		return true;
	#endif
	return false;
}

// end group loadbalance_bridge
/// \}

namespace bridge{
namespace LoadBalancing{

/// \addtogroup loadbalance_bridge
/// \{

/**
 * Class exporting the functionality. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */
struct Functionality
{

static void Common(Registry& reg, string grp) {
	reg.add_function("MetisIsAvailable", &MetisIsAvailable, grp);
	reg.add_function("ParmetisIsAvailable", &ParmetisIsAvailable, grp);
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
template <typename TDomain>
static void Domain(Registry& reg, string grp)
{
	string suffix = GetDomainSuffix<TDomain>();
	string tag = GetDomainTag<TDomain>();

	#ifdef UG_PARALLEL
		{
			typedef IPartitioner<TDomain::dim> T;
			string name = string("IPartitioner").append(suffix);
			reg.add_class_<T>(name, grp)
				.add_method("set_verbose", &T::set_verbose);
			reg.add_class_to_group(name, "IPartitioner", tag);
		}

		{
			typedef IPartitioner<TDomain::dim> TBase;
			typedef Partitioner_Bisection<TDomain::dim> T;
			string name = string("Partitioner_Bisection").append(suffix);
			reg.add_class_<T, TBase>(name, grp)
				.add_constructor()
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "Partitioner_Bisection", tag);
		}

		#ifdef UG_PARMETIS
		{
			typedef IPartitioner<TDomain::dim> TBase;
			typedef Partitioner_Parmetis<TDomain::dim> T;
			string name = string("Partitioner_Parmetis").append(suffix);
			reg.add_class_<T, TBase>(name, grp)
				.add_constructor()
				.add_method("set_regard_all_children", &T::set_regard_all_children)
				.add_method("set_child_weight", &T::set_child_weight)
				.add_method("set_sibling_weight", &T::set_sibling_weight)
				.add_method("set_itr_factor", &T::set_itr_factor)
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "Partitioner_Parmetis", tag);
		}
		#endif

		{
		//	Note that this class does not feature a constructor.
		//	One normally uses the derived class DomainLoadBalancer, registered
		//	in bridge/disc_bridges/domain_disc_bridge.cpp
			typedef LoadBalancer<TDomain::dim> T;
			string name = string("LoadBalancer").append(suffix);
			reg.add_class_<T>(name, grp)
					.add_method("add_distribution_level", &T::add_distribution_level)
					.add_method("rebalance", &T::rebalance)
					.add_method("set_balance_threshold", &T::set_balance_threshold)
					.add_method("set_element_threshold", &T::set_element_threshold)
					.add_method("set_partitioner", &T::set_partitioner)
					.add_method("create_quality_record", &T::create_quality_record)
					.add_method("print_quality_records", &T::print_quality_records);

			reg.add_class_to_group(name, "LoadBalancer", tag);
		}
	#endif
}

};// end of struct Functionality

// end group loadbalance_bridge
/// \}

}// end of namespace

/// \addtogroup loadbalance_bridge
void RegisterBridge_LoadBalancing(Registry& reg, string grp)
{
	grp.append("/LoadBalancing");

	typedef LoadBalancing::Functionality Functionality;

	try{
		RegisterCommon<Functionality>(reg, grp);
		RegisterDomainDependent<Functionality>(reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

}// end of namespace
}// end of namespace
