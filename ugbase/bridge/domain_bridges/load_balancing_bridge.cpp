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
	#include "lib_grid/parallelization/load_balancer_util.h"
	#include "lib_grid/parallelization/partitioner_dynamic_bisection.h"
	#include "lib_disc/parallelization/domain_load_balancer.h"
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

	#ifdef UG_PARALLEL
	{
		typedef ProcessHierarchy T;
		reg.add_class_<T>("ProcessHierarchy", grp)
			.add_constructor()
			.add_method("empty", &T::empty)
			.add_method("add_hierarchy_level", &T::add_hierarchy_level)
			.add_method("num_hierarchy_levels", &T::num_hierarchy_levels)
			.add_method("num_global_procs_involved", &T::num_global_procs_involved)
			.add_method("grid_base_level", &T::grid_base_level)
			.add_method("hierarchy_level_from_grid_level", &T::hierarchy_level_from_grid_level)
			.add_method("cluster_procs", &T::cluster_procs)
			.add_method("to_string", &T::to_string)
			.set_construct_as_smart_pointer(true);
	}

	{
		reg.add_class_<IBalanceWeights>("IBalanceWeights", grp);
	}

	{
		typedef IPartitioner T;
		reg.add_class_<T>("IPartitioner", grp)
			.add_method("set_verbose", &T::set_verbose)
			.add_method("partition", &T::partition)
			.add_method("set_next_process_hierarchy", &T::set_next_process_hierarchy);
	}

	{
	//	Note that this class does not feature a constructor.
	//	One normally uses the derived class DomainLoadBalancer
		typedef LoadBalancer T;
		reg.add_class_<T>("LoadBalancer", grp)
				//.add_method("add_distribution_level", &T::add_distribution_level)
				.add_method("set_next_process_hierarchy", &T::set_next_process_hierarchy)
				.add_method("rebalance", &T::rebalance)
				.add_method("set_balance_threshold", &T::set_balance_threshold)
				.add_method("set_element_threshold", &T::set_element_threshold)
				.add_method("set_partitioner", &T::set_partitioner)
				.add_method("create_quality_record", &T::create_quality_record)
				.add_method("print_quality_records", &T::print_quality_records)
				.add_method("set_balance_weights", &T::set_balance_weights);
	}

	#ifdef UG_DIM_1
	{
		typedef ug::Domain<1>	TDomain;
		string tag = GetDomainTag<TDomain>();
		typedef DomainPartitioner<TDomain, Partitioner_DynamicBisection<Edge, 1> > T;
		string name = string("EdgePartitioner_DynamicBisection1d");
		reg.add_class_<T, IPartitioner>(name, grp)
			.add_constructor<void (*)(TDomain&)>()
			.add_method("enable_static_partitioning", &T::enable_static_partitioning)
			.add_method("static_partitioning_enabled", &T::static_partitioning_enabled)
			.add_method("set_subset_handler", &T::set_subset_handler)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "EdgePartitioner_DynamicBisection", tag);
		reg.add_class_to_group(name, "Partitioner_DynamicBisection", tag);
	}
	#endif
	#ifdef UG_DIM_2
	{
		typedef ug::Domain<2>	TDomain;
		string tag = GetDomainTag<TDomain>();
		{
			typedef DomainPartitioner<TDomain, Partitioner_DynamicBisection<Edge, 2> >T;
			string name = string("EdgePartitioner_DynamicBisection2d");
			reg.add_class_<T, IPartitioner>(name, grp)
				.add_constructor<void (*)(TDomain&)>()
				.add_method("enable_static_partitioning", &T::enable_static_partitioning)
				.add_method("static_partitioning_enabled", &T::static_partitioning_enabled)
				.add_method("set_subset_handler", &T::set_subset_handler)
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "EdgePartitioner_DynamicBisection", tag);
		}
		{
			typedef DomainPartitioner<TDomain, Partitioner_DynamicBisection<Face, 2> > T;
			string name = string("FacePartitioner_DynamicBisection2d");
			reg.add_class_<T, IPartitioner>(name, grp)
				.add_constructor<void (*)(TDomain&)>()
				.add_method("enable_static_partitioning", &T::enable_static_partitioning)
				.add_method("static_partitioning_enabled", &T::static_partitioning_enabled)
				.add_method("set_subset_handler", &T::set_subset_handler)
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "FacePartitioner_DynamicBisection", tag);
			reg.add_class_to_group(name, "Partitioner_DynamicBisection", tag);
		}
	}
	#endif
	#ifdef UG_DIM_3
	{
		typedef ug::Domain<3>	TDomain;
		string tag = GetDomainTag<TDomain>();
		{
			typedef DomainPartitioner<TDomain, Partitioner_DynamicBisection<Edge, 3> > T;
			string name = string("EdgePartitioner_DynamicBisection3d");
			reg.add_class_<T, IPartitioner>(name, grp)
				.add_constructor<void (*)(TDomain&)>()
				.add_method("enable_static_partitioning", &T::enable_static_partitioning)
				.add_method("static_partitioning_enabled", &T::static_partitioning_enabled)
				.add_method("set_subset_handler", &T::set_subset_handler)
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "EdgePartitioner_DynamicBisection", tag);
		}
		{
			typedef DomainPartitioner<TDomain, Partitioner_DynamicBisection<Face, 3> > T;
			string name = string("FacePartitioner_DynamicBisection3d");
			reg.add_class_<T, IPartitioner>(name, grp)
				.add_constructor<void (*)(TDomain&)>()
				.add_method("enable_static_partitioning", &T::enable_static_partitioning)
				.add_method("static_partitioning_enabled", &T::static_partitioning_enabled)
				.add_method("set_subset_handler", &T::set_subset_handler)
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "FacePartitioner_DynamicBisection", tag);
		}
		{
			typedef DomainPartitioner<TDomain, Partitioner_DynamicBisection<Volume, 3> > T;
			string name = string("VolumePartitioner_DynamicBisection3d");
			reg.add_class_<T, IPartitioner>(name, grp)
				.add_constructor<void (*)(TDomain&)>()
				.add_method("enable_static_partitioning", &T::enable_static_partitioning)
				.add_method("static_partitioning_enabled", &T::static_partitioning_enabled)
				.add_method("set_subset_handler", &T::set_subset_handler)
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "VolumePartitioner_DynamicBisection", tag);
			reg.add_class_to_group(name, "Partitioner_DynamicBisection", tag);
		}
	}
	#endif
	#endif
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

		#ifdef UG_PARMETIS
		{
			typedef DomainPartitioner<TDomain, Partitioner_Parmetis<TDomain::dim> > T;
			string name = string("Partitioner_Parmetis").append(suffix);
			reg.add_class_<T, IPartitioner>(name, grp)
				.template add_constructor<void (*)(TDomain&)>()
				.add_method("set_child_weight", &T::set_child_weight)
				.add_method("set_sibling_weight", &T::set_sibling_weight)
				.add_method("set_itr_factor", &T::set_itr_factor)
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "Partitioner_Parmetis", tag);
		}
		#endif

		{
			typedef DomainBalanceWeights<TDomain, AnisotropicBalanceWeights<TDomain::dim> > T;
			string name = string("AnisotropicBalanceWeights").append(suffix);
			reg.add_class_<T, IBalanceWeights>(name, grp)
				.template add_constructor<void (*)(TDomain&)>()
				.add_method("set_weight_factor", &T::set_weight_factor)
				.add_method("weight_factor", &T::weight_factor)
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "AnisotropicBalanceWeights", tag);
		}

		{
			string name = string("DomainLoadBalancer").append(suffix);
			typedef DomainLoadBalancer<TDomain> T;
			typedef LoadBalancer TBase;
			reg.add_class_<T, TBase>(name, grp)
				.template add_constructor<void (*)(SmartPtr<TDomain>)>("Domain")
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "DomainLoadBalancer", tag);
		}

		reg.add_function("CreateProcessHierarchy", &CreateProcessHierarchy<TDomain>, grp,
						 "ProcessHierarchy", "Domain, minNumElemsPerProcPerLvl, "
						 "maxNumRedistProcs, maxNumProcs, minDistLvl, "
						 "maxLvlsWithoutRedist");

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
