/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
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
#include <vector>
#include <string>

// include bridge
#include "bridge/bridge.h"
#include "bridge/util.h"
#include "bridge/util_domain_dependent.h"

#include "bindings/lua/lua_user_data.h"

#ifdef UG_PARALLEL
	#include "lib_disc/parallelization/domain_load_balancer.h"
	#include "lib_grid/parallelization/load_balancer.h"
	#include "lib_grid/parallelization/load_balancer_util.h"
	#include "lib_grid/parallelization/partitioner_dynamic_bisection.h"
	#include "lib_grid/parallelization/balance_weights_ref_marks.h"
	#include "lib_grid/parallelization/partition_pre_processors/replace_coordinate.h"
	#include "lib_grid/parallelization/partition_post_processors/smooth_partition_bounds.h"
	#include "lib_grid/parallelization/partition_post_processors/cluster_element_stacks.h"
#endif

using namespace std;

namespace ug{

/**
 * \defgroup loadbalance_bridge Load Balancing Bridge
 * \ingroup domain_bridge
 * \{
 */

#ifdef UG_PARALLEL

	template <int dim>
	class PartPreProc_RasterProjectorCoordinates : public IPartitionPreProcessor
	{
	private:
		typedef PPP_ReplaceCoordinate <dim> 	PPP_RC;
		SmartPtr<PPP_RC>						m_ppp;

	public:
		typedef MathVector<dim>	vector_t;

		PartPreProc_RasterProjectorCoordinates (const Attachment<vector_t> aPos)
		{
			typedef RasterLayersProjector::rel_z_attachment_t attachment_t;

			static const std::string attName("RasterLayersProjector_ARelZ");
			if(!GlobalAttachments::is_declared(attName)){
				GlobalAttachments::declare_attachment<attachment_t>(attName);
			}

			attachment_t aRelZ = GlobalAttachments::attachment<attachment_t>(attName);

			typedef PPP_ReplaceCoordinate <dim> 	PPP_RC;

			m_ppp = make_sp(new PPP_RC (aPos, aRelZ, dim - 1));
		}

		virtual void partitioning_starts (	MultiGrid* mg,
		                                 	IPartitioner* partitioner)
		{
			m_ppp->partitioning_starts (mg, partitioner);
		}


		virtual void partitioning_done (	MultiGrid* mg,
		                                	IPartitioner* partitioner)
		{
			m_ppp->partitioning_done (mg, partitioner);
		}

	};


	template <class TDomain>
	class BalanceWeightsLuaCallback : public IBalanceWeights
	{
		public:
			BalanceWeightsLuaCallback(SmartPtr<TDomain> spDom, const char* luaCallbackName) :
				m_spDom(spDom),
				m_time(0)
			{
				m_pmg = spDom->grid().get();
				m_aaPos = spDom->position_accessor();
			//	we'll pass the following arguments: x, y, z, lvl, t
				m_callback.set_lua_callback(luaCallbackName, 5);
			}

			virtual ~BalanceWeightsLuaCallback()	{}


			void set_time(number time)	{m_time = time;}
			number time() const			{return m_time;}

			virtual number get_weight(Vertex* e)	{return	get_weight_impl(e);}
			virtual number get_weight(Edge* e)		{return	get_weight_impl(e);}
			virtual number get_weight(Face* e)		{return	get_weight_impl(e);}
			virtual number get_weight(Volume* e)	{return	get_weight_impl(e);}

		private:
			typedef typename TDomain::grid_type grid_t;
			typedef typename TDomain::position_type pos_t;
			typedef typename TDomain::position_accessor_type aapos_t;

			template <class TElem>
			number get_weight_impl(TElem* e)
			{
				pos_t c = CalculateCenter(e, m_aaPos);
				vector3 p;
				VecCopy(p, c, 0);
				number weight;
				m_callback(weight, 5, p.x(), p.y(), p.z(), (number)m_pmg->get_level(e), m_time);
				return weight;
			}

			SmartPtr<TDomain>			m_spDom;
			MultiGrid*					m_pmg;
			aapos_t						m_aaPos;
			number						m_time;
			LuaFunction<number, number>	m_callback;
	};
#endif

// end group loadbalance_bridge
/// \}

namespace bridge{
namespace LoadBalancing{

/// \addtogroup loadbalance_bridge
/// \{

#ifdef UG_PARALLEL

template <class TDomain, class TPartitioner>
static void RegisterDynamicBisectionPartitioner(
	Registry& reg,
	string name,
	string grpName,
	string clsGrpName)
{
	reg.add_class_<TPartitioner, IPartitioner>(name, grpName)
		.template add_constructor<void (*)(TDomain&)>()
		.add_method("set_subset_handler",
			&TPartitioner::set_subset_handler)
		.add_method("enable_longest_split_axis",
			&TPartitioner::enable_longest_split_axis)
		.add_method("enable_split_axis",
			&TPartitioner::enable_split_axis)
		.add_method("set_start_split_axis",
			&TPartitioner::set_start_split_axis)
		.add_method("num_split_improvement_iterations",
			&TPartitioner::num_split_improvement_iterations)
		.add_method("set_num_split_improvement_iterations",
			&TPartitioner::set_num_split_improvement_iterations)
		.add_method("enable_static_partitioning",
			&TPartitioner::enable_static_partitioning)
		.add_method("static_partitioning_enabled",
			&TPartitioner::static_partitioning_enabled)
		.set_construct_as_smart_pointer(true);

	reg.add_class_to_group(name, clsGrpName, GetDomainTag<TDomain>());
}

template <class TDomain, class elem_t>
static void RegisterSmoothPartitionBounds(
	Registry& reg,
	string name,
	string grpName,
	string clsGrpName)
{
	typedef SmoothPartitionBounds<elem_t>	T;
	reg.add_class_<T, IPartitionPostProcessor>(name, grpName)
		.add_constructor()
		.set_construct_as_smart_pointer(true);
	reg.add_class_to_group(name, clsGrpName, GetDomainTag<TDomain>());
}

template <class TDomain, class elem_t, class vector_t>
static void RegisterClusterElementStacks(
	Registry& reg,
	string name,
	string grpName,
	string clsGrpName)
{
	typedef Attachment<vector_t> apos_t;
	typedef ClusterElementStacks<elem_t, vector_t>	T;
	reg.add_class_<T, IPartitionPostProcessor>(name, grpName)
		.add_constructor()
		.template add_constructor<void (*)(const apos_t&, const vector_t&)>()
		.add_method("set_position_attachment", &T::set_position_attachment)
		.add_method("set_stacking_direction", &T::set_stacking_direction)
		.set_construct_as_smart_pointer(true);
	reg.add_class_to_group(name, clsGrpName, GetDomainTag<TDomain>());
}

#endif


/**
 * Class exporting the functionality. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */
struct Functionality
{

static void Common(Registry& reg, string grp) {
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
			.add_method("add_partition_hint", &T::add_partition_hint)
			.add_method("partition_hint", &T::partition_hint)
			.set_construct_as_smart_pointer(true);
	}

	{
		reg.add_class_<IBalanceWeights>("IBalanceWeights", grp);
	}

	{
		string name = string("ICommunicationWeights");
		reg.add_class_<ICommunicationWeights>(name, grp);
	}

	{
		string name("BalanceWeightsRefMarks");
		typedef BalanceWeightsRefMarks	T;
		reg.add_class_<T, IBalanceWeights>(name, grp)
			.add_constructor<void (*)(IRefiner*)>()
			.set_construct_as_smart_pointer(true);
	}

	{
		typedef IPartitionPreProcessor T;
		reg.add_class_<T>("IPartitionPreProcessor", grp);
	}

	{
		typedef IPartitionPostProcessor T;
		reg.add_class_<T>("IPartitionPostProcessor", grp);
	}

	{
		typedef IPartitioner T;
		reg.add_class_<T>("IPartitioner", grp)
			.add_method("set_verbose", &T::set_verbose)
			.add_method("partition", &T::partition)
			.add_method("set_balance_weights", &T::set_balance_weights)
			.add_method("supports_balance_weights", &T::supports_balance_weights)
			.add_method("set_communication_weights", &T::set_communication_weights)
			.add_method("supports_communication_weights", &T::supports_communication_weights)
			.add_method("set_next_process_hierarchy", &T::set_next_process_hierarchy)
			.add_method("enable_clustered_siblings", &T::enable_clustered_siblings)
			.add_method("clustered_siblings_enabled", &T::clustered_siblings_enabled)
			.add_method("set_partition_pre_processor", &T::set_partition_pre_processor)
			.add_method("set_partition_post_processor", &T::set_partition_post_processor);
	}

	{
	//	Note that this class does not feature a constructor.
	//	One normally uses the derived class DomainLoadBalancer
		typedef LoadBalancer T;
		reg.add_class_<T>("LoadBalancer", grp)
				//.add_method("add_distribution_level", &T::add_distribution_level)
				.add_method("enable_vertical_interface_creation", &T::enable_vertical_interface_creation)
				.add_method("set_next_process_hierarchy", &T::set_next_process_hierarchy)
				.add_method("rebalance", &T::rebalance)
				.add_method("set_balance_threshold", &T::set_balance_threshold)
				.add_method("set_element_threshold", &T::set_element_threshold)
				.add_method("set_partitioner", &T::set_partitioner)
				.add_method("create_quality_record", &T::create_quality_record)
				.add_method("print_quality_records", &T::print_quality_records)
				.add_method("estimate_distribution_quality", static_cast<number (T::*)()>(&T::estimate_distribution_quality))
				.add_method("set_balance_weights", &T::set_balance_weights)
				.add_method("problems_occurred", &T::problems_occurred);
	}

	#ifdef UG_DIM_1
	{
		typedef ug::Domain<1>	TDomain;
		RegisterDynamicBisectionPartitioner<
				TDomain,
				DomainPartitioner<TDomain, Partitioner_DynamicBisection<Edge, 1> > >(
			reg,
			"EdgePartitioner_DynamicBisection1d",
			grp,
			"Partitioner_DynamicBisection");


		RegisterSmoothPartitionBounds<TDomain, Edge>(
			reg,
			"SmoothPartitionBounds1d",
			grp,
			"SmoothPartitionBounds");

		// RegisterClusterElementStacks<TDomain, Edge, vector1>(
		// 	reg,
		// 	"ClusterElementStacks1d",
		// 	grp,
		// 	"ClusterElementStacks");
	}
	#endif
	#ifdef UG_DIM_2
	{
		typedef ug::Domain<2>	TDomain;

		RegisterDynamicBisectionPartitioner<
				TDomain,
				DomainPartitioner<TDomain, Partitioner_DynamicBisection<Edge, 2> > >(
			reg,
			"EdgePartitioner_DynamicBisection2d",
			grp,
			"ManifoldPartitioner_DynamicBisection");

		RegisterDynamicBisectionPartitioner<
				TDomain,
				DomainPartitioner<TDomain, Partitioner_DynamicBisection<Face, 2> > >(
			reg,
			"FacePartitioner_DynamicBisection2d",
			grp,
			"Partitioner_DynamicBisection");

		RegisterSmoothPartitionBounds<TDomain, Face>(
			reg,
			"SmoothPartitionBounds2d",
			grp,
			"SmoothPartitionBounds");

		RegisterClusterElementStacks<TDomain, Face, vector2>(
			reg,
			"ClusterElementStacks2d",
			grp,
			"ClusterElementStacks");
	}
	#endif
	#ifdef UG_DIM_3
	{
		typedef ug::Domain<3>	TDomain;

		RegisterDynamicBisectionPartitioner<
				TDomain,
				DomainPartitioner<TDomain, Partitioner_DynamicBisection<Edge, 3> > >(
			reg,
			"EdgePartitioner_DynamicBisection3d",
			grp,
			"HyperManifoldPartitioner_DynamicBisection");

		RegisterDynamicBisectionPartitioner<
				TDomain,
				DomainPartitioner<TDomain, Partitioner_DynamicBisection<Face, 3> > >(
			reg,
			"FacePartitioner_DynamicBisection3d",
			grp,
			"ManifoldPartitioner_DynamicBisection");

		RegisterDynamicBisectionPartitioner<
				TDomain,
				DomainPartitioner<TDomain, Partitioner_DynamicBisection<Volume, 3> > >(
			reg,
			"VolumePartitioner_DynamicBisection3d",
			grp,
			"Partitioner_DynamicBisection");

		RegisterSmoothPartitionBounds<TDomain, Volume>(
			reg,
			"SmoothPartitionBounds3d",
			grp,
			"SmoothPartitionBounds");

		RegisterClusterElementStacks<TDomain, Volume, vector3>(
			reg,
			"ClusterElementStacks3d",
			grp,
			"ClusterElementStacks");
	}
	#endif



	// #ifdef UG_DIM_1
	// {
	// 	typedef ug::Domain<1>	TDomain;
	// 	string tag = GetDomainTag<TDomain>();
	// 	typedef DomainPartitioner<TDomain, Partitioner_DynamicBisection<Edge, 1> > T;
	// 	string name = string("EdgePartitioner_DynamicBisection1d");
	// 	reg.add_class_<T, IPartitioner>(name, grp)
	// 		.add_constructor<void (*)(TDomain&)>()
	// 		.add_method("enable_static_partitioning", &T::enable_static_partitioning)
	// 		.add_method("static_partitioning_enabled", &T::static_partitioning_enabled)
	// 		.add_method("set_subset_handler", &T::set_subset_handler)
	// 		.add_method("num_split_improvement_iterations", &T::num_split_improvement_iterations)
	// 		.add_method("set_num_split_improvement_iterations", &T::set_num_split_improvement_iterations)
	// 		.set_construct_as_smart_pointer(true);
	// 	reg.add_class_to_group(name, "Partitioner_DynamicBisection", tag);
	// }
	// #endif
	// #ifdef UG_DIM_2
	// {
	// 	typedef ug::Domain<2>	TDomain;
	// 	string tag = GetDomainTag<TDomain>();
	// 	{
	// 		typedef DomainPartitioner<TDomain, Partitioner_DynamicBisection<Edge, 2> >T;
	// 		string name = string("EdgePartitioner_DynamicBisection2d");
	// 		reg.add_class_<T, IPartitioner>(name, grp)
	// 			.add_constructor<void (*)(TDomain&)>()
	// 			.add_method("enable_static_partitioning", &T::enable_static_partitioning)
	// 			.add_method("static_partitioning_enabled", &T::static_partitioning_enabled)
	// 			.add_method("set_subset_handler", &T::set_subset_handler)
	// 			.add_method("num_split_improvement_iterations", &T::num_split_improvement_iterations)
	// 			.add_method("set_num_split_improvement_iterations", &T::set_num_split_improvement_iterations)
	// 			.set_construct_as_smart_pointer(true);
	// 		reg.add_class_to_group(name, "ManifoldPartitioner_DynamicBisection", tag);
	// 	}
	// 	{
	// 		typedef DomainPartitioner<TDomain, Partitioner_DynamicBisection<Face, 2> > T;
	// 		string name = string("FacePartitioner_DynamicBisection2d");
	// 		reg.add_class_<T, IPartitioner>(name, grp)
	// 			.add_constructor<void (*)(TDomain&)>()
	// 			.add_method("enable_static_partitioning", &T::enable_static_partitioning)
	// 			.add_method("static_partitioning_enabled", &T::static_partitioning_enabled)
	// 			.add_method("set_subset_handler", &T::set_subset_handler)
	// 			.add_method("num_split_improvement_iterations", &T::num_split_improvement_iterations)
	// 			.add_method("set_num_split_improvement_iterations", &T::set_num_split_improvement_iterations)
	// 			.set_construct_as_smart_pointer(true);
	// 		reg.add_class_to_group(name, "Partitioner_DynamicBisection", tag);
	// 	}
	// }
	// #endif
	// #ifdef UG_DIM_3
	// {
	// 	typedef ug::Domain<3>	TDomain;
	// 	string tag = GetDomainTag<TDomain>();
	// 	{
	// 		typedef DomainPartitioner<TDomain, Partitioner_DynamicBisection<Edge, 3> > T;
	// 		string name = string("EdgePartitioner_DynamicBisection3d");
	// 		reg.add_class_<T, IPartitioner>(name, grp)
	// 			.add_constructor<void (*)(TDomain&)>()
	// 			.add_method("enable_static_partitioning", &T::enable_static_partitioning)
	// 			.add_method("static_partitioning_enabled", &T::static_partitioning_enabled)
	// 			.add_method("set_subset_handler", &T::set_subset_handler)
	// 			.add_method("num_split_improvement_iterations", &T::num_split_improvement_iterations)
	// 			.add_method("set_num_split_improvement_iterations", &T::set_num_split_improvement_iterations)
	// 			.set_construct_as_smart_pointer(true);
	// 		reg.add_class_to_group(name, "HyperManifoldPartitioner_DynamicBisection", tag);
	// 	}
	// 	{
	// 		typedef DomainPartitioner<TDomain, Partitioner_DynamicBisection<Face, 3> > T;
	// 		string name = string("FacePartitioner_DynamicBisection3d");
	// 		reg.add_class_<T, IPartitioner>(name, grp)
	// 			.add_constructor<void (*)(TDomain&)>()
	// 			.add_method("enable_static_partitioning", &T::enable_static_partitioning)
	// 			.add_method("static_partitioning_enabled", &T::static_partitioning_enabled)
	// 			.add_method("set_subset_handler", &T::set_subset_handler)
	// 			.add_method("num_split_improvement_iterations", &T::num_split_improvement_iterations)
	// 			.add_method("set_num_split_improvement_iterations", &T::set_num_split_improvement_iterations)
	// 			.set_construct_as_smart_pointer(true);
	// 		reg.add_class_to_group(name, "ManifoldPartitioner_DynamicBisection", tag);
	// 	}
	// 	{
	// 		typedef DomainPartitioner<TDomain, Partitioner_DynamicBisection<Volume, 3> > T;
	// 		string name = string("VolumePartitioner_DynamicBisection3d");
	// 		reg.add_class_<T, IPartitioner>(name, grp)
	// 			.add_constructor<void (*)(TDomain&)>()
	// 			.add_method("enable_static_partitioning", &T::enable_static_partitioning)
	// 			.add_method("static_partitioning_enabled", &T::static_partitioning_enabled)
	// 			.add_method("set_subset_handler", &T::set_subset_handler)
	// 			.add_method("num_split_improvement_iterations", &T::num_split_improvement_iterations)
	// 			.add_method("set_num_split_improvement_iterations", &T::set_num_split_improvement_iterations)
	// 			.set_construct_as_smart_pointer(true);
	// 		reg.add_class_to_group(name, "Partitioner_DynamicBisection", tag);
	// 	}
	// }
	// #endif
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

		{
			typedef PartPreProc_RasterProjectorCoordinates<TDomain::dim> T;
			string name = string("PartPreProc_RasterProjectorCoordinates").append(suffix);
			reg.add_class_<T, IPartitionPreProcessor>(name, grp)
				.template add_constructor<void (*)(const Attachment<MathVector<TDomain::dim> >&)>()
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "PartPreProc_RasterProjectorCoordinates", tag);
		}

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
			typedef BalanceWeightsLuaCallback<TDomain> T;
			string name = string("BalanceWeightsLuaCallback").append(suffix);
			reg.add_class_<T, IBalanceWeights>(name, grp)
				.template add_constructor<void (*)(SmartPtr<TDomain> spDom,
												   const char* luaCallbackName)>()
				.add_method("set_time", &T::set_time)
				.add_method("time", &T::time)
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "BalanceWeightsLuaCallback", tag);
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

		reg.add_function("CreateProcessHierarchy",
						 static_cast<SPProcessHierarchy (*)(TDomain&, size_t,
						 									size_t, size_t, int,
						 									int)>
						 	(&CreateProcessHierarchy<TDomain>),
						 grp, "ProcessHierarchy", "Domain # minNumElemsPerProcPerLvl # "
						 "maxNumRedistProcs # maxNumProcs # minDistLvl # "
						 "maxLvlsWithoutRedist");
		reg.add_function("CreateProcessHierarchy",
						 static_cast<SPProcessHierarchy (*)(TDomain&, size_t,
						 									size_t, size_t, int,
						 									int, IRefiner*)>
						 	(&CreateProcessHierarchy<TDomain>),
						 grp, "ProcessHierarchy", "Domain # minNumElemsPerProcPerLvl # "
						 "maxNumRedistProcs # maxNumProcs # minDistLvl # "
						 "maxLvlsWithoutRedist # refiner");

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
