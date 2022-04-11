/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
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

// extern headers
#include <iostream>
#include <sstream>
#include <string>

// include bridge
#include "bridge/bridge.h"
#include "bridge/util.h"
#include "bridge/util_domain_dependent.h"
#include "bridge/util_domain_algebra_dependent.h"

// lib_disc includes
#include "lib_disc/domain.h"
//#include "lib_disc/function_spaces/approximation_space.h"

// ordering algorithms
#include "lib_disc/ordering_strategies/algorithms/ordering_algorithms.cpp"

#include "lib_disc/ordering_strategies/io_grid_points_ordering.cpp"
#include "lib_disc/ordering_strategies/io_grid_function_ordering.cpp"
#include "lib_disc/ordering_strategies/io_sorted_grid_function_ordering.cpp"

using namespace std;

namespace ug{
namespace bridge{
namespace Ordering{

/**
 * \defgroup ordering_bridge Ordering Bridge
 * \ingroup disc_bridge
 * \{
 */

/**
 * Class exporting the functionality. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */
struct Functionality
{

/**
 * Function called for the registration of Domain and Algebra dependent parts.
 * All Functions and Classes depending on both Domain and Algebra
 * are to be placed here when registering. The method is called for all
 * available Domain and Algebra types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TDomain, typename TAlgebra>
static void DomainAlgebra(Registry& reg, string grp)
{
	string suffix = GetDomainAlgebraSuffix<TDomain,TAlgebra>();
	string tag = GetDomainAlgebraTag<TDomain,TAlgebra>();

	typedef ug::GridFunction<TDomain, TAlgebra> TFct;

	typedef std::vector<size_t> ordering_container_type;

	typedef ug::GridFunction<TDomain, TAlgebra> TFct;

	typedef SmartPtr<UserData<MathVector<TDomain::dim>, TDomain::dim> > TSpUserData;

	typedef MathVector<TDomain::dim> small_vec_t;

//	Lexicographic ordering
	{
		typedef LexOrdering<TAlgebra, TDomain, ordering_container_type> T;
		typedef IOrderingAlgorithm<TAlgebra, ordering_container_type> TBase;
		string name = string("LexOrdering").append(suffix);
		reg.add_class_<T, TBase>(name, grp, "LexOrdering")
			.add_constructor()
			.add_method("set_direction", &T::set_direction)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "LexOrdering", tag);
	}

	{
		typedef DirectionalOrdering<TAlgebra, TDomain, ordering_container_type> T;
		typedef IOrderingAlgorithm<TAlgebra, ordering_container_type> TBase;
		string name = string("DirectionalOrdering").append(suffix);
		reg.add_class_<T, TBase>(name, grp, "DirectionalOrdering")
			.add_constructor()
			//.add_method("set_direction", static_cast<void (T::*)(const char*)>(&T::set_direction))
			//.add_method("set_direction", static_cast<void (T::*)(TSpUserData)>(&T::set_direction))
			.add_method("set_direction", static_cast<void (T::*)(small_vec_t*)>(&T::set_direction))
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "DirectionalOrdering", tag);
	}

//	Boost Dirichlet Cuthill-McKee ordering
	{
		typedef BoostDirichletCuthillMcKeeOrdering<TAlgebra, TDomain, ordering_container_type> T;
		typedef IOrderingAlgorithm<TAlgebra, ordering_container_type> TBase;
		string name = string("BoostDirichletCuthillMcKeeOrdering").append(suffix);
		reg.add_class_<T, TBase>(name, grp, "LibDiscCuthillMcKeeOrdering")
			.add_constructor()
			.add_method("set_reverse", &T::set_reverse)
			.add_method("select_dirichlet_subset", static_cast<void (T::*)(const char*)>(&T::select_dirichlet_subset))
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "BoostDirichletCuthillMcKeeOrdering", tag);
	}

//	River ordering (topological ordering beginning at selected sources)
	{
		typedef RiverOrdering<TAlgebra, TDomain, ordering_container_type> T;
		typedef IOrderingAlgorithm<TAlgebra, ordering_container_type> TBase;
		string name = string("RiverOrdering").append(suffix);
		reg.add_class_<T, TBase>(name, grp, "RiverOrdering")
			.add_constructor()
			.add_method("select_sources", static_cast<void (T::*)(const char*)>(&T::select_sources))
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "RiverOrdering", tag);
	}

//	Weighted Cuthill-McKee ordering
	{
		typedef WeightedCuthillMcKeeOrdering<TAlgebra, TDomain, ordering_container_type> T;
		typedef IOrderingAlgorithm<TAlgebra, ordering_container_type> TBase;
		string name = string("WeightedCuthillMcKeeOrdering").append(suffix);
		reg.add_class_<T, TBase>(name, grp, "WeightedCuthillMcKeeOrdering")
			.add_constructor()
			.add_method("set_reverse", &T::set_reverse)
			.add_method("select_dirichlet_subset", static_cast<void (T::*)(int)>(&T::select_dirichlet_subset))
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "WeightedCuthillMcKeeOrdering", tag);
	}

//	FollowConvection ordering
	{
		typedef FollowConvectionOrdering<TAlgebra, TDomain, ordering_container_type> T;
		typedef IOrderingAlgorithm<TAlgebra, ordering_container_type> TBase;
		string name = string("FollowConvectionOrdering").append(suffix);
		reg.add_class_<T, TBase>(name, grp, "FollowConvectionOrdering")
			.add_constructor()
			.add_method("select_dirichlet_subset", static_cast<void (T::*)(const char*)>(&T::select_dirichlet_subset))
			.add_method("set_velocity", static_cast<void (T::*)(const char*)>(&T::set_velocity))
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "FollowConvectionOrdering", tag);
	}


	/* Preprocessor */

//	Angle preprocessor
	{
		typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
			boost::property<boost::vertex_color_t,
			boost::default_color_type,
			boost::property<boost::vertex_degree_t, int> > > G_t;
		//typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS> G_t;
		typedef AnglePreprocessor<TAlgebra, TDomain, G_t> T;
		typedef IOrderingPreprocessor<TAlgebra, G_t> TBase;
		string name = string("AnglePreprocessor").append(suffix);
		reg.add_class_<T, TBase>(name, grp, "AnglePreprocessor")
			.add_constructor()
			.add_method("set_velocity", static_cast<void (T::*)(const char*)>(&T::set_velocity))
			.add_method("set_threshold", static_cast<void (T::*)(number)>(&T::set_threshold))
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "AnglePreprocessor", tag);
	}

	/* IO */

//	GridPointsOrdering
	{
		typedef GridPointsOrdering<TDomain, TAlgebra> T;
		string name = string("GridPointsOrdering").append(suffix);
		reg.add_class_<T>(name, grp)
			.template add_constructor<void (*)(SmartPtr<TFct>, const char*)>("GridPointsOrdering")
			.add_method("get", &T::get)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "GridPointsOrdering", tag);
	}

//	GridFunctionOrdering
	{
		typedef GridFunctionOrdering<TDomain, TAlgebra> T;
		string name = string("GridFunctionOrdering").append(suffix);
		reg.add_class_<T>(name, grp)
			.template add_constructor<void (*)(SmartPtr<TFct>, const char*)>("GridFunctionOrdering")
			.add_method("get", &T::get)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "GridFunctionOrdering", tag);
	}

//	SortedGridFunctionOrdering
	{
		typedef SortedGridFunctionOrdering<TDomain, TAlgebra> T;
		typedef IOrderingAlgorithm<TAlgebra, ordering_container_type> TOrdAlgo;
		string name = string("SortedGridFunctionOrdering").append(suffix);
		reg.add_class_<T>(name, grp)
			.template add_constructor<void (*)(SmartPtr<TFct>, SmartPtr<TOrdAlgo>, const char*)>("SortedGridFunctionOrdering")
			.add_method("get", &T::get)
			.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "SortedGridFunctionOrdering", tag);
	}
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
	typedef ApproximationSpace<TDomain> approximation_space_type;

//	group string
	grp.append("/ApproximationSpace");

//	Order Cuthill-McKee
	{
		reg.add_function("OrderCuthillMcKee", static_cast<void (*)(approximation_space_type&, bool)>(&OrderCuthillMcKee), grp);
	}

//	Order lexicographically
	{
		reg.add_function("OrderLex", static_cast<void (*)(approximation_space_type&, const char*)>(&OrderLex<TDomain>), grp);
	}
//	Order in downwind direction
	{
		reg.add_function("OrderDownwind", static_cast<void (*)(approximation_space_type&, SmartPtr<UserData<MathVector<TDomain::dim>, TDomain::dim> >)> (&ug::OrderDownwind<TDomain>), grp);
		reg.add_function("OrderDownwind", static_cast<void (*)(approximation_space_type&, const std::vector<number>&)>(&ug::OrderDownwind<TDomain>), grp);

        reg.add_function("OrderDownwind", static_cast<void (*)(approximation_space_type&, SmartPtr<UserData<MathVector<TDomain::dim>, TDomain::dim> >, number)> (&ug::OrderDownwind<TDomain>), grp);
        reg.add_function("OrderDownwind", static_cast<void (*)(approximation_space_type&, const std::vector<number>&, number)>(&ug::OrderDownwind<TDomain>), grp);

#ifdef UG_FOR_LUA
		reg.add_function("OrderDownwind", static_cast<void (*)(approximation_space_type&, const char*)>(&ug::OrderDownwind<TDomain>), grp);
        reg.add_function("OrderDownwind", static_cast<void (*)(approximation_space_type&, const char*, number)>(&ug::OrderDownwind<TDomain>), grp);
#endif
	}

}

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
	string suffix = GetDimensionSuffix<dim>();
	string tag = GetDimensionTag<dim>();

}

/**
 * Function called for the registration of Algebra dependent parts.
 * All Functions and Classes depending on Algebra
 * are to be placed here when registering. The method is called for all
 * available Algebra types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TAlgebra>
static void Algebra(Registry& reg, string grp)
{
	string suffix = GetAlgebraSuffix<TAlgebra>();
	string tag = GetAlgebraTag<TAlgebra>();

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
}

}; // end Functionality

// end group ordering_bridge
/// \}

}// namespace Ordering

/// \addtogroup ordering_bridge
void RegisterBridge_Ordering(Registry& reg, string grp)
{
	grp.append("/Discretization");
	typedef Ordering::Functionality Functionality;

	try{
		//RegisterCommon<Functionality>(reg,grp);
//		RegisterDimensionDependent<Functionality>(reg,grp);
		RegisterDomainDependent<Functionality>(reg,grp);
//		RegisterAlgebraDependent<Functionality>(reg,grp);
		RegisterDomainAlgebraDependent<Functionality>(reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

}//	end of namespace bridge
}//	end of namespace ug
