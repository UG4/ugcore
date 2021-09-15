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
#include "bridge/util_algebra_dependent.h"


#include "lib_algebra/lib_algebra.h"

// ordering algorithms
#include "lib_algebra/ordering_strategies/algorithms/ordering_algorithms.cpp"

using namespace std;

namespace ug{
namespace bridge{
namespace Ordering{

/**
 * \defgroup ordering_bridge Ordering Bridge
 * \ingroup ordering_bridge
 * \{
 */

/**
 * Class exporting the functionality. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */
struct Functionality
{

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

//	typedefs for this algebra
	typedef typename TAlgebra::vector_type vector_type;
	typedef typename TAlgebra::matrix_type matrix_type;
	typedef std::vector<size_t> ordering_container_type;

	{
		typedef IOrderingAlgorithm<typename TAlgebra::matrix_type, ordering_container_type> TBase;
		reg.add_class_<TBase>("IOrderingAlgorithm", grp);
	}

//	Weighted Cuthill McKee
	{
		typedef WeightedCuthillMcKeeOrdering<typename TAlgebra::matrix_type, ordering_container_type> T;
		typedef IOrderingAlgorithm<typename TAlgebra::matrix_type, ordering_container_type> TBase;
		string name = string("WeightedCuthillMcKeeOrdering").append(suffix);
		reg.add_class_<T, TBase>(name, grp, "WeightedCuthillMcKeeOrdering")
			.add_constructor()
			.add_method("set_reverse", &T::set_reverse)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "WeightedCuthillMcKeeOrdering", tag);
	}

//	Boost Shortest Paths
	{
		typedef BoostShortestPathsOrdering<typename TAlgebra::matrix_type, ordering_container_type> T;
		typedef IOrderingAlgorithm<typename TAlgebra::matrix_type, ordering_container_type> TBase;
		string name = string("BoostShortestPathsOrdering").append(suffix);
		reg.add_class_<T, TBase>(name, grp, "BoostShortestPathsOrdering")
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "BoostShortestPathsOrdering", tag);
	}

//	Boost Cuthill McKee
	{
		typedef BoostCuthillMcKeeOrdering<typename TAlgebra::matrix_type, ordering_container_type> T;
		typedef IOrderingAlgorithm<typename TAlgebra::matrix_type, ordering_container_type> TBase;
		string name = string("BoostCuthillMcKeeOrdering").append(suffix);
		reg.add_class_<T, TBase>(name, grp, "BoostCuthillMcKeeOrdering")
			.add_constructor()
			.add_method("set_reverse", &T::set_reverse)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "BoostCuthillMcKeeOrdering", tag);
	}

//	Boost Minimum-Degree
	{
		typedef BoostMinimumDegreeOrdering<typename TAlgebra::matrix_type, ordering_container_type> T;
		typedef IOrderingAlgorithm<typename TAlgebra::matrix_type, ordering_container_type> TBase;
		string name = string("BoostMinimumDegreeOrdering").append(suffix);
		reg.add_class_<T, TBase>(name, grp, "BoostMinimumDegreeOrdering")
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "BoostMinimumDegreeOrdering", tag);
	}


//	Native Cuthill McKee
	{
		typedef NativeCuthillMcKeeOrdering<typename TAlgebra::matrix_type, ordering_container_type> T;
		typedef IOrderingAlgorithm<typename TAlgebra::matrix_type, ordering_container_type> TBase;
		string name = string("NativeCuthillMcKeeOrdering").append(suffix);
		reg.add_class_<T, TBase>(name, grp, "NativeCuthillMcKeeOrdering")
			.add_constructor()
			.add_method("set_reverse", &T::set_reverse)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NativeCuthillMcKeeOrdering", tag);
	}
}

}; // end Functionality

// end group ordering_bridge
/// \}

}// end Ordering

/// \addtogroup ordering_bridge
void RegisterBridge_AlgebraOrdering(Registry& reg, string grp)
{
	grp.append("/Algebra/Ordering");
	typedef Ordering::Functionality Functionality;

	try{
		RegisterAlgebraDependent<Functionality>(reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

} // namespace bridge
} // namespace ug
