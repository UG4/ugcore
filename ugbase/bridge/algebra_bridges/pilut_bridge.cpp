/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
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

#if 0

// preconditioner
#include "lib_algebra/lib_algebra.h"
#include "lib_algebra/operator/preconditioner/pilut.h"
using namespace std;

namespace ug{
namespace bridge{
namespace Preconditioner{

/**
 * \defgroup precond_bridge Preconditioner Bridge
 * \ingroup algebra_bridge
 * \{
 */

/**
 * Class exporting the functionality. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */
struct Functionality2
{

/**
 * Function called for the registration of Algebra dependent parts.
 * All Functions and Classes depending on Algebra
 * are to be placed here when registering. The method is called for all
 * available Algebra types, based on the current build options.
 *
 * @param reg				registry
 * @param grp				group for sorting of functionality
 */
template <typename TAlgebra>
static void Algebra(Registry& reg, string grp)
{
	string suffix = GetAlgebraSuffix<TAlgebra>();
	string tag = GetAlgebraTag<TAlgebra>();

	using vector_type = typename TAlgebra::vector_type;
	using matrix_type = typename TAlgebra::matrix_type;



//	PILU Threshold
	{
		using T = PILUTPreconditioner<TAlgebra>;
		using TBase = IPreconditioner<TAlgebra>;
		string name = string("PILUT").append(suffix);
		reg.add_class_<T,TBase>(name, grp, "Incomplete LU Decomposition with threshold")
			.add_constructor()
			.template add_constructor<void (*)(number)>("threshold parameter")
			.add_method("set_threshold", &T::set_threshold,
						"", "threshold", "sets threshold of incomplete LU factorisation")
			.add_method("set_info", &T::set_info,
						"", "info", "sets storage information output")
			.add_method("set_show_progress", &T::set_show_progress,
						"", "onoff", "switches the progress indicator on/off")
			//.add_method("set_group_size", &T::set_group_size,"", "groupSize", "sets group size")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "PILUT", tag);
	}
}
	

}; // end Functionality

// end group precond_bridge
/// \}

}// end Preconditioner

/// \addtogroup precond_bridge
void RegisterBridge_PILUT(Registry& reg, string grp)
{
	grp.append("/Algebra/Preconditioner");
	using Functionality2 = Preconditioner::Functionality2;

	try{
		RegisterAlgebraDependent<Functionality2>(reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

} // namespace bridge
} // namespace ug

#endif

namespace ug{
namespace bridge{
/// \addtogroup precond_bridge
void RegisterBridge_PILUT(Registry& reg, std::string grp)
{

}

} // namespace bridge
} // namespace ug
