/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
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

#include <string>

// include brigde
#include "bridge/bridge.h"
#include "bridge/util.h"

// includes of lib_discretization
#include "lib_disc/common/function_group.h"
#include "lib_disc/common/multi_index.h"
#include "lib_disc/dof_manager/function_pattern.h"
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"

using namespace std;

namespace ug
{
namespace bridge
{

/**
 * \defgroup disccommon_bridge Common Discretization Bridge
 * \ingroup disc_bridge
 * \{
 */
template <typename TRegistry=Registry>
void RegisterBridge_DiscCommon_(TRegistry& reg, string parentGroup)
{
//	get group string
	string grp = parentGroup; grp.append("/Discretization");

	try
	{
		//  MultiIndex2
		{
			typedef MultiIndex<2, size_t> T;
			string name = string("MultiIndex2");
			reg.template add_class_<T>(name, grp);
		}

#ifdef UG_PARALLEL
	//	IDomainDecompositionInfo, StandardDomainDecompositionInfo
		{
		typedef pcl::IDomainDecompositionInfo Tbase;
		reg.template add_class_<Tbase>("IDomainDecompositionInfo", grp);

		typedef pcl::StandardDomainDecompositionInfo T;
		reg.template add_class_<T, Tbase>("StandardDomainDecompositionInfo", grp)
			.add_constructor()
			.add_method("map_proc_id_to_subdomain_id", &T::map_proc_id_to_subdomain_id)
			.add_method("set_num_subdomains",          &T::set_num_subdomains)
			.add_method("get_num_subdomains",          &T::get_num_subdomains)
			.add_method("set_num_spatial_dimensions",  &T::set_num_spatial_dimensions)
			.add_method("get_num_spatial_dimensions",  &T::get_num_spatial_dimensions)
			.set_construct_as_smart_pointer(true);
		}
#endif

	}
	UG_REGISTRY_CATCH_THROW(grp);
}

// end group disccommon_bridge
/// \}

} // end namespace bridge

UG_REGISTRY_DEFINE(RegisterBridge_DiscCommon);
} // end namespace ug
