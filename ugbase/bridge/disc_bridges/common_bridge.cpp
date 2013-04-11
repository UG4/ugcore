/*
 * common_bridge.cpp
 *
 *  Created on: 14.02.2011
 *      Author: andreasvogel
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

void RegisterBridge_DiscCommon(Registry& reg, string parentGroup)
{
//	get group string
	string grp = parentGroup; grp.append("/Discretization");

	try
	{
		//  MultiIndex2
		{
			typedef MultiIndex<2, size_t> T;
			string name = string("MultiIndex2");
			reg.add_class_<T>(name, grp);
		}

#ifdef UG_PARALLEL
	//	IDomainDecompositionInfo, StandardDomainDecompositionInfo
		{
		typedef pcl::IDomainDecompositionInfo Tbase;
		reg.add_class_<Tbase>("IDomainDecompositionInfo", grp);
		typedef pcl::StandardDomainDecompositionInfo T;
		reg.add_class_<T, Tbase>("StandardDomainDecompositionInfo", grp)
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
} // end namespace ug
