/*
 * lib_disc_bridge_algebra_independent.cpp
 *
 *  Created on: 14.02.2011
 *      Author: andreasvogel
 */

#include <string>

// include brigde
#include "../bridge.h"

// includes of lib_discretization
#include "lib_disc/common/function_group.h"
#include "lib_disc/dof_manager/function_pattern.h"
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"

using namespace std;

namespace ug
{
namespace bridge
{

bool RegisterLibDisc_Common(Registry& reg, string parentGroup)
{
//	get group string
	string grp = parentGroup; grp.append("/Discretization");

	try
	{
//	FunctionPattern
	{
		typedef FunctionPattern T;
		string elemGrp = grp; elemGrp.append("/ApproximationSpace");
		reg.add_class_<T>("FunctionPattern", grp)
			.add_method("clear", &T::clear)
			.add_method("add_fct", static_cast<void (T::*)(const char*, const char*, int, const char*)>(&T::add_fct),
						"", "Name # Type|selection|value=[\"Lagrange\",\"DG\"] # Order # Subsets", "Adds a function to the Function Pattern",
						"currently no help available")
			.add_method("add_fct", static_cast<void (T::*)(const char*, const char*, int)>(&T::add_fct),
						"", "Name # Type|selection|value=[\"Lagrange\",\"DG\"] # Order", "Adds a function to the Function Pattern",
						"currently no help available");
	}

//	Elem Discs
	{
		string elemGrp = grp; elemGrp.append("/ElemDisc");
		typedef IElemDisc T;
		reg.add_class_<T>("IElemDisc", grp)
			.add_method("set_functions", &T::set_functions,
						"", "Functions (sep. by ',')")
			.add_method("set_subsets",  &T::set_subsets,
						"", "Subsets (sep. by ',')");
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
		.add_method("get_num_subdomains",          &T::get_num_subdomains);
	}
#endif

	}
	catch(UG_REGISTRY_ERROR_RegistrationFailed ex)
	{
		UG_LOG("### ERROR in RegisterLibDisc_Common: "
				"Registration failed (using name " << ex.name << ").\n");
		return false;
	}

	return true;
}

} // end namespace bridge
} // end namespace ug
