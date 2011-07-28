/*
 * lib_disc_bridge_algebra_independent.cpp
 *
 *  Created on: 14.02.2011
 *      Author: andreasvogel
 */

#include <string>

// include brigde
#include "../ug_bridge.h"

// includes of lib_discretization
#include "lib_discretization/common/function_group.h"
#include "lib_discretization/dof_manager/function_pattern.h"
#include "lib_discretization/spatial_discretization/elem_disc/elem_disc_interface.h"

using namespace std;

namespace ug
{
namespace bridge
{

bool RegisterStaticLibDiscInterface(Registry& reg, const char* parentGroup)
{
	try
	{
	//	get group string
		string grp = parentGroup; grp.append("/Discretization");

	//	FunctionGroup
		{
			reg.add_class_<FunctionGroup>("FunctionGroup", grp)
				.add_constructor()
				.add_method("clear", &FunctionGroup::clear)
				.add_method("set_function_pattern", &FunctionGroup::set_function_pattern)
				.add_method("add_function", static_cast<bool (FunctionGroup::*)(const char*)>(&FunctionGroup::add));
		}

	//	FunctionPattern
		{
			typedef FunctionPattern T;
			reg.add_class_<T>("FunctionPattern", grp)
				.add_method("clear", &T::clear)
				.add_method("add_fct", static_cast<bool (T::*)(const char*, const char*, int, const char*)>(&T::add_fct),
				            "Success", "Name # Type|selection|value=[\"Lagrange\",\"DG\"] # Order # Subsets", "Adds a function to the Function Pattern",
				            "currently no help available")
				.add_method("add_fct", static_cast<bool (T::*)(const char*, const char*, int)>(&T::add_fct),
				            "Success", "Name # Type|selection|value=[\"Lagrange\",\"DG\"] # Order", "Adds a function to the Function Pattern",
				            "currently no help available");
		}

	//	Elem Discs
		{
		//	Base class
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
		UG_LOG("### ERROR in RegisterLibDiscretizationInterface: "
				"Registration failed (using name " << ex.name << ").\n");
		return false;
	}

	return true;
}



} // end namespace bridge
} // end namespace ug
