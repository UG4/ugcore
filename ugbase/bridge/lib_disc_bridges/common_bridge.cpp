/*
 * common_bridge.cpp
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

//	Elem Discs
	{
		string elemGrp = grp; elemGrp.append("/ElemDisc");
		typedef IElemDisc T;
		reg.add_class_<T>("IElemDisc", grp);
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
