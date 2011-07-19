// created by Andreas Vogel, Sebastian Reiter
// s.b.reiter@googlemail.com
// 10.02.2011 (m,d,y)

#include <iostream>
#include <sstream>
#include <vector>

#include "registry/registry.h"
#include "ug_bridge.h"

#include "common/profiler/profiler.h"
#include "lib_grid/lib_grid.h"

#include "lib_discretization/domain.h"
#include "lib_discretization/domain_util.h"

#include "lib_discretization/parallelization/domain_distribution.h"

#ifdef UG_PARALLEL
	#include "lib_grid/parallelization/load_balancing.h"
	#include "lib_grid/parallelization/parallelization.h"
#endif

using namespace std;

namespace ug{

template <typename TDomain>
static bool LoadDomain(TDomain& domain, const char* filename)
{
#ifdef UG_PARALLEL
	if(pcl::GetProcRank() != 0)
		return true;
#endif

	const char * p = strstr(filename, ".ugx");
	if(p == NULL)
	{
		UG_LOG("Currently only '.ugx' format supported for domains.\n");
		return false;
	}

	if(!LoadGridFromUGX(domain.get_grid(), domain.get_subset_handler(), filename))
	{
		UG_LOG("Cannot load grid.\n");
		return false;
	}

	return true;
}

template <typename TDomain>
static bool SaveDomain(TDomain& domain, const char* filename)
{
	const char * p = strstr(filename, ".ugx");
	if(p == NULL)
	{
		UG_LOG("Currently only '.ugx' format supported for domains.\n");
		return false;
	}

	return SaveGridToUGX(domain.get_grid(), domain.get_subset_handler(), filename);
}

template <typename TDomain>
static bool SavePartitionMap(PartitionMap& pmap, TDomain& domain,
							 const char* filename)
{
	const char * p = strstr(filename, ".ugx");
	if(p == NULL)
	{
		UG_LOG("Currently only '.ugx' format supported for partition-maps.\n");
		return false;
	}

	if(&domain.get_grid() != pmap.get_partition_handler().get_assigned_grid())
	{
		UG_LOG("WARNING in SavePartitionMap: The given partition map was not"
				" created for the given domain. Aborting...\n");
		return false;
	}

	return SavePartitionMapToFile(pmap, filename, domain.get_position_attachment());
}


template <typename TDomain>
static bool TestDomainInterfaces(TDomain* dom)
{
	#ifdef UG_PARALLEL
		return TestGridLayoutMap(dom->get_grid(),
					dom->get_distributed_grid_manager()->grid_layout_map());
	#endif
	return true;
}


namespace bridge{

template <typename TDomain>
static bool RegisterDomainInterface_(Registry& reg, const char* parentGroup)
{
	typedef TDomain domain_type;
	static const int dim = domain_type::dim;

//	the dimension suffix
	std::stringstream ss;	ss << dim << "d";
	std::string dimSuffix = ss.str();

//	the dimension tag
	std::string dimTag = "dim=";
	dimTag.append(dimSuffix);

//	get group string
	std::string grp = parentGroup;
	grp.append("/").append(dimSuffix);

//	Domain
	{
		std::string name = string("Domain").append(dimSuffix);
		reg.add_class_<domain_type>(name.c_str(), grp.c_str())
			.add_constructor()
			.add_method("get_subset_handler|hide=true", static_cast<MGSubsetHandler& (domain_type::*)()>(&domain_type::get_subset_handler))
			.add_method("get_grid|hide=true", static_cast<MultiGrid& (domain_type::*)()>(&domain_type::get_grid))
			.add_method("get_dim|hide=true", static_cast<int (domain_type::*)() const>(&domain_type::get_dim));

		reg.add_class_to_group(name.c_str(), "Domain", dimTag.c_str());
	}

// 	LoadDomain
	reg.add_function("LoadDomain", &LoadDomain<domain_type>, grp.c_str(),
					"Success", "Domain # Filename | load-dialog | endings=[\"ugx\"]; description=\"*.ugx-Files\" # Number Refinements",
					"Loads a domain", "No help");

//	SaveDomain
	reg.add_function("SaveDomain", &SaveDomain<domain_type>, grp.c_str(),
					"Success", "Domain # Filename|save-dialog",
					"Saves a domain", "No help");

//	SavePartitionMap
	reg.add_function("SavePartitionMap", &SavePartitionMap<domain_type>, grp.c_str(),
					"Success", "PartitionMap # Domain # Filename|save-dialog",
					"Saves a partition map", "No help");

//	DistributeDomain
	reg.add_function("DistributeDomain", static_cast<bool (*)(TDomain&)>(
					 &DistributeDomain<domain_type>), grp.c_str());

	reg.add_function("DistributeDomain", static_cast<bool (*)(TDomain&, PartitionMap&)>(
					 &DistributeDomain<domain_type>), grp.c_str());

//	todo: remove this
	{
		std::string name = string("DistributeDomain").append(dimSuffix);
		reg.add_function(name.c_str(), static_cast<bool (*)(TDomain&)>(
						 &DistributeDomain<domain_type>), grp.c_str());
	}

	reg.add_function("PartitionDomain_Bisection",
					 &PartitionDomain_Bisection<domain_type>, grp.c_str());

	reg.add_function("PartitionDomain_MetisKWay",
					 &PartitionDomain_MetisKWay<domain_type>, grp.c_str());

	reg.add_function("RedistributeDomain",
					 &RedistributeDomain<domain_type>, grp.c_str());

//	debugging
	reg.add_function("TestDomainInterfaces", &TestDomainInterfaces<domain_type>, grp.c_str());

	return true;
}

///	methods that are only available for 2d and 3d are registered here
template <typename TDomain>
static bool RegisterDomainInterface_2d_3d(Registry& reg, const char* parentGroup)
{
	typedef TDomain domain_type;
	static const int dim = domain_type::dim;

//	get group string
	std::stringstream grpSS; grpSS << parentGroup << "/" << dim << "d";
	std::string grp = grpSS.str();

	reg.add_function("PartitionDomain_RegularGrid",
					 &PartitionDomain_RegularGrid<domain_type>, grp.c_str());

	return true;
}

bool RegisterDomainInterface(Registry& reg, const char* parentGroup)
{
	bool bSuccess = true;

#ifdef UG_DIM_1
	bSuccess &= RegisterDomainInterface_<Domain<1, MultiGrid, MGSubsetHandler> >(reg, parentGroup);
#endif
#ifdef UG_DIM_2
	bSuccess &= RegisterDomainInterface_<Domain<2, MultiGrid, MGSubsetHandler> >(reg, parentGroup);
	bSuccess &= RegisterDomainInterface_2d_3d<Domain<2, MultiGrid, MGSubsetHandler> >(reg, parentGroup);
#endif
#ifdef UG_DIM_3
	bSuccess &= RegisterDomainInterface_<Domain<3, MultiGrid, MGSubsetHandler> >(reg, parentGroup);
	bSuccess &= RegisterDomainInterface_2d_3d<Domain<3, MultiGrid, MGSubsetHandler> >(reg, parentGroup);
#endif
	return bSuccess;
}

}// end of namespace
}// end of namespace
