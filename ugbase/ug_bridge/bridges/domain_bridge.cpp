// created by Andreas Vogel, Sebastian Reiter
// s.b.reiter@googlemail.com
// 10.02.2011 (m,d,y)

#include <iostream>
#include <sstream>
#include <vector>

#include "../registry.h"
#include "../ug_bridge.h"

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
static bool LoadDomain(TDomain& domain, const char* filename, int numRefs)
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

	if(numRefs <= 0) return true;

	GlobalMultiGridRefiner refiner;
	refiner.assign_grid(domain.get_grid());
	for(int i = 0; i < numRefs; ++i)
	{
		refiner.refine();
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
static void GlobalRefineParallelDomain(TDomain& domain)
{
#ifdef UG_PARALLEL
//	get distributed grid manager
	typedef typename TDomain::distributed_grid_manager_type distributed_grid_manager_type;
	distributed_grid_manager_type* pDistGridMgr = domain.get_distributed_grid_manager();

//	check that manager exists
	if(!pDistGridMgr)
	{
		UG_LOG("GlobalRefineParallelDomain: Cannot find Distributed Grid Manager.\n");
		throw(int(1));
	}
	distributed_grid_manager_type& distGridMgr = *pDistGridMgr;

//	create Refiner
	ParallelGlobalRefiner_MultiGrid refiner(distGridMgr);
#else
	GlobalMultiGridRefiner refiner;
	refiner.assign_grid(domain.get_grid());
#endif

//	perform refinement.
	refiner.refine();
}


template <typename TDomain>
static IRefiner* GlobalDomainRefiner(TDomain* dom)
{
//todo: support normal grids, too!
	#ifdef UG_PARALLEL
		if(pcl::GetNumProcesses() > 1){
			return new ParallelGlobalRefiner_MultiGrid(*dom->get_distributed_grid_manager());
		}
	#endif

	return new GlobalMultiGridRefiner(dom->get_grid());
}


template <typename TDomain>
static IRefiner* HangingNodeDomainRefiner(TDomain* dom)
{
//todo: support normal grids, too!
	#ifdef UG_PARALLEL
		if(pcl::GetNumProcesses() > 1){
			return new ParallelHangingNodeRefiner_MultiGrid(*dom->get_distributed_grid_manager());
		}
	#endif

	return new HangingNodeRefiner_MultiGrid(dom->get_grid());
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

//	get group string
	std::stringstream grpSS; grpSS << parentGroup << "/" << dim << "d";
	std::string grp = grpSS.str();

//	Domain
	{
		std::stringstream ss; ss << "Domain" << dim << "d";
		reg.add_class_<domain_type>(ss.str().c_str(), grp.c_str())
			.add_constructor()
			.add_method("get_subset_handler|hide=true", (MGSubsetHandler& (domain_type::*)()) &domain_type::get_subset_handler)
			.add_method("get_grid|hide=true", (MultiGrid& (domain_type::*)()) &domain_type::get_grid)
			.add_method("get_dim|hide=true", (int (domain_type::*)()) &domain_type::get_dim);
	}

// 	LoadDomain
	reg.add_function("LoadDomain", &LoadDomain<domain_type>, grp.c_str(),
					"Success", "Domain # Filename | load-dialog | endings=[\"ugx\"]; description=\"*.ugx-Files\" # Number Refinements",
					"Loads a domain", "No help");
//	todo: remove this
	{
		std::stringstream ss; ss << "LoadDomain" << dim << "d";
		reg.add_function(ss.str().c_str(), &LoadDomain<domain_type>, grp.c_str(),
						"Success", "Domain # Filename | load-dialog | endings=[\"ugx\"]; description=\"*.ugx-Files\" # Number Refinements",
						"Loads a domain", "No help");
	}

//	SaveDomain
	reg.add_function("SaveDomain", &SaveDomain<domain_type>, grp.c_str(),
					"Success", "Domain # Filename|save-dialog",
					"Saves a domain", "No help");

//	SavePartitionMap
	reg.add_function("SavePartitionMap", &SavePartitionMap<domain_type>, grp.c_str(),
					"Success", "PartitionMap # Domain # Filename|save-dialog",
					"Saves a partition map", "No help");

//	todo: remove this
	{
		std::stringstream ss; ss << "SaveDomain" << dim << "d";
		reg.add_function(ss.str().c_str(), &SaveDomain<domain_type>, grp.c_str(),
						"Success", "Domain # Filename|save-dialog",
						"Saves a domain", "No help");
	}

//	DistributeDomain
	reg.add_function("DistributeDomain", (bool (*)(TDomain&))
					 &DistributeDomain<domain_type>, grp.c_str());

	reg.add_function("DistributeDomain", (bool (*)(TDomain&, PartitionMap&))
					 &DistributeDomain<domain_type>, grp.c_str());

//	todo: remove this
	{
		std::stringstream ss; ss << "DistributeDomain" << dim << "d";
		reg.add_function(ss.str().c_str(), (bool (*)(TDomain&))
						 &DistributeDomain<domain_type>, grp.c_str());
	}

	reg.add_function("PartitionDomain_Bisection",
					 &PartitionDomain_Bisection<domain_type>, grp.c_str());

	reg.add_function("PartitionDomain_MetisKWay",
					 &PartitionDomain_MetisKWay<domain_type>, grp.c_str());

	reg.add_function("RedistributeDomain",
					 &RedistributeDomain<domain_type>, grp.c_str());

//	GlobalRefineParallelDomain
//	todo: remove this
	{
		std::stringstream ss; ss << "GlobalRefineParallelDomain" << dim << "d";
		reg.add_function(ss.str().c_str(), &GlobalRefineParallelDomain<domain_type>, grp.c_str());
	}

//	refiner registration
	reg.add_function("GlobalDomainRefiner", &GlobalDomainRefiner<domain_type>, grp.c_str());
	reg.add_function("HangingNodeDomainRefiner", &HangingNodeDomainRefiner<domain_type>, grp.c_str());

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
