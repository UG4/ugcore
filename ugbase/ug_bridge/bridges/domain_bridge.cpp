// created by Andreas Vogel, Sebastian Reiter
// s.b.reiter@googlemail.com
// 10.02.2011 (m,d,y)

#include <iostream>
#include <sstream>

#include "../registry.h"
#include "../ug_bridge.h"

#include "common/profiler/profiler.h"
#include "lib_grid/lib_grid.h"

#include "lib_discretization/domain.h"
#include "lib_discretization/domain_util.h"

using namespace std;

namespace ug{
namespace bridge{

template <typename TDomain>
static bool DistributeDomain(TDomain& domainOut)
{
#ifdef UG_PARALLEL
//	typedefs
	typedef typename TDomain::subset_handler_type subset_handler_type;
	typedef typename TDomain::distributed_grid_manager_type distributed_grid_manager_type;

//	get distributed grid manager
	distributed_grid_manager_type* pDistGridMgr = domainOut.get_distributed_grid_manager();

//	check that manager exists
	if(!pDistGridMgr)
	{
		UG_LOG("DistibuteDomain: Cannot find Distributed Grid Manager.\n");
		return false;
	}
	distributed_grid_manager_type& distGridMgrOut = *pDistGridMgr;

//	get subset handler
	subset_handler_type& sh = domainOut.get_subset_handler();

//	get number of processes
	const int numProcs = pcl::GetNumProcesses();
	if(numProcs == 1) return true;

//	check, that grid is a multigrid
	MultiGrid* pMG = dynamic_cast<MultiGrid*>(distGridMgrOut.get_assigned_grid());
	if(pMG == NULL)
	{
		UG_LOG("DistibuteDomain: MultiGrid-Domain required in current implementation.\n");
		return false;
	}
	MultiGrid& mg = *pMG;

//	get Grid Layout
//	GridLayoutMap& glmOut = distGridMgrOut.grid_layout_map();

//	make sure that each grid has a position attachment - even if no data
//	will be received.
	typedef typename TDomain::position_attachment_type position_attachment_type;
	position_attachment_type& domPosition = domainOut.get_position_attachment();
	bool tmpPosAttachment = false;
	if(!mg.has_vertex_attachment(aPosition))
	{
	// convert to 3d positions (FVGeometry depends on PositionCoordinates)
       mg.attach_to_vertices(aPosition);
       ConvertMathVectorAttachmentValues<VertexBase>(mg, domPosition, aPosition);

       UG_LOG("DistributeDomain: temporarily adding Position Attachment.\n");
       tmpPosAttachment = true;
	}

//	AdjustFunctions
//	FuncAdjustGrid funcAdjustGrid = DefaultAdjustGrid;
	FuncAdjustGrid funcAdjustGrid = AdjustGrid_AutoAssignSubsetsAndRefine(0,1,0);
	FuncPartitionGrid funcPartitionGrid = PartitionGrid_Bisection;

//	perform distribution
	AdjustAndDistributeGrid(distGridMgrOut, sh, mg, sh, numProcs, true,
							funcAdjustGrid, funcPartitionGrid);

	if(tmpPosAttachment)
	{
	// convert to 3d positions (FVGeometry depends on PositionCoordinates)
       ConvertMathVectorAttachmentValues<VertexBase>(mg, aPosition, domPosition);
       mg.detach_from_vertices(aPosition);

       UG_LOG("DistributeDomain: removing temporary Position Attachment.\n");
 	}

#endif

//	in serial case: do nothing
	return true;
}

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
/*
Ich habe mir erlaubt das hier auszukommentieren. Ich denke man sollte den
Namen der Domain komplett im Skript spezifizieren. Das geht auch ohne weiteres:
domName = "mydom_" .. GetProcessRank() .. ".ugx"
Das kann man besser über eine util Funktion lösen.
Grüße,
Sebastian


	std::string name(filename);

#ifdef UG_PARALLEL
//	search for ending
	size_t found = name.find_first_of(".");

//	remove endings
	name.resize(found);

//	add new ending, containing process number
	int rank = pcl::GetProcRank();
	char ext[20];
	sprintf(ext, "_p%04d.ugx", rank);
	name.append(ext);
#endif

	return SaveGridToUGX(domain.get_grid(), domain.get_subset_handler(), name.c_str());
*/
	return SaveGridToUGX(domain.get_grid(), domain.get_subset_handler(), filename);
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
	ParallelGlobalMultiGridRefiner refiner(distGridMgr);
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
			return new ParallelGlobalMultiGridRefiner(*dom->get_distributed_grid_manager());
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
static bool RegisterDomainInterface_(Registry& reg, const char* parentGroup)
{
	typedef TDomain domain_type;
	static const int dim = domain_type::dim;

//	get group string
	std::stringstream grpSS; grpSS << parentGroup << "/" << dim << "d";
	std::string grp = grpSS.str();

	domain_type dom;

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
//	todo: remove this
	{
		std::stringstream ss; ss << "SaveDomain" << dim << "d";
		reg.add_function(ss.str().c_str(), &SaveDomain<domain_type>, grp.c_str(),
						"Success", "Domain # Filename|save-dialog",
						"Saves a domain", "No help");
	}

//	DistributeDomain
	reg.add_function("DistributeDomain", &DistributeDomain<domain_type>, grp.c_str());
//	todo: remove this
	{
		std::stringstream ss; ss << "DistributeDomain" << dim << "d";
		reg.add_function(ss.str().c_str(), &DistributeDomain<domain_type>, grp.c_str());
	}

//	GlobalRefineParallelDomain
//	todo: remove this
	{
		std::stringstream ss; ss << "GlobalRefineParallelDomain" << dim << "d";
		reg.add_function(ss.str().c_str(), &GlobalRefineParallelDomain<domain_type>, grp.c_str());
	}

//	refiner registration
	reg.add_function("GlobalDomainRefiner", &GlobalDomainRefiner<domain_type>, grp.c_str());
	reg.add_function("HangingNodeDomainRefiner", &HangingNodeDomainRefiner<domain_type>, grp.c_str());

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
#endif
#ifdef UG_DIM_3
	bSuccess &= RegisterDomainInterface_<Domain<3, MultiGrid, MGSubsetHandler> >(reg, parentGroup);
#endif
	return bSuccess;
}

}// end of namespace
}// end of namespace
