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

#ifdef UG_PARALLEL
	#include "lib_grid/parallelization/load_balancing.h"
#endif

using namespace std;

namespace ug{

///	Used to describe how a domain shall be distributed in a parallel environment.
class PartitionMap{
	public:
		void clear()
		{
			m_targetProcs.clear();
			m_shPartitions.clear();
		}

		void assign_grid(Grid& grid)
		{
			if(&grid != m_shPartitions.get_assigned_grid())
				m_shPartitions.assign_grid(grid);
		}

		SubsetHandler& get_partition_handler()
		{return m_shPartitions;}

		void add_target_proc(int targetProcRank)
		{m_targetProcs.push_back(targetProcRank);}

		void add_target_procs(int first, int num)
		{
			for(int i = 0; i < num; ++i)
				add_target_proc(first + i);
		}

		size_t num_target_procs()
		{return m_targetProcs.size();}

		int get_target_proc(size_t index)
		{
			if(index < m_targetProcs.size())
				return m_targetProcs[index];
			UG_LOG("BAD INDEX in PartitionMap::get_target_proc: " << index);
			if(num_target_procs() > 0){
				UG_LOG("    Max valid index: " << num_target_procs() - 1 << endl);
			}
			else{
				UG_LOG("    No target processes available.\n");
			}
			return -1;
		}

		int* get_target_procs()
		{return &m_targetProcs.front();}

		std::vector<int>& get_target_proc_vec()
		{return m_targetProcs;}

	///	returns the index at which the given process lies. -1 if it doesn't exist.
		int find_target_proc(int procRank)
		{
			for(size_t i = 0; i < m_targetProcs.size(); ++i){
				if(m_targetProcs[i] == procRank)
					return i;
			}
			return -1;
		}

	private:
		std::vector<int>	m_targetProcs;
		SubsetHandler		m_shPartitions;
};



///	partitions a domain by repeatedly cutting it along the different axis
template <typename TDomain>
static bool PartitionDomain_Bisection(TDomain& domain, PartitionMap& partitionMap,
									  int firstAxisToCut)
{
	MultiGrid& mg = domain.get_grid();
	partitionMap.assign_grid(mg);
	#ifdef UG_PARALLEL
		if(mg.num<Volume>() > 0)
			PartitionElementsByRepeatedIntersection<Volume, TDomain::dim>(
												partitionMap.get_partition_handler(),
												mg, mg.num_levels() - 1,
												partitionMap.num_target_procs(),
												domain.get_position_attachment(),
												firstAxisToCut);
		else if(mg.num<Face>() > 0)
			PartitionElementsByRepeatedIntersection<Face, TDomain::dim>(
												partitionMap.get_partition_handler(),
												mg, mg.num_levels() - 1,
												partitionMap.num_target_procs(),
												domain.get_position_attachment(),
												firstAxisToCut);
		else if(mg.num<EdgeBase>() > 0)
			PartitionElementsByRepeatedIntersection<EdgeBase, TDomain::dim>(
												partitionMap.get_partition_handler(),
												mg, mg.num_levels() - 1,
												partitionMap.num_target_procs(),
												domain.get_position_attachment(),
												firstAxisToCut);
		else if(mg.num<VertexBase>() > 0)
			PartitionElementsByRepeatedIntersection<VertexBase, TDomain::dim>(
												partitionMap.get_partition_handler(),
												mg, mg.num_levels() - 1,
												partitionMap.num_target_procs(),
												domain.get_position_attachment(),
												firstAxisToCut);
		else{
			LOG("partitioning could not be performed - "
				<< "grid doesn't contain any elements!\n");
			return false;
		}
	#endif

	return true;
}

///	partitions a domain by sorting all elements into a regular grid
template <typename TDomain>
static bool PartitionDomain_RegularGrid(TDomain& domain, PartitionMap& partitionMap,
										int numCellsX, int numCellsY)
{
//	prepare the partition map and a vertex position attachment accessor
	MultiGrid& mg = domain.get_grid();
	partitionMap.assign_grid(mg);

	#ifdef UG_PARALLEL
	//	a distributed grid manager is required
		if(!domain.get_distributed_grid_manager()){
			UG_LOG("A distributed grid manager is required in the given domain.\n");
			return false;
		}

		typedef typename TDomain::position_attachment_type TAPos;
		Grid::AttachmentAccessor<VertexBase, TAPos> aaPos(mg,
												domain.get_position_attachment());

	//	this callback allows us to only distribute surface elements, which are no ghosts
		IsRegularSurfaceElem cbConsiderElem(*domain.get_distributed_grid_manager());

	//	we need a process to which elements which are not considered will be send.
	//	Those elements should stay on the current process.
		int localProc = 0;
		localProc = pcl::GetProcRank();

		int bucketSubset = partitionMap.find_target_proc(localProc);
		if(bucketSubset == -1)
			bucketSubset = (int)partitionMap.num_target_procs();

	//	partition the grid
		if(mg.num<Volume>() > 0)
			PartitionElements_RegularGrid<Volume>(
										partitionMap.get_partition_handler(),
										mg.begin<Volume>(), mg.end<Volume>(),
										numCellsX, numCellsY, aaPos,
										cbConsiderElem, bucketSubset);
		else if(mg.num<Face>() > 0)
			PartitionElements_RegularGrid<Face>(
										partitionMap.get_partition_handler(),
										mg.begin<Face>(), mg.end<Face>(),
										numCellsX, numCellsY, aaPos,
										cbConsiderElem, bucketSubset);
		else if(mg.num<EdgeBase>() > 0)
			PartitionElements_RegularGrid<EdgeBase>(
										partitionMap.get_partition_handler(),
										mg.begin<EdgeBase>(), mg.end<EdgeBase>(),
										numCellsX, numCellsY, aaPos,
										cbConsiderElem, bucketSubset);
		else if(mg.num<VertexBase>() > 0)
			PartitionElements_RegularGrid<VertexBase>(
										partitionMap.get_partition_handler(),
										mg.begin<VertexBase>(), mg.end<VertexBase>(),
										numCellsX, numCellsY, aaPos,
										cbConsiderElem, bucketSubset);
		else{
			LOG("partitioning could not be performed - "
				<< "grid doesn't contain any elements!\n");
			return false;
		}

	//	if elements have been assigned to bucketProc, then we have to make sure,
	//	that it is also present in the process-map
		if(!partitionMap.get_partition_handler().empty(bucketSubset)){
			if(bucketSubset >= (int)partitionMap.num_target_procs())
				partitionMap.add_target_proc(localProc);
		}
	#endif

	return true;
}

template <typename TDomain>
static bool RedistributeDomain(TDomain& domainOut,
							   PartitionMap& partitionMap,
							   bool createVerticalInterfaces)
{
//todo	Use a process-communicator to restrict communication

	typedef typename TDomain::distributed_grid_manager_type distributed_grid_manager_type;
	typedef typename TDomain::position_attachment_type	position_attachment_type;
//	make sure that the input is fine
	typename TDomain::grid_type& grid = domainOut.get_grid();
	SubsetHandler& shPart = partitionMap.get_partition_handler();

	if(shPart.get_assigned_grid() != &grid){
		partitionMap.assign_grid(grid);
	}

//	used to check whether all processes are correctly prepared for redistribution
	bool performDistribution = true;

//	make sure that the number of subsets and target processes match
	const int numSubs = partitionMap.get_partition_handler().num_subsets();
	const int numTargetProcs = (int)partitionMap.num_target_procs();
	if(numSubs > numTargetProcs){
		UG_LOG("ERROR in RedistributeDomain: More partitions than target processes.\n");
		performDistribution = false;
	}
	else if(numSubs < numTargetProcs){
		UG_LOG("ERROR in RedistributeDomain: More target processes than partitions.\n");
		performDistribution = false;
	}

//todo:	check whether all target-processes in partitionMap are in the valid range.

#ifdef UG_PARALLEL
//	make sure that manager exists
	distributed_grid_manager_type* pDistGridMgr = domainOut.get_distributed_grid_manager();
	if(!pDistGridMgr)
	{
		UG_LOG("ERROR in RedistibuteDomain: Domain has to feature a Distributed Grid Manager.\n");
		performDistribution = false;
	}

//todo	Use a process-communicator to restrict communication
	if(!pcl::AllProcsTrue(performDistribution))
		return false;

	distributed_grid_manager_type& distGridMgr = *pDistGridMgr;

//	data serialization
	GeomObjAttachmentSerializer<VertexBase, position_attachment_type>
		posSerializer(grid, domainOut.get_position_attachment());
	SubsetHandlerSerializer shSerializer(domainOut.get_subset_handler());

	GridDataSerializationHandler serializer;
	serializer.add(&posSerializer);
	serializer.add(&shSerializer);

//	now call redistribution
	RedistributeGrid(distGridMgr, shPart, serializer, serializer,
					 createVerticalInterfaces, &partitionMap.get_target_proc_vec());

#endif

//	in the serial case there's nothing to do.
	return true;
}



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

	reg.add_function("PartitionDomain_Bisection",
					 &PartitionDomain_Bisection<domain_type>, grp.c_str());

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

	reg.add_class_<PartitionMap>("PartitionMap", "ug4")
		.add_constructor()
		.add_method("clear", &PartitionMap::clear)
		.add_method("get_partition_handler", &PartitionMap::get_partition_handler)
		.add_method("add_target_proc", &PartitionMap::add_target_proc)
		.add_method("add_target_procs", &PartitionMap::add_target_procs)
		.add_method("num_target_procs", &PartitionMap::num_target_procs)
		.add_method("get_target_proc", &PartitionMap::get_target_proc);


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
