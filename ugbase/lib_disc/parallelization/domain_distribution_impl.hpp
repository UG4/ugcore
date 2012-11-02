// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 25.05.2011 (m,d,y)

#ifndef __H__UG__domain_distribution_impl__
#define __H__UG__domain_distribution_impl__

#include "domain_distribution.h"
#include "lib_grid/algorithms/attachment_util.h"
#include "lib_grid/parallelization/load_balancing.h"

#ifdef UG_PARALLEL
	#include "pcl/pcl.h"
	#include "lib_grid/parallelization/parallelization.h"
#endif


namespace ug
{

///	partitions a domain by repeatedly cutting it along the different axis
template <typename TDomain>
static bool PartitionDomain_Bisection(TDomain& domain, PartitionMap& partitionMap,
									  int firstAxisToCut)
{
	PROFILE_FUNC_GROUP("parallelization")
	SmartPtr<MultiGrid> pMG = domain.grid();
	partitionMap.assign_grid(*pMG);
	#ifdef UG_PARALLEL
	//	we need a process to which elements which are not considered will be send.
	//	Those elements should stay on the current process.
		int localProc = 0;
		localProc = pcl::GetProcRank();

		int bucketSubset = partitionMap.find_target_proc(localProc);
		if(bucketSubset == -1)
			bucketSubset = (int)partitionMap.num_target_procs();

		if(pMG->num<Volume>() > 0){
			partitionMap.get_partition_handler().assign_subset(
							pMG->begin<Volume>(), pMG->end<Volume>(), bucketSubset);

			PartitionElementsByRepeatedIntersection<Volume, TDomain::dim>(
												partitionMap.get_partition_handler(),
												*pMG, pMG->num_levels() - 1,
												partitionMap.num_target_procs(),
												domain.position_attachment(),
												firstAxisToCut);
		}
		else if(pMG->num<Face>() > 0){
			partitionMap.get_partition_handler().assign_subset(
							pMG->begin<Face>(), pMG->end<Face>(), bucketSubset);

			PartitionElementsByRepeatedIntersection<Face, TDomain::dim>(
												partitionMap.get_partition_handler(),
												*pMG, pMG->num_levels() - 1,
												partitionMap.num_target_procs(),
												domain.position_attachment(),
												firstAxisToCut);
		}
		else if(pMG->num<EdgeBase>() > 0){
			partitionMap.get_partition_handler().assign_subset(
							pMG->begin<EdgeBase>(), pMG->end<EdgeBase>(), bucketSubset);

			PartitionElementsByRepeatedIntersection<EdgeBase, TDomain::dim>(
												partitionMap.get_partition_handler(),
												*pMG, pMG->num_levels() - 1,
												partitionMap.num_target_procs(),
												domain.position_attachment(),
												firstAxisToCut);
		}
		else if(pMG->num<VertexBase>() > 0){
			partitionMap.get_partition_handler().assign_subset(
							pMG->begin<VertexBase>(), pMG->end<VertexBase>(), bucketSubset);

			PartitionElementsByRepeatedIntersection<VertexBase, TDomain::dim>(
												partitionMap.get_partition_handler(),
												*pMG, pMG->num_levels() - 1,
												partitionMap.num_target_procs(),
												domain.position_attachment(),
												firstAxisToCut);
		}
		else{
			LOG("partitioning could not be performed - "
				<< "grid doesn't contain any elements!\n");
			return false;
		}

		return true;
	#endif

//	Assign all elements to partition 0
	UG_LOG("WARNING: Serial fallback implementation of PartitionDomain_Bisection is used.\n");
	partitionMap.get_partition_handler().assign_subset(pMG->begin<VertexBase>(),
													   pMG->end<VertexBase>(), 0);
	partitionMap.get_partition_handler().assign_subset(pMG->begin<EdgeBase>(),
													   pMG->end<EdgeBase>(), 0);
	partitionMap.get_partition_handler().assign_subset(pMG->begin<Face>(),
													   pMG->end<Face>(), 0);
	partitionMap.get_partition_handler().assign_subset(pMG->begin<Volume>(),
													   pMG->end<Volume>(), 0);
	return true;
}

///	partitions a domain by sorting all elements into a regular grid
template <typename TDomain>
static bool PartitionDomain_RegularGrid(TDomain& domain, PartitionMap& partitionMap,
										int numCellsX, int numCellsY,
										bool surfaceOnly)
{
	PROFILE_FUNC_GROUP("parallelization")
//	prepare the partition map and a vertex position attachment accessor
	SmartPtr<MultiGrid> pMG = domain.grid();
	partitionMap.assign_grid(*pMG);

	#ifdef UG_PARALLEL
	//	a distributed grid manager is required
		if(!domain.distributed_grid_manager()){
			UG_LOG("A distributed grid manager is required in the given domain.\n");
			return false;
		}

		typedef typename TDomain::position_attachment_type TAPos;
		Grid::AttachmentAccessor<VertexBase, TAPos> aaPos(*pMG,
												domain.position_attachment());

	//	this callback allows us to only distribute surface elements, which are no ghosts
		IsRegularSurfaceElem cbConsiderElem(*domain.distributed_grid_manager());

	//	we need a process to which elements which are not considered will be send.
	//	Those elements should stay on the current process.
		int localProc = 0;
		localProc = pcl::GetProcRank();

		int bucketSubset = partitionMap.find_target_proc(localProc);
		if(bucketSubset == -1)
			bucketSubset = (int)partitionMap.num_target_procs();

	//	partition the grid
		if(pMG->num<Volume>() > 0){
			if(!surfaceOnly)
				PartitionElements_RegularGrid<Volume>(
											partitionMap.get_partition_handler(),
											pMG->begin<Volume>(), pMG->end<Volume>(),
											numCellsX, numCellsY, aaPos,
											Grid::volume_traits::cb_consider_all, bucketSubset);
			else
				PartitionElements_RegularGrid<Volume>(
											partitionMap.get_partition_handler(),
											pMG->begin<Volume>(), pMG->end<Volume>(),
											numCellsX, numCellsY, aaPos,
											cbConsiderElem, bucketSubset);
		}
		else if(pMG->num<Face>() > 0){
			if(!surfaceOnly)
				PartitionElements_RegularGrid<Face>(
											partitionMap.get_partition_handler(),
											pMG->begin<Face>(), pMG->end<Face>(),
											numCellsX, numCellsY, aaPos,
											Grid::face_traits::cb_consider_all, bucketSubset);
			else
				PartitionElements_RegularGrid<Face>(
											partitionMap.get_partition_handler(),
											pMG->begin<Face>(), pMG->end<Face>(),
											numCellsX, numCellsY, aaPos,
											cbConsiderElem, bucketSubset);
		}
		else if(pMG->num<EdgeBase>() > 0){
			if(!surfaceOnly)
				PartitionElements_RegularGrid<EdgeBase>(
											partitionMap.get_partition_handler(),
											pMG->begin<EdgeBase>(), pMG->end<EdgeBase>(),
											numCellsX, numCellsY, aaPos,
											Grid::edge_traits::cb_consider_all, bucketSubset);
			else
				PartitionElements_RegularGrid<EdgeBase>(
											partitionMap.get_partition_handler(),
											pMG->begin<EdgeBase>(), pMG->end<EdgeBase>(),
											numCellsX, numCellsY, aaPos,
											cbConsiderElem, bucketSubset);
		}
		else if(pMG->num<VertexBase>() > 0){
			if(!surfaceOnly)
				PartitionElements_RegularGrid<VertexBase>(
											partitionMap.get_partition_handler(),
											pMG->begin<VertexBase>(), pMG->end<VertexBase>(),
											numCellsX, numCellsY, aaPos,
											Grid::vertex_traits::cb_consider_all, bucketSubset);
			else
				PartitionElements_RegularGrid<VertexBase>(
											partitionMap.get_partition_handler(),
											pMG->begin<VertexBase>(), pMG->end<VertexBase>(),
											numCellsX, numCellsY, aaPos,
											cbConsiderElem, bucketSubset);
		}
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

		return true;
	#endif

	UG_LOG("WARNING: PartitionDomain_RegularGrid is currently only implemented for");
	UG_LOG(" parallel environments.\n");
	return false;
}

template <typename TDomain>
static bool
PartitionDomain_MetisKWay(TDomain& domain, PartitionMap& partitionMap,
						  int numPartitions, size_t baseLevel,
						  int hWeight, int vWeight)
{
	PROFILE_FUNC_GROUP("parallelization")
//	prepare the partition map
	SmartPtr<MultiGrid> pMG = domain.grid();
	partitionMap.assign_grid(*pMG);

#ifdef UG_PARALLEL
//	we need a process to which elements which are not considered will be send.
//	Those elements should stay on the current process.
	int localProc = 0;
	localProc = pcl::GetProcRank();

	int bucketSubset = partitionMap.find_target_proc(localProc);
	if(bucketSubset == -1)
		bucketSubset = (int)partitionMap.num_target_procs();

	SubsetHandler& shPart = partitionMap.get_partition_handler();
//	call the actual partitioning routine
	if(pMG->num<Volume>() > 0){
		PartitionMultiGrid_MetisKway<Volume>(shPart, *pMG, numPartitions,
											baseLevel, hWeight, vWeight);
	//	assign all elements below baseLevel to bucketSubset
		for(size_t lvl = 0; lvl < baseLevel; ++lvl)
			shPart.assign_subset(pMG->begin<Volume>(lvl), pMG->end<Volume>(lvl),
								 bucketSubset);
	}
	else if(pMG->num<Face>() > 0){
		PartitionMultiGrid_MetisKway<Face>(shPart, *pMG, numPartitions,
											baseLevel, hWeight, vWeight);
	//	assign all elements below baseLevel to bucketSubset
		for(size_t lvl = 0; lvl < baseLevel; ++lvl)
			shPart.assign_subset(pMG->begin<Face>(lvl), pMG->end<Face>(lvl),
								 bucketSubset);
	}
	else if(pMG->num<EdgeBase>() > 0){
		PartitionMultiGrid_MetisKway<EdgeBase>(shPart, *pMG, numPartitions,
												baseLevel, hWeight, vWeight);
	//	assign all elements below baseLevel to bucketSubset
		for(size_t lvl = 0; lvl < baseLevel; ++lvl)
			shPart.assign_subset(pMG->begin<EdgeBase>(lvl), pMG->end<EdgeBase>(lvl),
								 bucketSubset);
	}

	if(!partitionMap.get_partition_handler().empty(bucketSubset)){
		if(bucketSubset >= (int)partitionMap.num_target_procs())
			partitionMap.add_target_proc(localProc);
	}

	return true;
#else
	UG_LOG("WARNING in PartitionDomain_MetisKWay: Only available in parallel builds.\n");
	return false;
#endif
}


template <typename TDomain>
static bool
PartitionDomain_LevelBased(TDomain& domain, PartitionMap& partitionMap,
						  	   int numPartitions, size_t level)
{
	PROFILE_FUNC_GROUP("parallelization")
	//	prepare the partition map
	SmartPtr<MultiGrid> pMG = domain.grid();
	partitionMap.assign_grid(*pMG);
	SubsetHandler& shPart = partitionMap.get_partition_handler();
//	call the actual partitioning routine
	if(pMG->num<Volume>() > 0)
		PartitionMultiGridLevel_MetisKway<Volume>(shPart, *pMG, numPartitions, level);

	else if(pMG->num<Face>() > 0)
		PartitionMultiGridLevel_MetisKway<Face>(shPart, *pMG, numPartitions, level);

	else if(pMG->num<EdgeBase>() > 0)
		PartitionMultiGridLevel_MetisKway<EdgeBase>(shPart, *pMG, numPartitions, level);

	return true;
}

template <typename TDomain>
static bool RedistributeDomain(TDomain& domainOut,
							   PartitionMap& partitionMap,
							   bool createVerticalInterfaces)
{
	PROFILE_FUNC_GROUP("parallelization")
//todo	Use a process-communicator to restrict communication

	typedef typename TDomain::position_attachment_type	position_attachment_type;
//	make sure that the input is fine
	typedef typename TDomain::grid_type GridType;
	SmartPtr<GridType> pGrid = domainOut.grid();
	SubsetHandler& shPart = partitionMap.get_partition_handler();

	if(shPart.grid() != pGrid.get()){
		partitionMap.assign_grid(*pGrid);
	}

#ifdef UG_PARALLEL
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

	PCL_PROFILE(RedistributeDomain);

//	make sure that manager exists
	DistributedGridManager* pDistGridMgr = domainOut.distributed_grid_manager();
	if(!pDistGridMgr)
	{
		UG_LOG("ERROR in RedistibuteDomain: Domain has to feature a Distributed Grid Manager.\n");
		performDistribution = false;
	}

//todo	Use a process-communicator to restrict communication
	if(!pcl::AllProcsTrue(performDistribution))
		return false;

	DistributedGridManager& distGridMgr = *pDistGridMgr;

//	data serialization
	GeomObjAttachmentSerializer<VertexBase, position_attachment_type>
		posSerializer(*pGrid, domainOut.position_attachment());
	SubsetHandlerSerializer shSerializer(*domainOut.subset_handler());

	GridDataSerializationHandler serializer;
	serializer.add(&posSerializer);
	serializer.add(&shSerializer);

//	now call redistribution
	RedistributeGrid(distGridMgr, shPart, serializer, serializer,
					 createVerticalInterfaces, &partitionMap.get_target_proc_vec());

	PCL_PROFILE_END();
#endif

//	in the serial case there's nothing to do.
	return true;
}



template <typename TDomain>
static bool DistributeDomain(TDomain& domainOut)
{
	PROFILE_FUNC_GROUP("parallelization")
#ifdef UG_PARALLEL
//	typedefs
	typedef typename TDomain::subset_handler_type subset_handler_type;

//	get distributed grid manager
	DistributedGridManager* pDistGridMgr = domainOut.distributed_grid_manager();

//	check that manager exists
	if(!pDistGridMgr)
	{
		UG_LOG("DistibuteDomain: Cannot find Distributed Grid Manager.\n");
		return false;
	}
	DistributedGridManager& distGridMgrOut = *pDistGridMgr;

//	get subset handler
	SmartPtr<subset_handler_type> pSH = domainOut.subset_handler();

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
	position_attachment_type& domPosition = domainOut.position_attachment();
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
	AdjustAndDistributeGrid(distGridMgrOut, *pSH, mg, *pSH, numProcs, true,
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
static bool DistributeDomain(TDomain& domainOut, PartitionMap& partitionMap)
{
	PROFILE_FUNC_GROUP("parallelization")
#ifdef UG_PARALLEL
//	typedefs
	typedef typename TDomain::subset_handler_type subset_handler_type;

//	get distributed grid manager
	DistributedGridManager* pDistGridMgr = domainOut.distributed_grid_manager();

//	check that manager exists
	if(!pDistGridMgr)
	{
		UG_LOG("DistibuteDomain: Cannot find Distributed Grid Manager.\n");
		return false;
	}
	DistributedGridManager& distGridMgrOut = *pDistGridMgr;

//	get subset handler
	SmartPtr<subset_handler_type> pSH = domainOut.subset_handler();

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
	GridLayoutMap& glmOut = distGridMgrOut.grid_layout_map();

//	make sure that each grid has a position attachment - even if no data
//	will be received.
	typedef typename TDomain::position_attachment_type position_attachment_type;
	position_attachment_type& domPosition = domainOut.position_attachment();
	bool tmpPosAttachment = false;
	if(!mg.has_vertex_attachment(aPosition))
	{
	// convert to 3d positions (FVGeometry depends on PositionCoordinates)
       mg.attach_to_vertices(aPosition);
       ConvertMathVectorAttachmentValues<VertexBase>(mg, domPosition, aPosition);

       UG_LOG("DistributeDomain: temporarily adding Position Attachment.\n");
       tmpPosAttachment = true;
	}

	if(pcl::GetProcRank() == 0){
		DistributeGrid_KeepSrcGrid(mg, *pSH, glmOut,
								   partitionMap.get_partition_handler(),
								   0, false);
	}
	else{
		UG_LOG("receiving grid...\n");
			ReceiveGrid(mg, *pSH, glmOut, 0, true);
		UG_LOG("receive done.\n");
	}

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

}//	end of namespace

#endif
