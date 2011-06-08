// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 25.05.2011 (m,d,y)

#ifndef __H__UG__domain_distribution_impl__
#define __H__UG__domain_distribution_impl__

#include "domain_distribution.h"

namespace ug
{

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
										int numCellsX, int numCellsY,
										bool surfaceOnly)
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
		if(mg.num<Volume>() > 0){
			if(!surfaceOnly)
				PartitionElements_RegularGrid<Volume>(
											partitionMap.get_partition_handler(),
											mg.begin<Volume>(), mg.end<Volume>(),
											numCellsX, numCellsY, aaPos,
											(bool(*)(Volume*))ConsiderAll, bucketSubset);
			else
				PartitionElements_RegularGrid<Volume>(
											partitionMap.get_partition_handler(),
											mg.begin<Volume>(), mg.end<Volume>(),
											numCellsX, numCellsY, aaPos,
											cbConsiderElem, bucketSubset);
		}
		else if(mg.num<Face>() > 0){
			if(!surfaceOnly)
				PartitionElements_RegularGrid<Face>(
											partitionMap.get_partition_handler(),
											mg.begin<Face>(), mg.end<Face>(),
											numCellsX, numCellsY, aaPos,
											(bool(*)(Face*))ConsiderAll, bucketSubset);
			else
				PartitionElements_RegularGrid<Face>(
											partitionMap.get_partition_handler(),
											mg.begin<Face>(), mg.end<Face>(),
											numCellsX, numCellsY, aaPos,
											cbConsiderElem, bucketSubset);
		}
		else if(mg.num<EdgeBase>() > 0){
			if(!surfaceOnly)
				PartitionElements_RegularGrid<EdgeBase>(
											partitionMap.get_partition_handler(),
											mg.begin<EdgeBase>(), mg.end<EdgeBase>(),
											numCellsX, numCellsY, aaPos,
											(bool(*)(EdgeBase*))ConsiderAll, bucketSubset);
			else
				PartitionElements_RegularGrid<EdgeBase>(
											partitionMap.get_partition_handler(),
											mg.begin<EdgeBase>(), mg.end<EdgeBase>(),
											numCellsX, numCellsY, aaPos,
											cbConsiderElem, bucketSubset);
		}
		else if(mg.num<VertexBase>() > 0){
			if(!surfaceOnly)
				PartitionElements_RegularGrid<VertexBase>(
											partitionMap.get_partition_handler(),
											mg.begin<VertexBase>(), mg.end<VertexBase>(),
											numCellsX, numCellsY, aaPos,
											(bool(*)(VertexBase*))ConsiderAll, bucketSubset);
			else
				PartitionElements_RegularGrid<VertexBase>(
											partitionMap.get_partition_handler(),
											mg.begin<VertexBase>(), mg.end<VertexBase>(),
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
	#endif

	return true;
}

template <typename TDomain>
static bool
PartitionDomain_MetisKWay(TDomain& domain, PartitionMap& partitionMap,
						  int numPartitions)
{
//	prepare the partition map
	MultiGrid& mg = domain.get_grid();
	partitionMap.assign_grid(mg);

#ifdef UG_PARALLEL
//	call the actual partitioning routine
	if(mg.num<Volume>() > 0){
		PartitionGrid_MetisKway<Volume>(partitionMap.get_partition_handler(),
										mg, numPartitions);
	}
	else if(mg.num<Face>() > 0){
		PartitionGrid_MetisKway<Face>(partitionMap.get_partition_handler(),
										mg, numPartitions);
	}
	else if(mg.num<EdgeBase>() > 0){
		PartitionGrid_MetisKway<EdgeBase>(partitionMap.get_partition_handler(),
										  mg, numPartitions);
	}
	return true;
#else
	UG_LOG("WARNING in PartitionDomain_MetisKWay: Only available in parallel builds.\n");
	return false;
#endif
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

}//	end of namespace

#endif
