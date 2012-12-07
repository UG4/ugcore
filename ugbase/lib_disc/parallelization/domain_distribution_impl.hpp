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
	switch(domain.domain_info().element_type()){
		case VOLUME:
			PartitionMultiGridLevel_MetisKway<Volume>(shPart, *pMG, numPartitions, level);
			break;

		case FACE:
			PartitionMultiGridLevel_MetisKway<Face>(shPart, *pMG, numPartitions, level);
			break;

		case EDGE:
			PartitionMultiGridLevel_MetisKway<EdgeBase>(shPart, *pMG, numPartitions, level);
			break;

		default:
			UG_THROW("Partitioning only works for element types EDGE, FACE, and VOLUME!");
			break;
	}

	return true;
}


template <typename TDomain>
static bool
PartitionDistributedDomain_LevelBased(TDomain& domain, PartitionMap& partitionMap,
						  	   	   	   	   int numPartitions, size_t level)
{
	PROFILE_FUNC_GROUP("parallelization")
	//	prepare the partition map
	SmartPtr<MultiGrid> pMG = domain.grid();
	partitionMap.assign_grid(*pMG);
	SubsetHandler& shPart = partitionMap.get_partition_handler();

//	call the actual partitioning routine
	switch(domain.domain_info().element_type()){
		case VOLUME:
			PartitionMultiGridLevel_ParmetisKway<Volume>(shPart, *pMG, numPartitions, level);
			break;

		case FACE:
			PartitionMultiGridLevel_ParmetisKway<Face>(shPart, *pMG, numPartitions, level);
			break;

		case EDGE:
			PartitionMultiGridLevel_ParmetisKway<EdgeBase>(shPart, *pMG, numPartitions, level);
			break;

		default:
			UG_THROW("Partitioning only works for element types EDGE, FACE, and VOLUME!");
			break;
	}

	return true;
}


template <typename TDomain>
static bool DistributeDomain(TDomain& domainOut,
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
	//bool performDistribution = true;

//	make sure that the number of subsets and target processes match
//	THIS MAKES NO SENSE FOR PARALLEL REDISTRIBUTION - IT IS CLEAR THAT SOME
//	PROCS WON'T DELIVER TO ALL PROCS IN THE MAP.
/*	const int numSubs = partitionMap.get_partition_handler().num_subsets();
	const int numTargetProcs = (int)partitionMap.num_target_procs();
	if(numSubs > numTargetProcs){
		UG_LOG("ERROR in RedistributeDomain: More partitions than target processes.\n");
		performDistribution = false;
	}
	else if(numSubs < numTargetProcs){
		UG_LOG("ERROR in RedistributeDomain: More target processes than partitions.\n");
		performDistribution = false;
	}
*/

//todo:	check whether all target-processes in partitionMap are in the valid range.

	PCL_PROFILE(RedistributeDomain);

//todo	Use a process-communicator to restrict communication
/*
	if(!pcl::AllProcsTrue(performDistribution))
		return false;
*/

//	data serialization
	GeomObjAttachmentSerializer<VertexBase, position_attachment_type>
		posSerializer(*pGrid, domainOut.position_attachment());
	SubsetHandlerSerializer shSerializer(*domainOut.subset_handler());

	GridDataSerializationHandler serializer;
	serializer.add(&posSerializer);
	serializer.add(&shSerializer);

//	now call redistribution
	DistributeGrid(*pGrid, shPart, serializer, serializer,
					 createVerticalInterfaces, &partitionMap.get_target_proc_vec());

	PCL_PROFILE_END();
#endif

//	in the serial case there's nothing to do.
	return true;
}

}//	end of namespace

#endif
