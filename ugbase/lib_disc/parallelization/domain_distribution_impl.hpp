// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 25.05.2011 (m,d,y)

#ifndef __H__UG__domain_distribution_impl__
#define __H__UG__domain_distribution_impl__

#include "domain_distribution.h"
#include "lib_grid/algorithms/attachment_util.h"
#include "lib_grid/parallelization/load_balancing.h"
#include "common/serialization.h"

#ifdef UG_PARALLEL
	#include "pcl/pcl.h"
	#include "lib_grid/parallelization/parallelization.h"
#endif


namespace ug
{

///	partitions a domain by sorting all elements into a regular grid
template <typename TDomain>
static bool PartitionDomain_RegularGrid(TDomain& domain, PartitionMap& partitionMap,
										int numCellsX, int numCellsY, int numCellsZ,
										bool surfaceOnly)
{
	PROFILE_FUNC_GROUP("parallelization");
//	prepare the partition map and a vertex position attachment accessor
	SmartPtr<MultiGrid> pMG = domain.grid();
	partitionMap.assign_grid(*pMG);

	#ifdef UG_PARALLEL
	
	SubsetHandler& partitionHandler = *partitionMap.get_partition_handler();
	
	//	a distributed grid manager is required
		if(!domain.distributed_grid_manager()){
			UG_LOG("A distributed grid manager is required in the given domain.\n");
			return false;
		}

		typedef typename TDomain::position_attachment_type TAPos;
		Grid::AttachmentAccessor<Vertex, TAPos> aaPos(*pMG,
												domain.position_attachment());

	//	this callback allows us to only distribute surface elements, which are no ghosts
		IsRegularSurfaceElem cbConsiderElem(*domain.distributed_grid_manager());

	//	we need a process to which elements which are not considered will be send.
	//	Those elements should stay on the current process.
		int localProc = 0;
		localProc = pcl::ProcRank();

		int bucketSubset = partitionMap.find_target_proc(localProc);
		if(bucketSubset == -1)
			bucketSubset = (int)partitionMap.num_target_procs();

	//	partition the grid
		if(pMG->num<Volume>() > 0){
			if(!surfaceOnly)
				PartitionElements_RegularGrid<Volume>(
											partitionHandler,
											pMG->begin<Volume>(), pMG->end<Volume>(),
											numCellsX, numCellsY, numCellsZ, aaPos,
											Grid::volume_traits::cb_consider_all, bucketSubset);
			else
				PartitionElements_RegularGrid<Volume>(
											partitionHandler,
											pMG->begin<Volume>(), pMG->end<Volume>(),
											numCellsX, numCellsY, numCellsZ, aaPos,
											cbConsiderElem, bucketSubset);
		}
		else if(pMG->num<Face>() > 0){
			if(!surfaceOnly)
				PartitionElements_RegularGrid<Face>(
											partitionHandler,
											pMG->begin<Face>(), pMG->end<Face>(),
											numCellsX, numCellsY, numCellsZ, aaPos,
											Grid::face_traits::cb_consider_all, bucketSubset);
			else
				PartitionElements_RegularGrid<Face>(
											partitionHandler,
											pMG->begin<Face>(), pMG->end<Face>(),
											numCellsX, numCellsY, numCellsZ, aaPos,
											cbConsiderElem, bucketSubset);
		}
		else if(pMG->num<Edge>() > 0){
			if(!surfaceOnly)
				PartitionElements_RegularGrid<Edge>(
											partitionHandler,
											pMG->begin<Edge>(), pMG->end<Edge>(),
											numCellsX, numCellsY, numCellsZ, aaPos,
											Grid::edge_traits::cb_consider_all, bucketSubset);
			else
				PartitionElements_RegularGrid<Edge>(
											partitionHandler,
											pMG->begin<Edge>(), pMG->end<Edge>(),
											numCellsX, numCellsY, numCellsZ, aaPos,
											cbConsiderElem, bucketSubset);
		}
		else if(pMG->num<Vertex>() > 0){
			if(!surfaceOnly)
				PartitionElements_RegularGrid<Vertex>(
											partitionHandler,
											pMG->begin<Vertex>(), pMG->end<Vertex>(),
											numCellsX, numCellsY, numCellsZ, aaPos,
											Grid::vertex_traits::cb_consider_all, bucketSubset);
			else
				PartitionElements_RegularGrid<Vertex>(
											partitionHandler,
											pMG->begin<Vertex>(), pMG->end<Vertex>(),
											numCellsX, numCellsY, numCellsZ, aaPos,
											cbConsiderElem, bucketSubset);
		}
		else{
			LOG("partitioning could not be performed - "
				<< "grid doesn't contain any elements!\n");
			return false;
		}

	//	if elements have been assigned to bucketProc, then we have to make sure,
	//	that it is also present in the process-map
		if(!partitionHandler.empty(bucketSubset)){
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
	PROFILE_FUNC_GROUP("parallelization");
//	prepare the partition map
	SmartPtr<MultiGrid> pMG = domain.grid();
	partitionMap.assign_grid(*pMG);

#ifdef UG_PARALLEL
	
	SubsetHandler& partitionHandler = *partitionMap.get_partition_handler();
	
//	we need a process to which elements which are not considered will be send.
//	Those elements should stay on the current process.
	int localProc = 0;
	localProc = pcl::ProcRank();

	int bucketSubset = partitionMap.find_target_proc(localProc);
	if(bucketSubset == -1)
		bucketSubset = (int)partitionMap.num_target_procs();

//	call the actual partitioning routine
	if(pMG->num<Volume>() > 0){
		PartitionMultiGrid_MetisKway<Volume>(partitionHandler, *pMG, numPartitions,
											baseLevel, hWeight, vWeight);
	//	assign all elements below baseLevel to bucketSubset
		for(size_t lvl = 0; lvl < baseLevel; ++lvl)
			partitionHandler.assign_subset(pMG->begin<Volume>(lvl), pMG->end<Volume>(lvl),
								 bucketSubset);
	}
	else if(pMG->num<Face>() > 0){
		PartitionMultiGrid_MetisKway<Face>(partitionHandler, *pMG, numPartitions,
											baseLevel, hWeight, vWeight);
	//	assign all elements below baseLevel to bucketSubset
		for(size_t lvl = 0; lvl < baseLevel; ++lvl)
			partitionHandler.assign_subset(pMG->begin<Face>(lvl), pMG->end<Face>(lvl),
								 bucketSubset);
	}
	else if(pMG->num<Edge>() > 0){
		PartitionMultiGrid_MetisKway<Edge>(partitionHandler, *pMG, numPartitions,
												baseLevel, hWeight, vWeight);
	//	assign all elements below baseLevel to bucketSubset
		for(size_t lvl = 0; lvl < baseLevel; ++lvl)
			partitionHandler.assign_subset(pMG->begin<Edge>(lvl), pMG->end<Edge>(lvl),
								 bucketSubset);
	}

	if(!partitionHandler.empty(bucketSubset)){
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
PartitionDomain_MetisKWay(TDomain& domain, PartitionMap& partitionMap,
						  int numPartitions, size_t baseLevel,
						  SmartPtr<PartitionWeighting> weightFct)
{
	PROFILE_FUNC_GROUP("parallelization");
//	prepare the partition map
	SmartPtr<MultiGrid> pMG = domain.grid();
	partitionMap.assign_grid(*pMG);

#ifdef UG_PARALLEL

	SubsetHandler& partitionHandler = *partitionMap.get_partition_handler();

	PartitionWeighting& wFct = *weightFct;
	wFct.set_subset_handler(domain.subset_handler().operator->());
//	we need a process to which elements which are not considered will be send.
//	Those elements should stay on the current process.
	int localProc = 0;
	localProc = pcl::ProcRank();

	int bucketSubset = partitionMap.find_target_proc(localProc);
	if(bucketSubset == -1)
		bucketSubset = (int)partitionMap.num_target_procs();

//	call the actual partitioning routine
	if(pMG->num<Volume>() > 0){
		// do not use boost::function<...> f = wFct, since this leads to slicing
		// of wFct and losing properties of derived objects
		boost::function<int (Volume*, Volume*)> f = boost::ref(wFct);
		PartitionMultiGrid_MetisKway<Volume>(partitionHandler, *pMG, numPartitions, baseLevel, f);
	//	assign all elements below baseLevel to bucketSubset
		for(size_t lvl = 0; lvl < baseLevel; ++lvl)
			partitionHandler.assign_subset(pMG->begin<Volume>(lvl), pMG->end<Volume>(lvl),
								 bucketSubset);
	}
	else if(pMG->num<Face>() > 0){
		boost::function<int (Face*, Face*)> f = boost::ref(wFct);
		PartitionMultiGrid_MetisKway<Face>(partitionHandler, *pMG, numPartitions, baseLevel, f);
	//	assign all elements below baseLevel to bucketSubset
		for(size_t lvl = 0; lvl < baseLevel; ++lvl)
			partitionHandler.assign_subset(pMG->begin<Face>(lvl), pMG->end<Face>(lvl),
								 bucketSubset);
	}
	else if(pMG->num<Edge>() > 0){
		boost::function<int (Edge*, Edge*)> f = boost::ref(wFct);
		PartitionMultiGrid_MetisKway<Edge>(partitionHandler, *pMG, numPartitions, baseLevel, f);
	//	assign all elements below baseLevel to bucketSubset
		for(size_t lvl = 0; lvl < baseLevel; ++lvl)
			partitionHandler.assign_subset(pMG->begin<Edge>(lvl), pMG->end<Edge>(lvl),
								 bucketSubset);
	}

	if(!partitionHandler.empty(bucketSubset)){
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
	PROFILE_FUNC_GROUP("parallelization");
	//	prepare the partition map
	SmartPtr<MultiGrid> pMG = domain.grid();
	partitionMap.assign_grid(*pMG);
	SubsetHandler& partitionHandler = *partitionMap.get_partition_handler();

//	call the actual partitioning routine
	switch(domain.domain_info().element_type()){
		case VOLUME:
			PartitionMultiGridLevel_MetisKway<Volume>(partitionHandler, *pMG, numPartitions, level);
			break;

		case FACE:
			PartitionMultiGridLevel_MetisKway<Face>(partitionHandler, *pMG, numPartitions, level);
			break;

		case EDGE:
			PartitionMultiGridLevel_MetisKway<Edge>(partitionHandler, *pMG, numPartitions, level);
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
	PROFILE_FUNC_GROUP("parallelization");
	//	prepare the partition map
	SmartPtr<MultiGrid> pMG = domain.grid();
	partitionMap.assign_grid(*pMG);
	SubsetHandler& partitionHandler = *partitionMap.get_partition_handler();

//	call the actual partitioning routine
	switch(domain.domain_info().element_type()){
		case VOLUME:
			PartitionMultiGridLevel_ParmetisKway<Volume>(partitionHandler, *pMG, numPartitions, level);
			break;

		case FACE:
			PartitionMultiGridLevel_ParmetisKway<Face>(partitionHandler, *pMG, numPartitions, level);
			break;

		case EDGE:
			PartitionMultiGridLevel_ParmetisKway<Edge>(partitionHandler, *pMG, numPartitions, level);
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
	PROFILE_FUNC_GROUP("parallelization");
//todo	Use a process-communicator to restrict communication

	typedef typename TDomain::position_attachment_type	position_attachment_type;
//	make sure that the input is fine
	typedef typename TDomain::grid_type GridType;
	SmartPtr<GridType> pGrid = domainOut.grid();
	SubsetHandler& partitionHandler = *partitionMap.get_partition_handler();

	if(partitionHandler.grid() != pGrid.get()){
		partitionMap.assign_grid(*pGrid);
	}

#ifdef UG_PARALLEL
//	used to check whether all processes are correctly prepared for redistribution
	//bool performDistribution = true;

//	make sure that the number of subsets and target processes match
//	THIS MAKES NO SENSE FOR PARALLEL REDISTRIBUTION - IT IS CLEAR THAT SOME
//	PROCS WON'T DELIVER TO ALL PROCS IN THE MAP.
/*	const int numSubs = partitionHandler.num_subsets();
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
	SPVertexDataSerializer posSerializer =
			GeomObjAttachmentSerializer<Vertex, position_attachment_type>::
								create(*pGrid, domainOut.position_attachment());

	SPGridDataSerializer shSerializer = SubsetHandlerSerializer::
											create(*domainOut.subset_handler());

	GridDataSerializationHandler serializer;
	serializer.add(posSerializer);
	serializer.add(shSerializer);

	std::vector<std::string> additionalSHNames = domainOut.additional_subset_handler_names();
	for(size_t i = 0; i < additionalSHNames.size(); ++i){
		SmartPtr<ISubsetHandler> sh = domainOut.additional_subset_handler(additionalSHNames[i]);
		if(sh.valid()){
			SPGridDataSerializer shSerializer = SubsetHandlerSerializer::create(*sh);
			serializer.add(shSerializer);
		}
	}

//	now call redistribution
	DistributeGrid(*pGrid, partitionHandler, serializer, createVerticalInterfaces,
				   &partitionMap.get_target_proc_vec());

	PCL_PROFILE_END();
#endif

//	in the serial case there's nothing to do.
	return true;
}

}//	end of namespace

#endif
