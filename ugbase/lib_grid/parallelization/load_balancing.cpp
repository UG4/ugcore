// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 09.06.2011 (m,d,y)
 
#include "load_balancing.h"
#include "lib_grid/algorithms/graph/dual_graph.h"

//	if we're using metis, then include the header now
#ifdef UG_PARALLEL
#include "pcl/pcl_base.h"
#include "pcl/pcl_interface_communicator.h"
#include "util/parallel_dual_graph.h"
#include "lib_grid/parallelization/util/compol_subset.h"
#endif

#ifdef UG_METIS
	extern "C" {
		#include "metis.h"
	}
#endif

#ifdef UG_PARMETIS
	extern "C" {
		#include "parmetis.h"
	}
#endif


using namespace std;


namespace ug{

////////////////////////////////////////////////////////////////////////////////
template <class TGeomBaseObj>
bool PartitionGrid_MetisKway(SubsetHandler& shPartitionOut,
							 Grid& grid, int numParts)
{
	#ifdef UG_METIS
		if(numParts == 1){
			shPartitionOut.assign_subset(grid.begin<TGeomBaseObj>(),
										grid.end<TGeomBaseObj>(), 0);
			return true;
		}

	//	construct the dual graph to the grid
		vector<idx_t> adjacencyMapStructure;
		vector<idx_t> adjacencyMap;
		ConstructDualGraph<TGeomBaseObj, idx_t>(adjacencyMapStructure,
												  adjacencyMap, grid);

	//	partition the graph using metis
	//	first define options for the partitioning.
		idx_t options[METIS_NOPTIONS];
		METIS_SetDefaultOptions(options);
		options[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY;
		options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
		options[METIS_OPTION_NUMBERING] = 0;
	  //request contiguous partitions
		//options[METIS_OPTION_CONTIG] = 1;
	  //note: using the option METIS_OPTION_DBGLVL could be useful for debugging.

		int nVrts = (int)adjacencyMapStructure.size() - 1;
		int nConstraints = 1;
		int edgeCut;
		vector<idx_t> partitionMap(nVrts);

		UG_DLOG(LIB_GRID, 1, "CALLING METIS\n");
		int metisRet =	METIS_PartGraphKway(&nVrts, &nConstraints,
											&adjacencyMapStructure.front(),
											&adjacencyMap.front(),
											NULL, NULL, NULL,
											&numParts, NULL, NULL, options,
											&edgeCut, &partitionMap.front());
		UG_DLOG(LIB_GRID, 1, "METIS DONE\n");

		if(metisRet != METIS_OK){
			UG_DLOG(LIB_GRID, 1, "METIS FAILED\n");
			return false;
		}

	//	assign the subsets to the subset-handler
		int counter = 0;
		typedef typename geometry_traits<TGeomBaseObj>::iterator iterator;
		for(iterator iter = grid.begin<TGeomBaseObj>();
			iter != grid.end<TGeomBaseObj>(); ++iter)
		{
			shPartitionOut.assign_subset(*iter, partitionMap[counter++]);
		}

		return true;
	#else
		UG_LOG("WARNING in PartitionGrid_MetisKway: METIS is not available in "
				"the current build. Please consider recompiling with METIS "
				"support enabled.\n");
		return false;
	#endif
}

//////////////////////////////
//	explicit instantiation
template bool PartitionGrid_MetisKway<EdgeBase>(SubsetHandler&, Grid&, int);
template bool PartitionGrid_MetisKway<Face>(SubsetHandler&, Grid&, int);
template bool PartitionGrid_MetisKway<Volume>(SubsetHandler&, Grid&, int);


////////////////////////////////////////////////////////////////////////////////
template <class TGeomBaseObj>
bool PartitionMultiGrid_MetisKway(SubsetHandler& shPartitionOut,
							 	  MultiGrid& mg, int numParts,
							 	  size_t baseLevel,
							 	  int hWeight, int vWeight)
{
#ifdef UG_METIS
	typedef TGeomBaseObj	TElem;
	typedef typename geometry_traits<TGeomBaseObj>::iterator	ElemIter;

//	only call metis if more than 1 part is required
	if(numParts > 1){
	//	here we'll store the dual graph
		vector<idx_t> adjacencyMapStructure;
		vector<idx_t> adjacencyMap;
		vector<idx_t> edgeWeightMap;

	//todo	add baseLevel to ConstructDualGraphMG.
		ConstructDualGraphMG<TGeomBaseObj, idx_t>(adjacencyMapStructure,
													adjacencyMap, &edgeWeightMap,
													mg, baseLevel, hWeight, vWeight);

	//	partition the graph using metis
	//	first define options for the partitioning.
		idx_t options[METIS_NOPTIONS];
		METIS_SetDefaultOptions(options);
		options[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY;
		options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
		options[METIS_OPTION_NUMBERING] = 0;
	  //request contiguous partitions
		//options[METIS_OPTION_CONTIG] = 1;
	  //note: using the option METIS_OPTION_DBGLVL could be useful for debugging.

		int nVrts = (int)adjacencyMapStructure.size() - 1;
		int nConstraints = 1;
		int edgeCut;
		vector<idx_t> partitionMap(nVrts);

		UG_DLOG(LIB_GRID, 1, "CALLING METIS\n");
		int metisRet =	METIS_PartGraphKway(&nVrts, &nConstraints,
											&adjacencyMapStructure.front(),
											&adjacencyMap.front(),
											NULL, NULL, &edgeWeightMap.front(),
											&numParts, NULL, NULL, options,
											&edgeCut, &partitionMap.front());
		UG_DLOG(LIB_GRID, 1, "METIS DONE\n");

		if(metisRet != METIS_OK){
			UG_DLOG(LIB_GRID, 1, "METIS FAILED\n");
			return false;
		}


	//	assign the subsets to the subset-handler
		int counter = 0;
		for(size_t lvl = baseLevel; lvl < mg.num_levels(); ++lvl){
			typedef typename geometry_traits<TGeomBaseObj>::iterator iterator;
			for(iterator iter = mg.begin<TGeomBaseObj>(lvl);
				iter != mg.end<TGeomBaseObj>(lvl); ++iter)
			{
				shPartitionOut.assign_subset(*iter, partitionMap[counter++]);
			}
		}
	}
	else{
	//	assign all elements to subset 0.
		for(size_t lvl = baseLevel; lvl < mg.num_levels(); ++lvl){
			shPartitionOut.assign_subset(mg.begin<TGeomBaseObj>(lvl),
										 mg.end<TGeomBaseObj>(lvl), 0);
		}
	}
	return true;
#else
	UG_LOG("WARNING in PartitionMultiGrid_MetisKway: METIS is not available in "
			"the current build. Please consider recompiling with METIS "
			"support enabled.\n");
	return false;
#endif
}

////////////////////////////////////////////////////////////////////////////////
template <class TGeomBaseObj>
bool PartitionMultiGrid_MetisKway(SubsetHandler& shPartitionOut,
							 	  MultiGrid& mg, int numParts, size_t baseLevel,
							 	  boost::function<int (TGeomBaseObj*, TGeomBaseObj*)>& weightFct)
{
#ifdef UG_METIS
	typedef TGeomBaseObj	TElem;
	typedef typename geometry_traits<TGeomBaseObj>::iterator	ElemIter;
//	only call metis if more than 1 part is required
	//if(numParts > 1){
	//	here we'll store the dual graph
		vector<idx_t> adjacencyMapStructure;
		vector<idx_t> adjacencyMap;
		vector<idx_t> edgeWeightMap;

		if (weightFct)
		{
			// find out how many elems there are
			size_t size = 0;
			for (size_t lvl = baseLevel; lvl < mg.num_levels(); lvl++)
				size += mg.num<TElem>(lvl);

			// prepare list of all elements
			vector<TElem*> elems(size);

			// construct dual graph with standard edge weights
			ConstructDualGraphMG<TGeomBaseObj, idx_t>(adjacencyMapStructure,
														adjacencyMap, &edgeWeightMap,
														mg, baseLevel, 1, 1,
														NULL, &elems.front());

			// correct edge weights by calling weightFct for all pairs of TElems
			// found by ConstructDualGraphMG
			UG_ASSERT(size+1 == adjacencyMapStructure.size(), "index lengths not matching");
			TElem* elem1;
			TElem* elem2;
			for (size_t el = 0; el < size; el++)
			{
				elem1 = elems[el];
				for (size_t edgeInd = (size_t) adjacencyMapStructure[el]; edgeInd < (size_t) adjacencyMapStructure[el+1]; edgeInd++)
				{
					elem2 = elems[(size_t) adjacencyMap[edgeInd]];
					edgeWeightMap[edgeInd] = weightFct(elem1, elem2);
				}
			}
		}
		else
		{
		//todo	add baseLevel to ConstructDualGraphMG.
			ConstructDualGraphMG<TGeomBaseObj, idx_t>(adjacencyMapStructure,
													  adjacencyMap, &edgeWeightMap,
													  mg, baseLevel, 1, 1);
		}
	if(numParts > 1){
	//	partition the graph using metis
	//	first define options for the partitioning.
		idx_t options[METIS_NOPTIONS];
		METIS_SetDefaultOptions(options);
		options[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY;
		options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
		options[METIS_OPTION_NUMBERING] = 0;
	  //request contiguous partitions
		//options[METIS_OPTION_CONTIG] = 1;
	  //note: using the option METIS_OPTION_DBGLVL could be useful for debugging.

		int nVrts = (int)adjacencyMapStructure.size() - 1;
		int nConstraints = 1;
		int edgeCut;
		vector<idx_t> partitionMap(nVrts);

		UG_DLOG(LIB_GRID, 1, "CALLING METIS\n");
		int metisRet =	METIS_PartGraphKway(&nVrts, &nConstraints,
											&adjacencyMapStructure.front(),
											&adjacencyMap.front(),
											NULL, NULL, &edgeWeightMap.front(),
											&numParts, NULL, NULL, options,
											&edgeCut, &partitionMap.front());
		UG_DLOG(LIB_GRID, 1, "METIS DONE\n");

		if(metisRet != METIS_OK){
			UG_DLOG(LIB_GRID, 1, "METIS FAILED\n");
			return false;
		}


	//	assign the subsets to the subset-handler
		int counter = 0;
		for(size_t lvl = baseLevel; lvl < mg.num_levels(); ++lvl){
			typedef typename geometry_traits<TGeomBaseObj>::iterator iterator;
			for(iterator iter = mg.begin<TGeomBaseObj>(lvl);
				iter != mg.end<TGeomBaseObj>(lvl); ++iter)
			{
				shPartitionOut.assign_subset(*iter, partitionMap[counter++]);
			}
		}
	}
	else{
	//	assign all elements to subset 0.
		for(size_t lvl = baseLevel; lvl < mg.num_levels(); ++lvl){
			shPartitionOut.assign_subset(mg.begin<TGeomBaseObj>(lvl),
										 mg.end<TGeomBaseObj>(lvl), 0);
		}
	}
	return true;
#else
	UG_LOG("WARNING in PartitionMultiGrid_MetisKway: METIS is not available in "
			"the current build. Please consider recompiling with METIS "
			"support enabled.\n");
	return false;
#endif
}

////////////////////////////////////////////////////////////////////////////////
template <class TGeomBaseObj>
bool PartitionMultiGridLevel_MetisKway(SubsetHandler& shPartitionOut,
							 	  MultiGrid& mg, int numParts, size_t level)
{
	if(level > mg.top_level()){
		UG_LOG("  WARNING in PartitionMultiGridLevel_MetisKway:"
				"Specified level too high: toplevel=" << mg.top_level()
				<< ", specified-level=" << level <<"\n");
		return true;
	}

	if(mg.num<TGeomBaseObj>() == 0)
		return true;

#ifdef UG_METIS
	typedef TGeomBaseObj	TElem;
	typedef typename geometry_traits<TGeomBaseObj>::iterator	ElemIter;

//	only call metis if more than 1 part is required
	int rootProc = 0;
	#ifdef UG_PARALLEL
		rootProc = pcl::GetProcRank();
	#endif

	if(numParts > 1){
	//	here we'll store the dual graph
		vector<idx_t> adjacencyMapStructure;
		vector<idx_t> adjacencyMap;


	//	we can optionally use a higher edge weight on siblings.
	//	note - higher memory and processing requirements if enabled.
		const bool useEdgeWeights = true;
		const idx_t siblingWeight = 2;

	//	we'll reuse the index later on during assignment of the edge weights
		typedef Attachment<idx_t> AIndex;
		AIndex aIndex;
		vector<TElem*> elems;

		if(useEdgeWeights){
			mg.attach_to<TElem>(aIndex);
			Grid::AttachmentAccessor<TElem, AIndex> aaInd;
			elems.resize(mg.num<TElem>(level));

		//	we also need a vector of elements, so that we can access them by index
		//	(also for the caclulation of edge weights)
			ConstructDualGraphMGLevel<TElem, idx_t>(adjacencyMapStructure, adjacencyMap,
													mg, level, &aIndex, &elems.front());
		}
		else{
			ConstructDualGraphMGLevel<TElem, idx_t>(adjacencyMapStructure, adjacencyMap,
												mg, level);
		}

	//	partition the graph using metis
	//	first define options for the partitioning.
		idx_t options[METIS_NOPTIONS];
		METIS_SetDefaultOptions(options);
		options[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY;
		options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
		options[METIS_OPTION_NUMBERING] = 0;
	  //request contiguous partitions
		//options[METIS_OPTION_CONTIG] = 1;
	  //note: using the option METIS_OPTION_DBGLVL could be useful for debugging.

		idx_t nVrts = (idx_t)adjacencyMapStructure.size() - 1;
		idx_t nConstraints = 1;
		idx_t edgeCut;
		vector<idx_t> partitionMap(nVrts);

	//	create a weight map for the vertices based on the number of children+1
	//	for each graph-vertex. This is not necessary, if we're already on the top level
		idx_t* pVrtSizeMap = NULL;
		vector<idx_t> vrtSizeMap;
		if(level < mg.top_level()){
			vrtSizeMap.reserve(nVrts);
			for(ElemIter iter = mg.begin<TElem>(level);
				iter != mg.end<TElem>(level); ++iter)
			{
				vrtSizeMap.push_back((mg.num_children_total(*iter) + 1));
			}
			assert((int)vrtSizeMap.size() == nVrts);
			pVrtSizeMap = &vrtSizeMap.front();
		}

	//	we'll also create weights for the edges, since we want to cluster elements,
	//	which have the same parent (this reduces several problems later on, like
	//	additional vertical interfaces)
		idx_t* pEdgeWeights = NULL;
		vector<idx_t> edgeWeights;
		if(useEdgeWeights){
			edgeWeights.reserve(adjacencyMap.size());
			for(int i_vrt = 0; i_vrt < nVrts; ++i_vrt){
				int start = adjacencyMapStructure[i_vrt];
				int end = adjacencyMapStructure[i_vrt + 1];

				GeometricObject* parent = mg.get_parent(elems[i_vrt]);
				for(int i_edge = start; i_edge < end; ++i_edge){
					if(parent == mg.get_parent(elems[adjacencyMap[i_edge]]))
						edgeWeights.push_back(siblingWeight);
					else
						edgeWeights.push_back(1);
				}
			}
			pEdgeWeights = &edgeWeights.front();
		}

		UG_DLOG(LIB_GRID, 1, "CALLING METIS\n");
		int metisRet =	METIS_PartGraphKway(&nVrts, &nConstraints,
											&adjacencyMapStructure.front(),
											&adjacencyMap.front(),
											NULL, pVrtSizeMap, pEdgeWeights,
											&numParts, NULL, NULL, options,
											&edgeCut, &partitionMap.front());
		UG_DLOG(LIB_GRID, 1, "METIS DONE\n");

		if(metisRet != METIS_OK){
			UG_DLOG(LIB_GRID, 1, "METIS FAILED\n");
			return false;
		}


	//	assign the subsets to the subset-handler
	//	all subsets below the specified level go to the rootProc.
	//	The ones on the specified level go as METIS assigned them.
	//	All children in levels above copy their from their parent.
		for(size_t lvl = 0; lvl < level; ++lvl)
			for(ElemIter iter = mg.begin<TElem>(lvl); iter != mg.end<TElem>(lvl); ++iter)
				shPartitionOut.assign_subset(*iter, rootProc);

		int counter = 0;
		for(ElemIter iter = mg.begin<TElem>(level); iter != mg.end<TElem>(level); ++iter)
			shPartitionOut.assign_subset(*iter, partitionMap[counter++]);

	//	currently there is a restriction in the functionality of the surface view.
	//	Because of that we have to make sure, that all siblings in the specified level
	//	are sent to the same process... we thus adjust the partition slightly.
	//todo:	Not all siblings should have to be sent to the same process...
	//		simply remove the following code block - make sure that surface-view supports this!
	//		However, problems with discretizations and solvers would occur
		if(level > 0){
		//	put all children in the subset of the first one.
			for(ElemIter iter = mg.begin<TElem>(level-1);
				iter != mg.end<TElem>(level-1); ++iter)
			{
				TElem* e = *iter;
				size_t numChildren = mg.num_children<TElem>(e);
				if(numChildren > 1){
					int partition = shPartitionOut.get_subset_index(mg.get_child<TElem>(e, 0));
					for(size_t i = 1; i < numChildren; ++i)
						shPartitionOut.assign_subset(mg.get_child<TElem>(e, i), partition);
				}
			}
		}

		for(size_t lvl = level; lvl < mg.top_level(); ++lvl){
			for(ElemIter iter = mg.begin<TElem>(lvl); iter != mg.end<TElem>(lvl); ++iter)
			{
				size_t numChildren = mg.num_children<TElem>(*iter);
				for(size_t i = 0; i < numChildren; ++i){
					shPartitionOut.assign_subset(mg.get_child<TElem>(*iter, i),
												shPartitionOut.get_subset_index(*iter));
				}
			}
		}

		if(useEdgeWeights){
			mg.detach_from<TElem>(aIndex);
		}
	}
	else{
	//	assign all elements to subset 0.
		for(size_t lvl = 0; lvl < mg.num_levels(); ++lvl){
			shPartitionOut.assign_subset(mg.begin<TGeomBaseObj>(lvl),
										 mg.end<TGeomBaseObj>(lvl), rootProc);
		}
	}
	return true;
#else
	UG_LOG("WARNING in PartitionMultiGrid_MetisKway: METIS is not available in "
			"the current build. Please consider recompiling with METIS "
			"support enabled.\n");
	return false;
#endif
}

////////////////////////////////////////////////////////////////////////////////
template <class TGeomBaseObj>
bool PartitionMultiGridLevel_ParmetisKway(SubsetHandler& shPartitionOut,
							 	  	  MultiGrid& mg, int numParts, size_t level)
{
	UG_DLOG(LIB_GRID, 1, "start - PartitionMultiGridLevel_ParmetisKway\n");

#if defined UG_PARMETIS && defined UG_PARALLEL
	typedef TGeomBaseObj	TElem;
	typedef typename geometry_traits<TGeomBaseObj>::iterator	ElemIter;

	int localProc = pcl::GetProcRank();

	pcl::ProcessCommunicator procCom;

//	here we'll store the dual graph
	vector<idx_t> adjacencyMapStructure;
	vector<idx_t> adjacencyMap;
	vector<idx_t> nodeOffsetMap;

	ConstructParallelDualGraphMGLevel<TElem, idx_t>(adjacencyMapStructure,
											adjacencyMap, nodeOffsetMap,
											mg, level, procCom);

	UG_DLOG(LIB_GRID, 2, "  parallel dual graph #vrts: " << (int)adjacencyMapStructure.size() - 1
						<< ", #edges: " << adjacencyMap.size() / 2 << "\n");

	if(!pcl::AllProcsTrue(adjacencyMap.size() > 1)){
		UG_THROW("ParMetis may only be executed if all processes contain non-ghost elements"
				" on the given level (at least two neighboring).")
	}

//	partition the graph using parmetis
	idx_t options[3]; options[0] = 0;//default values
	idx_t nVrts = (idx_t)adjacencyMapStructure.size() - 1;
	idx_t nConstraints = 1;
	idx_t edgeCut;
	idx_t wgtFlag = 2;//only vertices are weighted
	idx_t numFlag = 0;
	vector<idx_t> partitionMap(nVrts);
	vector<real_t> tpwgts(numParts, 1. / (number)numParts);
	real_t ubvec = 1.05;
//	create a weight map for the vertices based on the number of children+1
//	for each graph-vertex. This is not necessary, if we're already on the top level
	idx_t* pVrtSizeMap = NULL;
	vector<idx_t> vrtSizeMap;
	{
		vrtSizeMap.reserve(nVrts);
		for(ElemIter iter = mg.begin<TElem>(level);
			iter != mg.end<TElem>(level); ++iter)
		{
			vrtSizeMap.push_back((mg.num_children_total(*iter) + 1));
		}
		assert((int)vrtSizeMap.size() == nVrts);
		pVrtSizeMap = &vrtSizeMap.front();
	}

	UG_DLOG(LIB_GRID, 1, "CALLING PARMETIS\n");
	MPI_Comm mpiCom = procCom.get_mpi_communicator();
	int metisRet =	ParMETIS_V3_PartKway(&nodeOffsetMap.front(),
										&adjacencyMapStructure.front(),
										&adjacencyMap.front(),
										pVrtSizeMap, NULL, &wgtFlag,
										&numFlag, &nConstraints,
										&numParts, &tpwgts.front(), &ubvec, options,
										&edgeCut, &partitionMap.front(),
										&mpiCom);
	UG_DLOG(LIB_GRID, 1, "PARMETIS DONE\n");

	if(metisRet != METIS_OK){
		UG_THROW("PARMETIS FAILED on process " << localProc);
	}

//	assign the subsets to the subset-handler
//	all subsets below the specified level go to the rootProc.
//	The ones on the specified level go as METIS assigned them.
//	All children in levels above copy their from their parent.
	for(size_t lvl = 0; lvl < level; ++lvl)
		for(ElemIter iter = mg.begin<TElem>(lvl); iter != mg.end<TElem>(lvl); ++iter)
			shPartitionOut.assign_subset(*iter, localProc);

	int counter = 0;
	for(ElemIter iter = mg.begin<TElem>(level); iter != mg.end<TElem>(level); ++iter)
		shPartitionOut.assign_subset(*iter, partitionMap[counter++]);



	typedef typename GridLayoutMap::Types<TElem>::Layout::LevelLayout	ElemLayout;
	GridLayoutMap& glm = mg.distributed_grid_manager()->grid_layout_map();
	pcl::InterfaceCommunicator<ElemLayout>	com;
	ComPol_Subset<ElemLayout>	compolSHCopy(shPartitionOut, true);

//	copy subset indices from vertical slaves to vertical masters,
//	since partitioning was only performed on vslaves
	if(glm.has_layout<TElem>(INT_V_SLAVE))
		com.send_data(glm.get_layout<TElem>(INT_V_SLAVE).layout_on_level(level),
					  compolSHCopy);
	if(glm.has_layout<TElem>(INT_V_MASTER))
		com.receive_data(glm.get_layout<TElem>(INT_V_MASTER).layout_on_level(level),
						 compolSHCopy);
	com.communicate();

//todo: make sure that there are no problems at vertical interfaces
//	currently there is a restriction in the functionality of the surface view.
//	Because of that we have to make sure, that all siblings in the specified level
//	are sent to the same process... we thus adjust the partition slightly.
//todo:	Not all siblings should have to be sent to the same process...
//		simply remove the following code block - make sure that surface-view supports this!
//		However, problems with discretizations and solvers would occur.
//		If you remove the following block, consider reducing v-communication.
	if(level > 0){
	//	put all children in the subset of the first one.
		for(ElemIter iter = mg.begin<TElem>(level-1);
			iter != mg.end<TElem>(level-1); ++iter)
		{
			TElem* e = *iter;
			size_t numChildren = mg.num_children<TElem>(e);
			if(numChildren > 1){
				int partition = shPartitionOut.get_subset_index(mg.get_child<TElem>(e, 0));
				for(size_t i = 1; i < numChildren; ++i)
					shPartitionOut.assign_subset(mg.get_child<TElem>(e, i), partition);
			}
		}
	}


	for(size_t lvl = level; lvl < mg.top_level(); ++lvl){
	//	copy subset indices from vertical masters to vertical slaves
		if(glm.has_layout<TElem>(INT_V_MASTER))
			com.send_data(glm.get_layout<TElem>(INT_V_MASTER).layout_on_level(lvl),
						  compolSHCopy);
		if(glm.has_layout<TElem>(INT_V_SLAVE))
			com.receive_data(glm.get_layout<TElem>(INT_V_SLAVE).layout_on_level(lvl),
							 compolSHCopy);
		com.communicate();

		for(ElemIter iter = mg.begin<TElem>(lvl); iter != mg.end<TElem>(lvl); ++iter)
		{
			size_t numChildren = mg.num_children<TElem>(*iter);
			for(size_t i = 0; i < numChildren; ++i){
				shPartitionOut.assign_subset(mg.get_child<TElem>(*iter, i),
											shPartitionOut.get_subset_index(*iter));
			}
		}
	}

//	and a final copy on the top level...
	if(mg.num_levels() > 1){
		if(glm.has_layout<TElem>(INT_V_MASTER))
			com.send_data(glm.get_layout<TElem>(INT_V_MASTER).layout_on_level(mg.top_level()),
						  compolSHCopy);
		if(glm.has_layout<TElem>(INT_V_SLAVE))
			com.receive_data(glm.get_layout<TElem>(INT_V_SLAVE).layout_on_level(mg.top_level()),
							 compolSHCopy);
		com.communicate();
	}

	UG_DLOG(LIB_GRID, 1, "stop - PartitionMultiGridLevel_ParmetisKway\n");
	return true;
#else
	UG_LOG("WARNING in PartitionMultiGrid_MetisKway: METIS is not available in "
			"the current build. Please consider recompiling with METIS "
			"support enabled.\n");
	return false;
#endif
}


//////////////////////////////
//	explicit instantiation
template bool PartitionMultiGrid_MetisKway<EdgeBase>(SubsetHandler&, MultiGrid&,
													 int, size_t, int, int);
template bool PartitionMultiGrid_MetisKway<Face>(SubsetHandler&, MultiGrid&,
												 int, size_t, int, int);
template bool PartitionMultiGrid_MetisKway<Volume>(SubsetHandler&, MultiGrid&,
												   int, size_t, int, int);

template bool PartitionMultiGrid_MetisKway<EdgeBase>(SubsetHandler&, MultiGrid&, int, size_t,
													 boost::function<int (EdgeBase*, EdgeBase*)>&);
template bool PartitionMultiGrid_MetisKway<Face>(SubsetHandler&, MultiGrid&, int, size_t,
		 	 	 	 	 	 	 	 	 	 	 boost::function<int (Face*, Face*)>&);
template bool PartitionMultiGrid_MetisKway<Volume>(SubsetHandler&, MultiGrid&, int, size_t,
												   boost::function<int (Volume*, Volume*)>&);

template bool PartitionMultiGridLevel_MetisKway<EdgeBase>(SubsetHandler&,
													MultiGrid&, int, size_t);
template bool PartitionMultiGridLevel_MetisKway<Face>(SubsetHandler&,
													MultiGrid&, int, size_t);
template bool PartitionMultiGridLevel_MetisKway<Volume>(SubsetHandler&,
													MultiGrid&, int, size_t);

template bool PartitionMultiGridLevel_ParmetisKway<EdgeBase>(SubsetHandler&,
													MultiGrid&, int, size_t);
template bool PartitionMultiGridLevel_ParmetisKway<Face>(SubsetHandler&,
													MultiGrid&, int, size_t);
template bool PartitionMultiGridLevel_ParmetisKway<Volume>(SubsetHandler&,
													MultiGrid&, int, size_t);

}// end of namespace
