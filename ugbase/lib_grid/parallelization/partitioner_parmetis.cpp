// created by Sebastian Reiter
// s.b.reiter@gmail.com
// Feb 25, 2013 (d,m,y)

#include "partitioner_parmetis.h"
#include "load_balancer_util.h"
#include "distributed_grid.h"
#include "lib_grid/parallelization/util/compol_copy_attachment.h"
#include "lib_grid/parallelization/util/compol_subset.h"
#include "lib_grid/algorithms/attachment_util.h"
#include "lib_grid/algorithms/graph/dual_graph.h"
#include "lib_grid/algorithms/geom_obj_util/misc_util.h"
#include "common/util/table.h"

using namespace std;

namespace ug{

static const int CHILD_WEIGHT = 1;
static const int SIBLING_WEIGHT = 100;


template <int dim>
Partitioner_Parmetis<dim>::
Partitioner_Parmetis() :
	m_mg(NULL)
{
	m_processHierarchy = SPProcessHierarchy(new ProcessHierarchy);
	m_balanceWeights = SPBalanceWeights(new StdBalanceWeights<dim>);
	m_connectionWeights = SPConnectionWeights(new StdConnectionWeights<dim>);
}

template<int dim>
Partitioner_Parmetis<dim>::
~Partitioner_Parmetis()
{
}

template<int dim>
void Partitioner_Parmetis<dim>::
set_grid(MultiGrid* mg, Attachment<MathVector<dim> >)
{
	if(mg == m_mg)
		return;

	if(m_mg){
		m_mg->detach_from<elem_t>(m_aNumChildren);
		m_aaNumChildren.invalidate();
		m_sh.assign_grid(NULL);
		m_mg = NULL;
	}

	if(mg){
		m_mg = mg;
		m_sh.assign_grid(m_mg);
		m_mg->attach_to<elem_t>(m_aNumChildren);
		m_aaNumChildren.access(*m_mg, m_aNumChildren);
	}
}

template<int dim>
void Partitioner_Parmetis<dim>::
set_process_hierarchy(SPProcessHierarchy procHierarchy)
{
	m_processHierarchy = procHierarchy;
}

template<int dim>
void Partitioner_Parmetis<dim>::
set_balance_weights(SmartPtr<BalanceWeights<dim> > balanceWeights)
{
	m_balanceWeights = balanceWeights;
}

template<int dim>
void Partitioner_Parmetis<dim>::
set_connection_weights(SmartPtr<ConnectionWeights<dim> > conWeights)
{
	m_connectionWeights = conWeights;
}

template<int dim>
bool Partitioner_Parmetis<dim>::
supports_balance_weights() const
{
	return true;
}

template<int dim>
bool Partitioner_Parmetis<dim>::
supports_connection_weights() const
{
	return true;
}

//template<int dim>
//number Partitioner_Parmetis<dim>::
//estimate_distribution_quality()
//{
////todo	Consider connection weights in the final quality!
//	typedef typename Grid::traits<elem_t>::iterator ElemIter;
//	using std::min;
//
//	MultiGrid& mg = *m_mg;
//	DistributedGridManager& distGridMgr = *mg.distributed_grid_manager();
//
//	number minQuality = 1;
//
//	accumulate_child_counts(0, mg.top_level(), m_aNumChildren);
//
//	Table<stringstream> qualityOut;
//
////	calculate the quality estimate.
////todo The quality of a level could be weighted by the total amount of elements
////		in each level.
//	for(size_t hlvl = 0; hlvl < m_processHierarchy->num_hierarchy_levels(); ++hlvl)
//	{
//		number quality = 1;
//		size_t minLvl = m_processHierarchy->grid_base_level(hlvl);
//
//		if(minLvl > mg.top_level())
//			break;
//
//		int numProcs = m_processHierarchy->num_global_procs_involved(hlvl);
//		bool processParticipates = false;
//		pcl::ProcessCommunicator procComAll = m_processHierarchy->global_proc_com(hlvl);
//		if(!procComAll.empty()){
//			processParticipates = true;
//			if(numProcs > 1){
//				int localWeight = 0;
//				for(ElemIter iter = mg.begin<elem_t>(minLvl);
//					iter != mg.end<elem_t>(minLvl); ++iter)
//				{
//					if(!distGridMgr.is_ghost(*iter))
//						localWeight += (1 + m_aaNumChildren[*iter]);//todo: use balance weights
//				}
//
//				int minWeight = procComAll.allreduce(localWeight, PCL_RO_MIN);
//				int maxWeight = procComAll.allreduce(localWeight, PCL_RO_MAX);
//
//				if(maxWeight > 0)
//					quality = (number)minWeight / (number)maxWeight;
//				else
//					processParticipates = false;
//			}
//		}
//
//		minQuality = min(minQuality, quality);
//
//		if(verbose()){
//			qualityOut(hlvl + 1, 1) << minLvl;
//			qualityOut(hlvl + 1, 2) << numProcs;
//			if(processParticipates)
//				qualityOut(hlvl + 1, 3) << quality;
//			else
//				qualityOut(hlvl + 1, 3) << "idle";
//		}
//	}
//
//	if(verbose()){
//		qualityOut(0, 1) << "grid level";
//		qualityOut(0, 2) << "num procs";
//		qualityOut(0, 3) << "estimated quality";
//
//		UG_LOG("Estimated distribution quality:\n" << qualityOut << "\n");
//	}
//
////	the quality is a global property - we thus have to use the global minimum
//	pcl::ProcessCommunicator comGlobal;
//	return comGlobal.allreduce(minQuality, PCL_RO_MIN);
//}

template<int dim>
number Partitioner_Parmetis<dim>::
estimate_distribution_quality()
{
//todo	Consider connection weights in the final quality!
	typedef typename Grid::traits<elem_t>::iterator ElemIter;
	using std::min;

	MultiGrid& mg = *m_mg;
	DistributedGridManager& distGridMgr = *mg.distributed_grid_manager();

	number minQuality = 1;

	Table<stringstream> qualityOut;

//	calculate the quality in each level
	for(size_t lvl = 0; lvl < mg.num_levels(); ++lvl){
		size_t hlvl = m_processHierarchy->hierarchy_level_from_grid_level(lvl);
		int numProcs = m_processHierarchy->num_global_procs_involved(hlvl);
		bool processParticipates = false;
		pcl::ProcessCommunicator procComAll = m_processHierarchy->global_proc_com(hlvl);
		number quality = 1;
		if(!procComAll.empty()){
			processParticipates = true;
			if(numProcs > 1){
				int localWeight = 0;
				for(ElemIter iter = mg.begin<elem_t>(lvl);
					iter != mg.end<elem_t>(lvl); ++iter)
				{
					if(!distGridMgr.is_ghost(*iter))
						localWeight++;
				}

				int minWeight = procComAll.allreduce(localWeight, PCL_RO_MIN);
				int maxWeight = procComAll.allreduce(localWeight, PCL_RO_MAX);

				if(maxWeight > 0)
					quality = (number)minWeight / (number)maxWeight;
				else
					processParticipates = false;
			}
		}

		minQuality = min(minQuality, quality);

		if(verbose()){
			qualityOut(lvl + 1, 1) << lvl;
			qualityOut(lvl + 1, 2) << numProcs;
			if(processParticipates)
				qualityOut(lvl + 1, 3) << quality;
			else
				qualityOut(lvl + 1, 3) << "idle";
		}
	}

	if(verbose()){
		qualityOut(0, 1) << "grid level";
		qualityOut(0, 2) << "num procs";
		qualityOut(0, 3) << "estimated quality";

		UG_LOG("Estimated distribution quality:\n" << qualityOut << "\n");
	}

//	the quality is a global property - we thus have to use the global minimum
	pcl::ProcessCommunicator comGlobal;
	return comGlobal.allreduce(minQuality, PCL_RO_MIN);
}

template<int dim>
void Partitioner_Parmetis<dim>::
accumulate_child_counts(int baseLvl, int topLvl, AInt aInt,
						bool copyToVMastersOnBaseLvl)
{
	UG_DLOG(LIB_GRID, 1, "Partitioner_Parmetis-start accumulate_child_counts\n");
	typedef typename Grid::traits<elem_t>::iterator ElemIter;

	assert(m_mg);
	assert(m_mg->is_parallel());
	assert(m_mg->has_attachment<elem_t>(aInt));

	MultiGrid& mg = *m_mg;

	if((topLvl < baseLvl) || (baseLvl < 0) || topLvl >= (int)mg.num_levels()){
		UG_THROW("Bad levels supplied: baseLvl = " << baseLvl << ", topLvl = "
				 << topLvl << ", mg.num_levels() = " << mg.num_levels());
	}

	Grid::AttachmentAccessor<elem_t, AInt> aaNumChildren(mg, aInt);

	SetAttachmentValues(aaNumChildren, mg.begin<elem_t>(topLvl),
						mg.end<elem_t>(topLvl), 0);


	for(int lvl = topLvl - 1; lvl >= baseLvl; --lvl){
		for(ElemIter iter = mg.begin<elem_t>(lvl);
			iter != mg.end<elem_t>(lvl); ++iter)
		{
			elem_t* e = *iter;
			size_t numChildren = mg.num_children<elem_t>(e);
			aaNumChildren[e] = numChildren;
			for(size_t i_child = 0; i_child < numChildren; ++i_child){
				aaNumChildren[e] += aaNumChildren[mg.get_child<elem_t>(e, i_child)];
			}
		}

	//	since we perform load balancing on v-slaves, we don't have to inform
	//	v-masters on the base-lvl
		if(mg.is_parallel() && ((lvl > baseLvl) || copyToVMastersOnBaseLvl)){
			GridLayoutMap& glm = mg.distributed_grid_manager()->grid_layout_map();
		//	communicate the child counts from v-slaves to v-masters, since the
		//	latter havn't got children on their local processes.
			ComPol_CopyAttachment<layout_t, AInt> compolCopy(mg, aInt);
			if(glm.has_layout<elem_t>(INT_V_SLAVE)){
				m_intfcCom.send_data(glm.get_layout<elem_t>(INT_V_SLAVE).layout_on_level(lvl),
									 compolCopy);
			}
			if(glm.has_layout<elem_t>(INT_V_MASTER)){
				m_intfcCom.receive_data(glm.get_layout<elem_t>(INT_V_MASTER).layout_on_level(lvl),
										compolCopy);
			}
			m_intfcCom.communicate();
		}
	}
	UG_DLOG(LIB_GRID, 1, "Partitioner_Parmetis-stop accumulate_child_counts\n");
}

template<int dim>
void Partitioner_Parmetis<dim>::
partition(size_t baseLvl, size_t elementThreshold)
{
	UG_DLOG(LIB_GRID, 1, "Partitioner_Parmetis-start rebalance\n");

	typedef typename Grid::traits<elem_t>::iterator ElemIter;

	assert(m_mg);
	MultiGrid& mg = *m_mg;

//	The parallel dual graph is used by parmetis partitioning.
//	We don't initialize it now, since it only has to be initialized if
//	load-partitioning via parmetis is performed.
//	Declaring the graph here is for performance reasons only! It avoids
//	repeated attachment and detachment of required data in the graph.
	ParallelDualGraph<elem_t, idx_t> pdg;

	int localProc = pcl::GetProcRank();

	m_sh.clear();

//	assign all elements below baseLvl to the local process
	for(int i = 0; i < (int)baseLvl; ++i)
		m_sh.assign_subset(mg.begin<elem_t>(i), mg.end<elem_t>(i), localProc);

	//accumulate_child_counts(baseLvl, mg.top_level(), m_aNumChildren);

//	iterate through all hierarchy levels and perform rebalancing for all
//	hierarchy-sections which contain levels higher than baseLvl
	for(size_t hlevel = 0; hlevel < m_processHierarchy->num_hierarchy_levels(); ++ hlevel)
	{
		int minLvl = m_processHierarchy->grid_base_level(hlevel);
		int maxLvl = (int)mg.num_levels() - 1;
		if(hlevel + 1 < m_processHierarchy->num_hierarchy_levels()){
			maxLvl = min<int>(maxLvl,
						(int)m_processHierarchy->grid_base_level(hlevel + 1) - 1);
		}

		if(minLvl < (int)baseLvl)
			minLvl = (int)baseLvl;

		if(maxLvl < minLvl)
			continue;

		int numProcs = m_processHierarchy->num_global_procs_involved(hlevel);

		if(numProcs <= 1){
			for(int i = minLvl; i <= maxLvl; ++i)
				m_sh.assign_subset(mg.begin<elem_t>(i), mg.end<elem_t>(i), localProc);
			continue;
		}

	//	note that even if a process is yet used on a given hierarchy level, it may
	//	still contain some low dimensional dummy elements in h-interfaces. We thus
	//	continue execution on a process even if it is not originally involved in
	//	redistribution on this hlevel. Note that this is only a problem if
	//	maxLvl = minLvl + 1.
		pcl::ProcessCommunicator procComAll = m_processHierarchy->global_proc_com(hlevel);
		pcl::ProcessCommunicator globalCom;

	//	check whether there are enough elements to perform partitioning
		if(elementThreshold > 0){
			int numLocalElems = mg.num<elem_t>(minLvl);
			int numGlobalElems = globalCom.allreduce(numLocalElems, PCL_RO_SUM);

			if(numGlobalElems / numProcs < (int)elementThreshold){
				if(verbose()){
					UG_LOG("No partitioning performed for level " << minLvl
							<< ": Not enough elements.\n");
				}
			//	we can't perform partitioning on this hierarchy level.
			//	Simply assign all elements of this hierarchy level to the local proc.
				for(int i = minLvl; i <= maxLvl; ++i)
					m_sh.assign_subset(mg.begin<elem_t>(i), mg.end<elem_t>(i), localProc);
				continue;
			}
		}

		accumulate_child_counts(minLvl, maxLvl, m_aNumChildren);

	//	we have to find out how many of the target processes already contain a grid.
		int gotGrid = 0;
		if(mg.num<elem_t>(minLvl) > 0)
			gotGrid = 1;

		int numProcsWithGrid = globalCom.allreduce(gotGrid, PCL_RO_SUM);

		if(numProcsWithGrid == 0)
			continue;

		if(numProcsWithGrid == 1){
			if(gotGrid)
				partition_level_metis(minLvl, numProcs);
		}
		else{
			partition_level_parmetis(minLvl, numProcs, procComAll, pdg);
		}

	//	assign partitions to all children in this hierarchy level
		for(int lvl = minLvl; lvl < maxLvl; ++lvl){
			for(ElemIter iter = mg.begin<elem_t>(lvl); iter != mg.end<elem_t>(lvl); ++iter)
			{
				size_t numChildren = mg.num_children<elem_t>(*iter);
				int si = m_sh.get_subset_index(*iter);
				for(size_t i = 0; i < numChildren; ++i)
					m_sh.assign_subset(mg.get_child<elem_t>(*iter, i), si);
			}

			if(mg.is_parallel()){
				GridLayoutMap& glm = mg.distributed_grid_manager()->grid_layout_map();
			//	communicate partitions from v-masters to v-slaves, since v-slaves
			//	havn't got no parents on their procs.
				ComPol_Subset<layout_t>	compolSHCopy(m_sh, true);
				if(glm.has_layout<elem_t>(INT_V_MASTER)){
					m_intfcCom.send_data(glm.get_layout<elem_t>(INT_V_MASTER).layout_on_level(lvl),
										 compolSHCopy);
				}
				if(glm.has_layout<elem_t>(INT_V_SLAVE)){
					m_intfcCom.receive_data(glm.get_layout<elem_t>(INT_V_SLAVE).layout_on_level(lvl),
											compolSHCopy);
				}
				m_intfcCom.communicate();
			}
		}
	}

	UG_DLOG(LIB_GRID, 1, "Partitioner_Parmetis-stop rebalance\n");
}


template<int dim>
void Partitioner_Parmetis<dim>::
partition_level_metis(int lvl, int numTargetProcs)
{
	UG_DLOG(LIB_GRID, 1, "Partitioner_Parmetis-start partition_level_metis\n");
	typedef typename Grid::traits<elem_t>::iterator ElemIter;
	assert(m_mg);
	MultiGrid& mg = *m_mg;

//	here we'll store the dual graph
	vector<idx_t> adjacencyMapStructure;
	vector<idx_t> adjacencyMap;

	ConstructDualGraphMGLevel<elem_t, idx_t>(adjacencyMapStructure, adjacencyMap,
											 mg, lvl);

	if(adjacencyMap.empty())
		return;

 //note: using the option METIS_OPTION_DBGLVL could be useful for debugging.
	idx_t options[METIS_NOPTIONS];
	METIS_SetDefaultOptions(options);
	options[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY;
	options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
	options[METIS_OPTION_NUMBERING] = 0;
	//options[METIS_OPTION_CONTIG] = 1;	 //	request contiguous partitions

	idx_t nVrts = (idx_t)adjacencyMapStructure.size() - 1;
	idx_t nConstraints = 1;
	idx_t edgeCut;
	idx_t numParts = (idx_t)numTargetProcs;
	vector<idx_t> partitionMap(nVrts);

//todo: consider specified balance and connection weights!

//	create a weight map for the vertices based on the number of children+1
//	for each graph-vertex. This is not necessary, if we're already on the top level
	idx_t* pVrtSizeMap = NULL;
	vector<idx_t> vrtSizeMap;
	if(lvl < (int)mg.top_level()){
		vrtSizeMap.reserve(nVrts);
		for(ElemIter iter = mg.begin<elem_t>(lvl); iter != mg.end<elem_t>(lvl); ++iter)
			vrtSizeMap.push_back(CHILD_WEIGHT * m_aaNumChildren[*iter] + 1);

		assert((int)vrtSizeMap.size() == nVrts);
		pVrtSizeMap = &vrtSizeMap.front();
	}

//	under some circumstances, we have to cluster siblings. This requires special
//	edge weights...
	idx_t* pAdjWgts = NULL;
	vector<idx_t> adjwgts;
	if(base_class::clustered_siblings_enabled()){
	//	we have to weighten edges between siblings higher than between non-siblings
		vector<elem_t*> elems;
		for(ElemIter iter = mg.begin<elem_t>(lvl); iter != mg.end<elem_t>(lvl); ++iter)
			elems.push_back(*iter);

		adjwgts.resize(adjacencyMap.size(), 1);

		UG_ASSERT(adjacencyMapStructure.size() >= 2,
				  "There have to be at least 2 entries in the structure map if 1 vertex exists!");

		for(size_t ivrt = 0; ivrt < adjacencyMapStructure.size() - 1; ++ivrt){
			elem_t* e = elems[ivrt];
			size_t edgesBegin = adjacencyMapStructure[ivrt];
			size_t edgesEnd = adjacencyMapStructure[ivrt+1];

			for(size_t ie = edgesBegin; ie != edgesEnd; ++ie){
				elem_t* ce = elems[adjacencyMap[ie]];
				side_t* side = GetSharedSide(mg, e, ce);
				if(side && (mg.parent_type(side) == elem_t::BASE_OBJECT_ID)){
				//	e and ce should share the same parent, since the side shared
				//	by e and ce has a parent of the same type as e and ce
					adjwgts[ie] = SIBLING_WEIGHT;
				}
			}
		}

		pAdjWgts = &adjwgts.front();
	}

	UG_DLOG(LIB_GRID, 1, "CALLING METIS\n");
	int metisRet =	METIS_PartGraphKway(&nVrts, &nConstraints,
										&adjacencyMapStructure.front(),
										&adjacencyMap.front(),
										pVrtSizeMap, NULL, pAdjWgts,
										&numParts, NULL, NULL, options,
										&edgeCut, &partitionMap.front());
	UG_DLOG(LIB_GRID, 1, "METIS DONE\n");

	if(metisRet != METIS_OK){
		UG_THROW("METIS FAILED while partitioning the grid on level " << lvl);
	}

//	assign partition-subsets from graph-colors
	int counter = 0;
	for(ElemIter iter = mg.begin<elem_t>(lvl); iter != mg.end<elem_t>(lvl); ++iter)
		m_sh.assign_subset(*iter, partitionMap[counter++]);


//	clustered siblings help to ensure that all vertices which are connected to
//	a constrained vertex through are on the same process as the constrained vertex.
//	If only refinement is performed, it would be sufficient to only cluster
//	constrained siblings. However, coarsening would be rather complicated in that
//	case, since it is rather complicated to introduce constrained sibling elements if a
//	previously unconstrained sibling is not located on the same process...
	if(base_class::clustered_siblings_enabled()){
		//UG_LOG("NOTE: Clustering siblings during partitioning.\n");
		if(lvl > 0){
		//	put all children in the subset of the first one.
			for(ElemIter iter = mg.begin<elem_t>(lvl-1);
				iter != mg.end<elem_t>(lvl-1); ++iter)
			{
				elem_t* e = *iter;
				size_t numChildren = mg.num_children<elem_t>(e);
				if(numChildren > 1){
					int partition = m_sh.get_subset_index(mg.get_child<elem_t>(e, 0));
					for(size_t i = 1; i < numChildren; ++i){
						m_sh.assign_subset(mg.get_child<elem_t>(e, i), partition);
					}
				}
			}
		}
	}

	UG_DLOG(LIB_GRID, 1, "Partitioner_Parmetis-stop partition_level_metis\n");
}

template<int dim>
void Partitioner_Parmetis<dim>::
partition_level_parmetis(int lvl, int numTargetProcs,
						 const pcl::ProcessCommunicator& procComAll,
						 ParallelDualGraph<elem_t, idx_t>& pdg)
{
	UG_DLOG(LIB_GRID, 1, "Partitioner_Parmetis-start partition_level_parmetis\n");
	typedef typename Grid::traits<elem_t>::iterator ElemIter;
	assert(m_mg);
	MultiGrid& mg = *m_mg;

	int localProc = pcl::GetProcRank();

//	generate the parallel graph. H-Interface communication involved.
	pdg.set_grid(m_mg);
	pdg.generate_graph(lvl, procComAll);

	if(!procComAll.empty()){
		UG_DLOG(LIB_GRID, 2, "  parallel dual graph #vrts: " << (int)pdg.num_graph_vertices()
							<< ", #edges: " << (int)pdg.num_graph_edges() << "\n");

		pcl::ProcessCommunicator procCom = procComAll.
									create_sub_communicator(pdg.num_graph_vertices() > 0);
	//	parmetis would have problems otherwise. The node-offset-map now matches procCom.
		pdg.remove_empty_procs_from_node_offset_map();

		if(!procCom.empty()){
		//	partition the graph using parmetis
			//idx_t partOptions[3] = {1, 1, 0};
			idx_t partOptions[3]; partOptions[0] = 0;//default values
			//idx_t refineOptions[4] = {1, 0, 0, PARMETIS_PSR_UNCOUPLED};
			idx_t refineOptions[4]; refineOptions[0] = 0;
			idx_t nVrts = pdg.num_graph_vertices();
			idx_t nConstraints = 1;
			idx_t edgeCut;
			idx_t wgtFlag = 3;//vertices and edges are weighted
			idx_t numFlag = 0;
			idx_t numParts = (idx_t)numTargetProcs;
			vector<idx_t> partitionMap(nVrts);
			vector<real_t> tpwgts(numParts, 1. / (number)numParts);
			real_t ubvec = 1.05;
			real_t comVsRedistRation = 1000;

		//todo: consider specified balance and connection weights!

		//	create a weight map for the vertices based on the number of children+1
			idx_t* pVrtSizeMap = NULL;
			//vector<idx_t> vrtSizeMap(nVrts, 1);
			vector<idx_t> vrtSizeMap;
			vrtSizeMap.reserve(nVrts);
			for(ElemIter iter = mg.begin<elem_t>(lvl); iter != mg.end<elem_t>(lvl); ++iter){
				if(pdg.was_considered(*iter))
					vrtSizeMap.push_back(CHILD_WEIGHT * m_aaNumChildren[*iter] + 1);
			}

			idx_t* pAdjWgts = NULL;
			vector<idx_t> adjwgts;
			if(base_class::clustered_siblings_enabled()){
			//	we have to weighten edges between siblings higher than between non-siblings
				adjwgts.resize(pdg.num_graph_edges(), 1);
				idx_t* adjacencyMapStructure = pdg.adjacency_map_structure();
				idx_t* adjacencyMap = pdg.adjacency_map();//global
				int localOffset = pdg.parallel_offset_map()[procCom.get_local_proc_id()];

				for(idx_t ivrt = 0; ivrt < nVrts; ++ivrt){
					elem_t* e = pdg.get_element(ivrt);
					size_t edgesBegin = adjacencyMapStructure[ivrt];
					size_t edgesEnd = adjacencyMapStructure[ivrt+1];

					for(size_t ie = edgesBegin; ie != edgesEnd; ++ie){
					//	if the connected element doesn't lie on this process it
					//	shouldn't be a sibling anyways...
						int localInd = (int)adjacencyMap[ie] - localOffset;
						if((localInd < 0) || (localInd >= (int)pdg.num_graph_vertices()))
							continue;

						elem_t* ce = pdg.get_element(localInd);
						side_t* side = GetSharedSide(mg, e, ce);
						if(side && (mg.parent_type(side) == elem_t::BASE_OBJECT_ID)){
						//	e and ce should share the same parent, since the side shared
						//	by e and ce has a parent of the same type as e and ce
							adjwgts[ie] = SIBLING_WEIGHT;
						}
					}
				}

				pAdjWgts = &adjwgts.front();
			}

			assert((int)vrtSizeMap.size() == nVrts);
			pVrtSizeMap = &vrtSizeMap.front();

			MPI_Comm mpiCom = procCom.get_mpi_communicator();
			if((int)procCom.size() != numTargetProcs){
				UG_DLOG(LIB_GRID, 1, "Calling Parmetis_V3_PartKWay...");
				int metisRet =	ParMETIS_V3_PartKway(pdg.parallel_offset_map(),
													pdg.adjacency_map_structure(),
													pdg.adjacency_map(),
													pVrtSizeMap, pAdjWgts, &wgtFlag,
													&numFlag, &nConstraints,
													&numParts, &tpwgts.front(), &ubvec, partOptions,
													&edgeCut, &partitionMap.front(),
													&mpiCom);
				UG_DLOG(LIB_GRID, 1, "done\n");

				if(metisRet != METIS_OK){
					UG_THROW("ParMETIS_V3_PartKway failed on process " << localProc
							 << " while partitioning level " << lvl);
				}
			}
			else{
				UG_DLOG(LIB_GRID, 1, "Calling ParMETIS_V3_AdaptiveRepart...");
				int metisRet =	ParMETIS_V3_AdaptiveRepart(pdg.parallel_offset_map(),
													pdg.adjacency_map_structure(),
													pdg.adjacency_map(),
													pVrtSizeMap, pVrtSizeMap, pAdjWgts, &wgtFlag,
													&numFlag, &nConstraints,
													&numParts, &tpwgts.front(), &ubvec,
													&comVsRedistRation, refineOptions,
													&edgeCut, &partitionMap.front(),
													&mpiCom);
				UG_DLOG(LIB_GRID, 1, "done\n");
				if(metisRet != METIS_OK){
					UG_THROW("Parmetis_V3_RefineKWay failed on process " << localProc
							 << " while partitioning level " << lvl);
				}
			}

		//	assign partition-subsets from graph-colors
			int counter = 0;
			for(ElemIter iter = mg.begin<elem_t>(lvl); iter != mg.end<elem_t>(lvl); ++iter){
				if(pdg.was_considered(*iter))
					m_sh.assign_subset(*iter, partitionMap[counter++]);
			}
		}
	}

//	copy subset indices from vertical slaves to vertical masters,
//	since partitioning was only performed on vslaves
	GridLayoutMap& glm = mg.distributed_grid_manager()->grid_layout_map();
	ComPol_Subset<layout_t>	compolSHCopy(m_sh, true);

	if(glm.has_layout<elem_t>(INT_V_SLAVE))
		m_intfcCom.send_data(glm.get_layout<elem_t>(INT_V_SLAVE).layout_on_level(lvl),
							 compolSHCopy);
	if(glm.has_layout<elem_t>(INT_V_MASTER))
		m_intfcCom.receive_data(glm.get_layout<elem_t>(INT_V_MASTER).layout_on_level(lvl),
						 	 	compolSHCopy);
	m_intfcCom.communicate();

//	clustered siblings help to ensure that all vertices which are connected to
//	a constrained vertex through are on the same process as the constrained vertex.
//	If only refinement is performed, it would be sufficient to only cluster
//	constrained siblings. However, coarsening would be rather complicated in that
//	case, since it is rather complicated to introduce constrained sibling elements if a
//	previously unconstrained sibling is not located on the same process...
//todo: clustering should already be considered during graph-partitioning.
	if(base_class::clustered_siblings_enabled()){
		//UG_LOG("NOTE: Clustering siblings during partitioning.\n");
		if(lvl > 0){
		//	put all children in the subset of the first one.
			for(ElemIter iter = mg.begin<elem_t>(lvl-1);
				iter != mg.end<elem_t>(lvl-1); ++iter)
			{
				elem_t* e = *iter;
				size_t numChildren = mg.num_children<elem_t>(e);
				if(numChildren > 1){
					int partition = m_sh.get_subset_index(mg.get_child<elem_t>(e, 0));
					for(size_t i = 1; i < numChildren; ++i){
						m_sh.assign_subset(mg.get_child<elem_t>(e, i), partition);
					}
				}
			}
		}

	//	communicate the subset id's back from masters to slaves, since slaves don't
	//	have parents and thus can't adjust their sibling id's
		if(glm.has_layout<elem_t>(INT_V_MASTER))
			m_intfcCom.send_data(glm.get_layout<elem_t>(INT_V_MASTER).layout_on_level(lvl),
								 compolSHCopy);

		if(glm.has_layout<elem_t>(INT_V_SLAVE))
			m_intfcCom.receive_data(glm.get_layout<elem_t>(INT_V_SLAVE).layout_on_level(lvl),
									compolSHCopy);
		m_intfcCom.communicate();
	}

	UG_DLOG(LIB_GRID, 1, "Partitioner_Parmetis-stop partition_level_parmetis\n");
}

template<int dim>
SubsetHandler& Partitioner_Parmetis<dim>::
get_partitions()
{
	return m_sh;
}

template<int dim>
const std::vector<int>* Partitioner_Parmetis<dim>::
get_process_map() const
{
	return NULL;
}

template class Partitioner_Parmetis<1>;
template class Partitioner_Parmetis<2>;
template class Partitioner_Parmetis<3>;

}// end of namespace
