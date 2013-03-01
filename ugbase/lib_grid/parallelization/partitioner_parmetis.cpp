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
#include "util/parallel_dual_graph.h"

extern "C" {
	#include "metis.h"
	#include "parmetis.h"
}

using namespace std;

namespace ug{

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

template<int dim>
void Partitioner_Parmetis<dim>::
accumulate_child_counts(int baseLvl, int topLvl, AInt aInt)
{
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

		if(mg.is_parallel()){
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
}

template<int dim>
void Partitioner_Parmetis<dim>::
partition(size_t baseLvl, size_t elementThreshold)
{
	typedef typename Grid::traits<elem_t>::iterator ElemIter;

	assert(m_mg);
	MultiGrid& mg = *m_mg;

	int localProc = pcl::GetProcRank();

//	assign all elements below baseLvl to the local process
	for(int i = 0; i < (int)baseLvl; ++i)
		m_sh.assign_subset(mg.begin<elem_t>(i), mg.end<elem_t>(i), localProc);

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

		accumulate_child_counts(minLvl, maxLvl, m_aNumChildren);


		pcl::ProcessCommunicator procComAll = m_processHierarchy->global_proc_com(hlevel);

		if(procComAll.empty())
			continue;

	//	check whether there are enough elements to perform partitioning
		if(elementThreshold > 0){
			int numLocalElems = mg.num<elem_t>(minLvl);
			int numGlobalElems = procComAll.allreduce(numLocalElems, PCL_RO_SUM);

			if(numGlobalElems / numProcs < (int)elementThreshold){
			//	we can't perform partitioning on this hierarchy level.
			//	Simply assign all elements of this hierarchy level to the local proc.
				for(int i = minLvl; i <= maxLvl; ++i)
					m_sh.assign_subset(mg.begin<elem_t>(i), mg.end<elem_t>(i), localProc);
				continue;
			}
		}

	//	we have to find out how many of the target processes already contain a grid.
		int gotGrid = 0;
		if(mg.num<elem_t>(minLvl) > 0)
			gotGrid = 1;

		int numProcsWithGrid = procComAll.allreduce(gotGrid, PCL_RO_SUM);

		if(numProcsWithGrid == 0)
			continue;

		if(numProcsWithGrid == 1){
			partition_level_metis(minLvl, numProcs);
		}
		else{
			partition_level_parmetis(minLvl, numProcs, procComAll);
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
}


template<int dim>
void Partitioner_Parmetis<dim>::
partition_level_metis(int lvl, int numTargetProcs)
{
	typedef typename Grid::traits<elem_t>::iterator ElemIter;
	assert(m_mg);
	MultiGrid& mg = *m_mg;

//	here we'll store the dual graph
	vector<idx_t> adjacencyMapStructure;
	vector<idx_t> adjacencyMap;

	ConstructDualGraphMGLevel<elem_t, idx_t>(adjacencyMapStructure, adjacencyMap,
											 mg, lvl);

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
			vrtSizeMap.push_back(m_aaNumChildren[*iter] + 1);

		assert((int)vrtSizeMap.size() == nVrts);
		pVrtSizeMap = &vrtSizeMap.front();
	}

	UG_DLOG(LIB_GRID, 1, "CALLING METIS\n");
	int metisRet =	METIS_PartGraphKway(&nVrts, &nConstraints,
										&adjacencyMapStructure.front(),
										&adjacencyMap.front(),
										NULL, pVrtSizeMap, NULL,
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
}

template<int dim>
void Partitioner_Parmetis<dim>::
partition_level_parmetis(int lvl, int numTargetProcs,
						 const pcl::ProcessCommunicator& procComAll)
{
	typedef typename Grid::traits<elem_t>::iterator ElemIter;
	assert(m_mg);
	MultiGrid& mg = *m_mg;

	int localProc = pcl::GetProcRank();

//	here we'll store the dual graph
	vector<idx_t> adjacencyMapStructure;
	vector<idx_t> adjacencyMap;
	vector<idx_t> nodeOffsetMap;

	ConstructParallelDualGraphMGLevel<elem_t, idx_t>(adjacencyMapStructure,
											adjacencyMap, nodeOffsetMap,
											mg, lvl, procComAll);

	UG_DLOG(LIB_GRID, 2, "  parallel dual graph #vrts: " << (int)adjacencyMapStructure.size() - 1
						<< ", #edges: " << adjacencyMap.size() / 2 << "\n");

	pcl::ProcessCommunicator procCom = procComAll.
								create_sub_communicator(adjacencyMap.size() > 1);

//	partition the graph using parmetis
	idx_t options[3]; options[0] = 0;//default values
	idx_t nVrts = (idx_t)adjacencyMapStructure.size() - 1;
	idx_t nConstraints = 1;
	idx_t edgeCut;
	idx_t wgtFlag = 2;//only vertices are weighted
	idx_t numFlag = 0;
	idx_t numParts = (idx_t)numTargetProcs;
	vector<idx_t> partitionMap(nVrts);
	vector<real_t> tpwgts(numParts, 1. / (number)numParts);
	real_t ubvec = 1.05;

//todo: consider specified balance and connection weights!

//	create a weight map for the vertices based on the number of children+1
	idx_t* pVrtSizeMap = NULL;
	vector<idx_t> vrtSizeMap;
	vrtSizeMap.reserve(nVrts);
	for(ElemIter iter = mg.begin<elem_t>(lvl); iter != mg.end<elem_t>(lvl); ++iter)
		vrtSizeMap.push_back(m_aaNumChildren[*iter] + 1);

	assert((int)vrtSizeMap.size() == nVrts);
	pVrtSizeMap = &vrtSizeMap.front();

	if(!procCom.empty()){
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
			UG_THROW("PARMETIS FAILED on process " << localProc
					 << " while partitioning level " << lvl);
		}
	}

//	assign partition-subsets from graph-colors
//	make sure to ignore ghosts!
	assert(mg.is_parallel());
	DistributedGridManager& distGridMgr = *mg.distributed_grid_manager();
	int counter = 0;
	for(ElemIter iter = mg.begin<elem_t>(lvl); iter != mg.end<elem_t>(lvl); ++iter){
		if(!distGridMgr.is_ghost(*iter))
			m_sh.assign_subset(*iter, partitionMap[counter++]);
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
