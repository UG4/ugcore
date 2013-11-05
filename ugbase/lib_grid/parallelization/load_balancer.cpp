// created by Sebastian Reiter
// s.b.reiter@gmail.com
// Feb 25, 2013 (d,m,y)

#include <algorithm>
#include "load_balancer.h"
#include "load_balancer_util.h"
#include "distribution.h"
#include "partitioner_bisection.h"
#include "distributed_grid.h"
#include "lib_grid/parallelization/parallelization_util.h"
#include "common/util/table.h"

#ifdef UG_PARMETIS
#include "partitioner_parmetis.h"
#endif

using namespace std;

namespace ug{

ProcessHierarchy::~ProcessHierarchy()
{
}

void ProcessHierarchy::
add_hierarchy_level(size_t gridLvl, size_t numProcsPerProc)
{
	if(numProcsPerProc == 0){
		UG_THROW("a hierarchy level has to contain at least 1 process");
	}

	if(!m_levels.empty()){
		if(gridLvl <= m_levels.back().gridLvl){
			UG_THROW("A grid level in higher hierarchy levels has to be higher"
					 " than in the levels below.");
		}
	}
//	make sure that the lowest hierarchy level cares for grid-level 0
	if(m_levels.empty() && (gridLvl != 0)){
		add_hierarchy_level(0, 1);
	}

	size_t curLvl = m_levels.size();
	m_levels.resize(m_levels.size() + 1);
	HLevelInfo& hlevel = m_levels.back();

	hlevel.gridLvl = gridLvl;
	hlevel.numGlobalProcsInUse = numProcsPerProc;
	if(curLvl >= 1){
		hlevel.numGlobalProcsInUse *= m_levels[curLvl - 1].numGlobalProcsInUse;
	}
	//hlevel.clusterCom = create_cluster_communicator(curLvl, gridLvl, numProcsPerProc);
	hlevel.globalCom =  pcl::ProcessCommunicator::create_communicator(0, hlevel.numGlobalProcsInUse);
	init_cluster_procs(hlevel.clusterProcs, curLvl, numProcsPerProc);
}

/*
pcl::ProcessCommunicator ProcessHierarchy::
create_cluster_communicator(size_t hlvl, size_t gridLvl, size_t numProcsPerProc)
{
	if(hlvl == 0){
		return pcl::ProcessCommunicator::create_communicator(0, numProcsPerProc);
	}

	const HLevelInfo& parentLvl = get_hlevel_info(hlvl - 1);

	if(numProcsPerProc == 1){
		return parentLvl.clusterCom;
	}

//	calculate the root process for this cluster and create the group based on rootProc
	int rootProc = pcl::GetProcRank() / (int)parentLvl.numGlobalProcsInUse;

	std::vector<int> procs;
	procs.reserve(numProcsPerProc);
	procs.push_back(rootProc);
	int firstNewProc = (int)parentLvl.numGlobalProcsInUse + rootProc * ((int)numProcsPerProc - 1);

	for(int i = 0; i < (int)numProcsPerProc - 1; ++i){
		procs.push_back(firstNewProc + i);
	}

	return pcl::ProcessCommunicator::create_communicator(procs);
}*/

void ProcessHierarchy::
init_cluster_procs(std::vector<int>& clusterProcs, size_t hlvl, size_t numProcsPerProc)
{
	UG_ASSERT(numProcsPerProc > 0, "each proc has to distribute to at least 1 new process");

	clusterProcs.clear();
	if(hlvl == 0){
		clusterProcs.reserve(numProcsPerProc);
		for(size_t i = 0; i < numProcsPerProc; ++i){
			clusterProcs.push_back(i);
		}
		return;
	}

	const HLevelInfo& parentLvl = get_hlevel_info(hlvl - 1);

	if(numProcsPerProc == 1){
		clusterProcs = parentLvl.clusterProcs;
		return;
	}


//	calculate the root process for this cluster and create the group based on rootProc
	int localProc = pcl::GetProcRank();
	int rootProc = localProc;
	if(localProc >= (int)parentLvl.numGlobalProcsInUse)
		rootProc = (localProc - (int)parentLvl.numGlobalProcsInUse) / ((int)numProcsPerProc - 1);

	clusterProcs.push_back(rootProc);
	int firstNewProc = (int)parentLvl.numGlobalProcsInUse + rootProc * ((int)numProcsPerProc - 1);

	for(int i = 0; i < (int)numProcsPerProc - 1; ++i){
		clusterProcs.push_back(firstNewProc + i);
	}
}

bool ProcessHierarchy::
empty() const
{
	return m_levels.empty();
}

size_t ProcessHierarchy::
num_hierarchy_levels() const
{
	return m_levels.size();
}

size_t ProcessHierarchy::
num_global_procs_involved(size_t hierarchyLevel) const
{
	return get_hlevel_info(hierarchyLevel).numGlobalProcsInUse;
}

size_t ProcessHierarchy::
grid_base_level(size_t hierarchyLevel) const
{
	return get_hlevel_info(hierarchyLevel).gridLvl;
}

size_t ProcessHierarchy::
hierarchy_level_from_grid_level(size_t gridLevel) const
{
	if(m_levels.empty()){
		UG_THROW("No hierarchy levels exist.");
	}
	if(m_levels[0].gridLvl != 0){
		UG_THROW("The lowest grid level has to be 0!");
	}

	for(size_t i = 0; i + 1 < m_levels.size(); ++i){
		if(m_levels[i+1].gridLvl > gridLevel)
			return i;
	}
	return m_levels.size() - 1;
}

pcl::ProcessCommunicator ProcessHierarchy::
global_proc_com(size_t hierarchyLevel) const
{
	return get_hlevel_info(hierarchyLevel).globalCom;
}

const std::vector<int>& ProcessHierarchy::
cluster_procs(size_t hierarchyLevel) const
{
	return get_hlevel_info(hierarchyLevel).clusterProcs;
}

std::string ProcessHierarchy::
to_string() const
{
	StringStreamTable table;
	table(0, 0) << "lvl:";
	table(1, 0) << "procs:";
	for(size_t i = 0; i < m_levels.size(); ++i){
		table(0, i + 1) << m_levels[i].gridLvl;
		table(1, i + 1) << m_levels[i].numGlobalProcsInUse;
	}

	return table.to_string();
}

//pcl::ProcessCommunicator ProcessHierarchy::
//cluster_proc_com(size_t hierarchyLevel)
//{
//	return get_hlevel_info(hierarchyLevel).clusterCom;
//}



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template<int dim>
LoadBalancer<dim>::
LoadBalancer() :
	m_mg(NULL),
	m_balanceThreshold(0.9),
	m_elementThreshold(1),
	m_createVerticalInterfaces(true)
{
	m_processHierarchy = ProcessHierarchy::create();
	m_balanceWeights = SPBalanceWeights(new StdBalanceWeights<dim>);
	m_connectionWeights = SPConnectionWeights(new StdConnectionWeights<dim>);


	#ifdef UG_PARMETIS
		m_partitioner = SPPartitioner(new Partitioner_Parmetis<dim>);
	#else
		m_partitioner = SPPartitioner(new Partitioner_Bisection<dim>);
	#endif
}

template<int dim>
LoadBalancer<dim>::
~LoadBalancer()
{
}

template<int dim>
void LoadBalancer<dim>::
set_grid(MultiGrid* mg, Attachment<MathVector<dim> > aPos)
{
	m_mg = mg;
	m_aPos = aPos;
}

template<int dim>
void LoadBalancer<dim>::
enable_vertical_interface_creation(bool enable)
{
	m_createVerticalInterfaces = enable;
}

template<int dim>
void LoadBalancer<dim>::
set_partitioner(SmartPtr<IPartitioner<dim> > partitioner)
{
	m_partitioner = partitioner;
}

template<int dim>
void LoadBalancer<dim>::
set_balance_weights(SmartPtr<BalanceWeights<dim> > balanceWeights)
{
	m_balanceWeights = balanceWeights;
}

template<int dim>
void LoadBalancer<dim>::
set_connection_weights(SmartPtr<ConnectionWeights<dim> > conWeights)
{
	m_connectionWeights = conWeights;
}

//template<int dim>
//void LoadBalancer<dim>::
//add_distribution_level(size_t lvl, size_t numProcsPerProc)
//{
//	m_processHierarchy->add_hierarchy_level(lvl, numProcsPerProc);
//}

template<int dim>
void LoadBalancer<dim>::
set_next_process_hierarchy(SPProcessHierarchy procHierarchy)
{
	m_processHierarchy = procHierarchy;
}

template<int dim>
void LoadBalancer<dim>::
set_balance_threshold(number threshold)
{
	m_balanceThreshold = threshold;
}

template<int dim>
void LoadBalancer<dim>::
set_element_threshold(size_t threshold)
{
	m_elementThreshold = threshold;
}

//template<int dim>
//number LoadBalancer<dim>::
//distribution_quality()
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
////	calculate the quality estimate.
////todo The quality of a level could be weighted by the total amount of elements
////		in each level.
//	for(size_t lvl = mg.top_level(); lvl < mg.num_levels(); ++lvl){
//		size_t hlvl = m_processHierarchy->hierarchy_level_from_grid_level(lvl);
//		int numProcs = m_processHierarchy->num_global_procs_involved(hlvl);
//		if(numProcs <= 1)
//			continue;
//
//		pcl::ProcessCommunicator procComAll = m_processHierarchy->global_proc_com(hlvl);
//
//		int localWeight = 0;
//		for(ElemIter iter = mg.begin<elem_t>(lvl);
//			iter != mg.end<elem_t>(lvl); ++iter)
//		{
//			if(!distGridMgr.is_ghost(*iter))
//				localWeight += 1;//todo: use balance weights
//		}
//
//		int minWeight = procComAll.allreduce(localWeight, PCL_RO_MIN);
//		int maxWeight = procComAll.allreduce(localWeight, PCL_RO_MAX);
//
//		if(maxWeight <= 0)
//			continue;
//
//		number quality = (number)minWeight / (number)maxWeight;
//
//		minQuality = min(minQuality, quality);
//	}
//
//	return minQuality;
//}


template<int dim>
void LoadBalancer<dim>::
rebalance()
{
	GDIST_PROFILE_FUNC();

	UG_DLOG(LIB_GRID, 1, "LoadBalancer-start rebalance\n");
	UG_COND_THROW(!m_partitioner.valid(),
				  "LoadBalancer::rebalance can only be performed with a valid partitioner!");

	UG_COND_THROW(m_processHierarchy->empty(),
				  "A Process-Hierarchy has to be specifed for rebalancing");

	m_balanceWeights->set_grid(m_mg, m_aPos);
	m_connectionWeights->set_grid(m_mg, m_aPos);
	m_partitioner->set_grid(m_mg, m_aPos);
	m_partitioner->set_next_process_hierarchy(m_processHierarchy);
	m_partitioner->set_connection_weights(m_connectionWeights);
	m_partitioner->set_balance_weights(m_balanceWeights);

//todo:	check imbalance and find base-level on which to partition!
	m_balanceWeights->refresh_weights(0);
	m_connectionWeights->refresh_weights(0);

//	distribution quality is only interesting if repartitioning is supported.
//	If it is not we'll set it to -1, thus calling partition anyways
	number distQuality = -1;
	if(m_partitioner->supports_repartitioning()){
		distQuality = m_partitioner->estimate_distribution_quality();
		if(!m_partitioner->verbose()){
			UG_LOG("Current estimated distribution quality: " << distQuality << "\n");
		}
	}

	if(m_balanceThreshold > distQuality)
	{
		UG_LOG("Redistributing...\n");
		UG_DLOG(LIB_GRID, 1, "LoadBalancer-rebalance: partitioning...\n");
		m_partitioner->partition(0, m_elementThreshold);

		SubsetHandler& sh = m_partitioner->get_partitions();
		if(sh.num<elem_t>() != m_mg->num<elem_t>()){
			UG_THROW("All elements have to be assigned to subsets during partitioning! "
					 << "Please check your partitioner!");
		}

		const std::vector<int>* procMap = m_partitioner->get_process_map();

		UG_DLOG(LIB_GRID, 1, "LoadBalancer-rebalance: distributing...\n");
		if(!DistributeGrid(*m_mg, sh, m_serializer, m_createVerticalInterfaces, procMap))
		{
			UG_THROW("DistributeGrid failed!");
		}

		UG_LOG("Redistribution done\n");
		number newDistQuality = m_partitioner->estimate_distribution_quality();
		if(!m_partitioner->verbose()){
			UG_LOG("Estimated distribution quality after redistribution: " << newDistQuality << "\n");
		}
	}
	else{
		UG_LOG("No redistribution necessary.\n");
	}
	UG_DLOG(LIB_GRID, 1, "LoadBalancer-stop rebalance\n");
}

template<int dim>
void LoadBalancer<dim>::
create_quality_record(const char* label)
{
	std::vector<number>	lvlQualities;
	bool isVerbose = m_partitioner->verbose();
	m_partitioner->set_verbose(false);
	m_partitioner->estimate_distribution_quality(&lvlQualities);
	m_partitioner->set_verbose(isVerbose);
	ConstSPProcessHierarchy procH = m_partitioner->current_process_hierarchy();

//	fill header first
	if(m_qualityRecords(0, 0).str().empty()){
		m_qualityRecords(0, 0) << "level:";
	}
	for(size_t i = 0; i < lvlQualities.size(); ++i){
		if(m_qualityRecords(0, 2*i + 1).str().empty())
			m_qualityRecords(0, 2*i + 1) << i;
	}

	if(procH.valid()){
		size_t ri = max<size_t>(1, m_qualityRecords.num_rows());
		m_qualityRecords(ri, 0) << label;
		for(size_t i = 0; i < lvlQualities.size(); ++i){
			size_t hlvl = procH->hierarchy_level_from_grid_level(i);

			if(i == procH->grid_base_level(hlvl)){
			//	redistribution takes place on this level
				m_qualityRecords(ri, 2*i+1) << "(p" << procH->num_global_procs_involved(hlvl) << ")";
			}
			else
				m_qualityRecords(ri, 2*i+1) << "-";
			m_qualityRecords(ri, 2*i+2) << lvlQualities[i];
		}
	}
}

template<int dim>
void LoadBalancer<dim>::
print_quality_records() const
{
	UG_LOG(m_qualityRecords << "\n");
}

template class LoadBalancer<1>;
template class LoadBalancer<2>;
template class LoadBalancer<3>;

} // end of namespace
