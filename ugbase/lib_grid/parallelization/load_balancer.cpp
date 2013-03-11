// created by Sebastian Reiter
// s.b.reiter@gmail.com
// Feb 25, 2013 (d,m,y)

#include "load_balancer.h"
#include "load_balancer_util.h"
#include "distribution.h"
#include "partitioner_bisection.h"

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
		rootProc = localProc / (int)parentLvl.numGlobalProcsInUse;

	clusterProcs.push_back(rootProc);
	int firstNewProc = (int)parentLvl.numGlobalProcsInUse + rootProc * ((int)numProcsPerProc - 1);

	for(int i = 0; i < (int)numProcsPerProc - 1; ++i){
		clusterProcs.push_back(firstNewProc + i);
	}
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
	m_processHierarchy = SPProcessHierarchy(new ProcessHierarchy);
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

template<int dim>
void LoadBalancer<dim>::
add_distribution_level(size_t lvl, size_t numProcsPerProc)
{
	m_processHierarchy->add_hierarchy_level(lvl, numProcsPerProc);
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

template<int dim>
void LoadBalancer<dim>::
rebalance()
{
	UG_DLOG(LIB_GRID, 1, "LoadBalancer-start rebalance\n");
	if(!m_partitioner.valid()){
		UG_THROW("LoadBalancer::rebalance can only be performed with a valid partitioner!");
	}

	m_balanceWeights->set_grid(m_mg, m_aPos);
	m_connectionWeights->set_grid(m_mg, m_aPos);
	m_partitioner->set_grid(m_mg, m_aPos);
	m_partitioner->set_process_hierarchy(m_processHierarchy);
	m_partitioner->set_connection_weights(m_connectionWeights);
	m_partitioner->set_balance_weights(m_balanceWeights);

//todo:	check imbalance and find base-level on which to partition!
	m_balanceWeights->refresh_weights(0);
	m_connectionWeights->refresh_weights(0);
	m_partitioner->partition(0, m_elementThreshold);

	SubsetHandler& sh = m_partitioner->get_partitions();
	if(sh.num<elem_t>() != m_mg->num<elem_t>()){
		UG_THROW("All elements have to be assigned to subsets during partitioning! "
				 << "Please check your partitioner!");
	}

	const std::vector<int>* procMap = m_partitioner->get_process_map();

	if(!DistributeGrid(*m_mg, sh, m_serializer, m_createVerticalInterfaces, procMap))
	{
		UG_THROW("DistributeGrid failed!");
	}
	UG_DLOG(LIB_GRID, 1, "LoadBalancer-stop rebalance\n");
}


template class LoadBalancer<1>;
template class LoadBalancer<2>;
template class LoadBalancer<3>;

} // end of namespace
