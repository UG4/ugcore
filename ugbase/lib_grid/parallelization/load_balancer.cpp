// created by Sebastian Reiter
// s.b.reiter@gmail.com
// Feb 25, 2013 (d,m,y)

#include <algorithm>
#include "load_balancer.h"
#include "load_balancer_util.h"
#include "distribution.h"
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
	int rootProc = pcl::ProcRank() / (int)parentLvl.numGlobalProcsInUse;

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
		//clusterProcs = parentLvl.clusterProcs;
		clusterProcs.push_back(pcl::ProcRank());
		return;
	}


//	calculate the root process for this cluster and create the group based on rootProc
	int localProc = pcl::ProcRank();
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
LoadBalancer::
LoadBalancer() :
	m_mg(NULL),
	m_balanceThreshold(0.9),
	m_elementThreshold(1),
	m_createVerticalInterfaces(true)
{
	m_processHierarchy = ProcessHierarchy::create();
	m_balanceWeights = make_sp(new StdBalanceWeights());
}

LoadBalancer::
~LoadBalancer()
{
}

void LoadBalancer::
set_grid(MultiGrid* mg)
{
	m_mg = mg;
}

void LoadBalancer::
enable_vertical_interface_creation(bool enable)
{
	m_createVerticalInterfaces = enable;
}

void LoadBalancer::
set_partitioner(SmartPtr<IPartitioner> partitioner)
{
	m_partitioner = partitioner;
}

void LoadBalancer::
set_balance_weights(SmartPtr<IBalanceWeights> balanceWeights)
{
	m_balanceWeights = balanceWeights;
}
//
//template<int dim>
//void LoadBalancer::
//set_connection_weights(SmartPtr<ConnectionWeights<dim> > conWeights)
//{
//	m_connectionWeights = conWeights;
//}

//template<int dim>
//void LoadBalancer::
//add_distribution_level(size_t lvl, size_t numProcsPerProc)
//{
//	m_processHierarchy->add_hierarchy_level(lvl, numProcsPerProc);
//}

void LoadBalancer::
set_next_process_hierarchy(SPProcessHierarchy procHierarchy)
{
	m_processHierarchy = procHierarchy;
}

void LoadBalancer::
set_balance_threshold(number threshold)
{
	m_balanceThreshold = threshold;
}

void LoadBalancer::
set_element_threshold(size_t threshold)
{
	m_elementThreshold = threshold;
}

bool LoadBalancer::
problems_occurred()
{
	if(m_partitioner.valid())
		return m_partitioner->problems_occurred();
	return false;
}

number LoadBalancer::
estimate_distribution_quality(std::vector<number>* pLvlQualitiesOut)
{
	if(m_mg){
		int highestElem = VERTEX;
		if(m_mg->num<Volume>() > 0)		highestElem = VOLUME;
		else if(m_mg->num<Face>() > 0)	highestElem = FACE;
		else if(m_mg->num<Edge>() > 0)	highestElem = EDGE;

		pcl::ProcessCommunicator procCom;
		highestElem = procCom.allreduce(highestElem, PCL_RO_MAX);

		switch(highestElem){
		case VERTEX:
			return estimate_distribution_quality_impl<Vertex>(pLvlQualitiesOut);
		case EDGE:
			return estimate_distribution_quality_impl<Edge>(pLvlQualitiesOut);
		case FACE:
			return estimate_distribution_quality_impl<Face>(pLvlQualitiesOut);
		case VOLUME:
			return estimate_distribution_quality_impl<Volume>(pLvlQualitiesOut);
		}
	}

	if(pLvlQualitiesOut)
		pLvlQualitiesOut->clear();

	return 0;
}

template <class TElem>
number LoadBalancer::
estimate_distribution_quality_impl(std::vector<number>* pLvlQualitiesOut)
{
//todo	Consider connection weights in the final quality!
	typedef TElem elem_t;
	typedef typename Grid::traits<elem_t>::iterator ElemIter;
	using std::min;

	if(m_balanceWeights.invalid())
		m_balanceWeights = make_sp(new IBalanceWeights());

	MultiGrid& mg = *m_mg;
	DistributedGridManager& distGridMgr = *mg.distributed_grid_manager();

	number minQuality = 1;

	if(pLvlQualitiesOut)
		pLvlQualitiesOut->clear();

//	calculate the quality estimate.
//todo The quality of a level could be weighted by the total amount of elements
//		in each level.
	const ProcessHierarchy* procH = m_processHierarchy.get();

	for(size_t lvl = 0; lvl < mg.num_levels(); ++lvl){
		size_t hlvl = procH->hierarchy_level_from_grid_level(lvl);
		int numProcs = procH->num_global_procs_involved(hlvl);
		if(numProcs <= 1){
			if(pLvlQualitiesOut)
				pLvlQualitiesOut->push_back(1.0);
			continue;
		}

		pcl::ProcessCommunicator procComAll = procH->global_proc_com(hlvl);
		if(procComAll.size() == 0){
			pLvlQualitiesOut->push_back(-1);
		}
		else if(procComAll.size() == 1){
			pLvlQualitiesOut->push_back(1);
		}
		else{
			number localWeight = 0;
			IBalanceWeights& wgts = *m_balanceWeights;
			for(ElemIter iter = mg.begin<elem_t>(lvl);
				iter != mg.end<elem_t>(lvl); ++iter)
			{
				if(!distGridMgr.is_ghost(*iter))
					localWeight += wgts.get_weight(*iter);
			}

			number maxW = procComAll.allreduce(localWeight, PCL_RO_MAX);
			//number minW = procComAll.allreduce(localWeight, PCL_RO_MIN);
			number totalW = procComAll.allreduce(localWeight, PCL_RO_SUM);
			size_t numProcs = procComAll.size();

			if(maxW > 0){
				number quality = (totalW - maxW) / (maxW * number(numProcs - 1));
				minQuality = min(minQuality, quality);
				if(pLvlQualitiesOut)
					pLvlQualitiesOut->push_back(quality);
			}
			else if(pLvlQualitiesOut)
				pLvlQualitiesOut->push_back(-1);
		}
	}

	pcl::ProcessCommunicator comGlobal;
	return comGlobal.allreduce(minQuality, PCL_RO_MIN);
}

number LoadBalancer::
estimate_distribution_quality()
{
	std::vector<number> v;
	return estimate_distribution_quality(&v);
}

bool LoadBalancer::
rebalance()
{
	GDIST_PROFILE_FUNC();

	UG_DLOG(LIB_GRID, 1, "LoadBalancer-start rebalance\n");
	UG_COND_THROW(!m_partitioner.valid(),
				  "LoadBalancer::rebalance can only be performed with a valid partitioner!");

	UG_COND_THROW(m_processHierarchy->empty(),
				  "A Process-Hierarchy has to be specifed for rebalancing");

	m_partitioner->set_next_process_hierarchy(m_processHierarchy);
	//m_partitioner->set_connection_weights(m_connectionWeights);
	m_partitioner->set_balance_weights(m_balanceWeights);

//todo:	check imbalance and find base-level on which to partition!
	m_balanceWeights->refresh_weights(0);
	//m_connectionWeights->refresh_weights(0);

//	distribution quality is only interesting if repartitioning is supported.
//	If it is not we'll set it to -1, thus calling partition anyways
	number distQuality = -1;
	if(m_partitioner->supports_repartitioning()){
		distQuality = estimate_distribution_quality();
		if(!m_partitioner->verbose()){
			UG_LOG("Current estimated distribution quality: " << distQuality << "\n");
		}
	}

	if(m_balanceThreshold > distQuality)
	{
		UG_DLOG(LIB_GRID, 1, "LoadBalancer-rebalance: partitioning...\n");
		if(m_partitioner->partition(0, m_elementThreshold)){
			UG_LOG("Redistributing...\n");
			SubsetHandler& sh = m_partitioner->get_partitions();
//			if(sh.num<elem_t>() != m_mg->num<elem_t>()){
//				UG_THROW("All elements have to be assigned to subsets during partitioning! "
//						 << "Please check your partitioner!");
//			}

			const std::vector<int>* procMap = m_partitioner->get_process_map();

			UG_DLOG(LIB_GRID, 1, "LoadBalancer-rebalance: distributing...\n");
			if(!DistributeGrid(*m_mg, sh, m_serializer, m_createVerticalInterfaces, procMap))
			{
				UG_THROW("DistributeGrid failed!");
			}

			UG_LOG("Redistribution done\n");
			UG_DLOG(LIB_GRID, 1, "LoadBalancer-stop rebalance\n");
			return true;
		}
	}
	else{
		UG_LOG("No redistribution necessary.\n");
		UG_DLOG(LIB_GRID, 1, "LoadBalancer-stop rebalance\n");
		return true;
	}

	UG_DLOG(LIB_GRID, 1, "LoadBalancer-stop rebalance (FAILED)\n");
	return false;
}

void LoadBalancer::
create_quality_record(const char* label)
{
	std::vector<number>	lvlQualities;
	estimate_distribution_quality(&lvlQualities);
	ConstSPProcessHierarchy procH = m_processHierarchy;

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

void LoadBalancer::
print_quality_records() const
{
	UG_LOG(m_qualityRecords << "\n");
}

} // end of namespace
