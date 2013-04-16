// created by Sebastian Reiter
// s.b.reiter@gmail.com
// Mar 1, 2013 (d,m,y)

#include "partitioner_bisection.h"
#include "load_balancing.h"
#include "distributed_grid.h"

using namespace std;

namespace ug{

template <int dim>
Partitioner_Bisection<dim>::
Partitioner_Bisection() :
	m_mg(NULL),
	m_highestRedistLevel(-1)
{
	m_processHierarchy = SPProcessHierarchy(new ProcessHierarchy);
}

template<int dim>
Partitioner_Bisection<dim>::
~Partitioner_Bisection()
{
}

template<int dim>
void Partitioner_Bisection<dim>::
set_grid(MultiGrid* mg, Attachment<MathVector<dim> > aPos)
{
	m_mg = mg;
	m_sh.assign_grid(m_mg);
	m_aPos = aPos;
}

template<int dim>
void Partitioner_Bisection<dim>::
set_process_hierarchy(SPProcessHierarchy procHierarchy)
{
	m_processHierarchy = procHierarchy;
}

template<int dim>
void Partitioner_Bisection<dim>::
set_balance_weights(SmartPtr<BalanceWeights<dim> >)
{
}

template<int dim>
void Partitioner_Bisection<dim>::
set_connection_weights(SmartPtr<ConnectionWeights<dim> >)
{
}

template<int dim>
bool Partitioner_Bisection<dim>::
supports_balance_weights() const
{
	return false;
}

template<int dim>
bool Partitioner_Bisection<dim>::
supports_connection_weights() const
{
	return false;
}

template<int dim>
number Partitioner_Bisection<dim>::
estimate_distribution_quality()
{
//todo	Consider connection weights in the final quality!
	typedef typename Grid::traits<elem_t>::iterator ElemIter;
	using std::min;

	MultiGrid& mg = *m_mg;
	DistributedGridManager& distGridMgr = *mg.distributed_grid_manager();

	number minQuality = 1;

//	calculate the quality estimate.
//todo The quality of a level could be weighted by the total amount of elements
//		in each level.
	for(size_t lvl = 0; lvl < mg.num_levels(); ++lvl){
		size_t hlvl = m_processHierarchy->hierarchy_level_from_grid_level(lvl);
		int numProcs = m_processHierarchy->num_global_procs_involved(hlvl);
		if(numProcs <= 1)
			continue;

		pcl::ProcessCommunicator procComAll = m_processHierarchy->global_proc_com(hlvl);
		if(!procComAll.empty()){
			int localWeight = 0;
			for(ElemIter iter = mg.begin<elem_t>(lvl);
				iter != mg.end<elem_t>(lvl); ++iter)
			{
				if(!distGridMgr.is_ghost(*iter))
					localWeight += 1;
			}

			int minWeight = procComAll.allreduce(localWeight, PCL_RO_MIN);
			int maxWeight = procComAll.allreduce(localWeight, PCL_RO_MAX);

			if(maxWeight > 0){
				number quality = (number)minWeight / (number)maxWeight;
				minQuality = min(minQuality, quality);
			}
		}
	}

	pcl::ProcessCommunicator comGlobal;
	return comGlobal.allreduce(minQuality, PCL_RO_MIN);
}

template<int dim>
void Partitioner_Bisection<dim>::
partition(size_t baseLvl, size_t elementThreshold)
{
	typedef typename Grid::traits<elem_t>::iterator ElemIter;

	assert(m_mg);
	MultiGrid& mg = *m_mg;

	m_sh.clear();

//	assign all elements below baseLvl to the local process
	for(int i = 0; i < (int)baseLvl; ++i)
		m_sh.assign_subset(mg.begin<elem_t>(i), mg.end<elem_t>(i), 0);

//	iterate through all hierarchy levels and perform rebalancing for all
//	hierarchy-sections which contain levels higher than baseLvl
	m_procMap.clear();
	for(size_t hlevel = 0; hlevel < m_processHierarchy->num_hierarchy_levels(); ++ hlevel)
	{
		int minLvl = m_processHierarchy->grid_base_level(hlevel);
		int maxLvl = (int)mg.top_level();
		if(hlevel + 1 < m_processHierarchy->num_hierarchy_levels()){
			maxLvl = min<int>(maxLvl,
						(int)m_processHierarchy->grid_base_level(hlevel + 1) - 1);
		}

		if(minLvl < (int)baseLvl)
			minLvl = (int)baseLvl;

		if(maxLvl < minLvl)
			continue;

		const std::vector<int>& clusterProcs = m_processHierarchy->cluster_procs(hlevel);

		int numProcs = (int)clusterProcs.size();

		if(numProcs <= 1){
			for(int i = minLvl; i <= maxLvl; ++i)
				m_sh.assign_subset(mg.begin<elem_t>(i), mg.end<elem_t>(i), 0);
			continue;
		}

	//	we'll only perform partitioning, if we havn't done so on this level already.
	//	this is of course a severe restriction, however, currently no parallel
	//	bisection-partitioning algorithm exist in ug.
		if((int)hlevel <= m_highestRedistLevel){
			/*static bool warningPrinted = false;
			if(!warningPrinted){
				UG_ERR_LOG("WARNING: Partitioner_Bisection not used as intended. You should "
						   "call rebalance on your LoadBalancer whenever you reach a new "
						   "distribution-level during refinement. This is due to the fact "
						   "that a serial bisection partitioner is used even for parallel "
						   "bisections - which of course imposes some restrictions on the "
						   "distribution process.\n "
						   "NOTE: you may alternatively use the Partitioner_Parmetis, which "
						   "does not have those restrictions. To do so, make sure to compile "
						   "ug with Parmetis support (-DPARMETIS=ON).\n");
				warningPrinted = true;
			}*/

			for(int i = minLvl; i <= maxLvl; ++i)
				m_sh.assign_subset( mg.begin<elem_t>(i), mg.end<elem_t>(i), 0);
			continue;
		}

	//	check whether there are enough elements to perform partitioning
		if(elementThreshold > 0){
			int numLocalElems = mg.num<elem_t>(minLvl);
			if(numLocalElems / numProcs < (int)elementThreshold){
			//	we can't perform partitioning on this hierarchy level.
			//	Simply assign all elements of this hierarchy level to the local proc.
				for(int i = minLvl; i <= maxLvl; ++i)
					m_sh.assign_subset(mg.begin<elem_t>(i), mg.end<elem_t>(i), 0);
				continue;
			}
		}

		PartitionElementsByRepeatedIntersection<elem_t, dim>(m_sh, mg, minLvl,
															 numProcs, m_aPos);

	//	assign partitions to all children in this hierarchy level
		for(int lvl = minLvl; lvl < maxLvl; ++lvl){
			for(ElemIter iter = mg.begin<elem_t>(lvl); iter != mg.end<elem_t>(lvl); ++iter)
			{
				size_t numChildren = mg.num_children<elem_t>(*iter);
				int si = m_sh.get_subset_index(*iter);
				for(size_t i = 0; i < numChildren; ++i)
					m_sh.assign_subset(mg.get_child<elem_t>(*iter, i), si);
			}
		}

		m_procMap = clusterProcs;
		m_highestRedistLevel = hlevel;
	}

	if(m_procMap.empty() && (m_sh.num_subsets() > 0)){
		if(m_sh.num_subsets() != 1){
			UG_THROW("Something went wrong during partitioning. At this point"
					" either exactly one subset or a filled process map should exist.");
		}
		m_procMap.push_back(pcl::GetProcRank());
	}

//	make sure that everybody knows about the highestRedistLevel!
	pcl::ProcessCommunicator com;
	m_highestRedistLevel = com.allreduce(m_highestRedistLevel, PCL_RO_MAX);
}


template<int dim>
SubsetHandler& Partitioner_Bisection<dim>::
get_partitions()
{
	return m_sh;
}

template<int dim>
const std::vector<int>* Partitioner_Bisection<dim>::
get_process_map() const
{
	return &m_procMap;
}

template class Partitioner_Bisection<1>;
template class Partitioner_Bisection<2>;
template class Partitioner_Bisection<3>;

}// end of namespace
