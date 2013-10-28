// created by Sebastian Reiter
// s.b.reiter@gmail.com
// Feb 25, 2013 (d,m,y)

#include "partitioner_parmetis.h"
#include "load_balancer_util.h"
#include "distributed_grid.h"
#include "lib_grid/parallelization/util/compol_copy_attachment.h"
#include "lib_grid/parallelization/util/compol_subset.h"
#include "lib_grid/parallelization/parallelization_util.h"
#include "lib_grid/algorithms/attachment_util.h"
#include "lib_grid/algorithms/graph/dual_graph.h"
#include "lib_grid/algorithms/geom_obj_util/misc_util.h"
#include "common/util/table.h"
using namespace std;

namespace ug{

template <int dim>
Partitioner_Parmetis<dim>::
Partitioner_Parmetis() :
	m_mg(NULL),
	m_childWeight(1),
	m_siblingWeight(2),
	m_comVsRedistRatio(1000)
{
	m_processHierarchy = SPProcessHierarchy(new ProcessHierarchy);
	m_processHierarchy->add_hierarchy_level(0, 1);
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
set_next_process_hierarchy(SPProcessHierarchy procHierarchy)
{
	m_nextProcessHierarchy = procHierarchy;
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
ConstSPProcessHierarchy Partitioner_Parmetis<dim>::
current_process_hierarchy() const
{
	return m_processHierarchy;
}


template<int dim>
ConstSPProcessHierarchy Partitioner_Parmetis<dim>::
next_process_hierarchy() const
{
	return m_nextProcessHierarchy;
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
estimate_distribution_quality(std::vector<number>* pLvlQualitiesOut)
{
	GDIST_PROFILE_FUNC();
//todo	Consider connection weights in the final quality!
	typedef typename Grid::traits<elem_t>::iterator ElemIter;
	using std::min;

	MultiGrid& mg = *m_mg;
	DistributedGridManager& distGridMgr = *mg.distributed_grid_manager();

	number minQuality = 1;

	Table<stringstream> qualityOut;

	if(pLvlQualitiesOut)
		pLvlQualitiesOut->clear();

	const ProcessHierarchy* procH;
	if(m_nextProcessHierarchy.valid())
		procH = m_nextProcessHierarchy.get();
	else
		procH = m_processHierarchy.get();

//	calculate the quality in each level
	for(size_t lvl = 0; lvl < mg.num_levels(); ++lvl){
		size_t hlvl = procH->hierarchy_level_from_grid_level(lvl);
		int numProcs = procH->num_global_procs_involved(hlvl);
		bool processParticipates = false;
		pcl::ProcessCommunicator procComAll = procH->global_proc_com(hlvl);
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

		if(pLvlQualitiesOut){
			if(processParticipates)
				pLvlQualitiesOut->push_back(quality);
			else
				pLvlQualitiesOut->push_back(-1);
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
	GDIST_PROFILE_FUNC();
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


	for(int lvl = topLvl-1; lvl >= baseLvl; --lvl){
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


template <int dim>
void Partitioner_Parmetis<dim>::
gather_child_counts_from_level(int baseLvl, int childLvl, AInt aInt,
								bool copyToVMastersOnBaseLvl)
{
	GDIST_PROFILE_FUNC();
	UG_DLOG(LIB_GRID, 1, "Partitioner_Parmetis-start gather_child_counts_from_level\n");
	typedef typename Grid::traits<elem_t>::iterator ElemIter;

	if(childLvl <= baseLvl)
		return;

	assert(m_mg);
	assert(m_mg->is_parallel());
	assert(m_mg->has_attachment<elem_t>(aInt));

	MultiGrid& mg = *m_mg;
	Grid::AttachmentAccessor<elem_t, AInt> aaNumChildren(mg, aInt);
	GridLayoutMap& glm = mg.distributed_grid_manager()->grid_layout_map();
	ComPol_CopyAttachment<layout_t, AInt> compolCopy(mg, aInt);

//	first write child-counts to elements on childLvl - 1
	for(ElemIter iter = mg.begin<elem_t>(childLvl - 1);
		iter != mg.end<elem_t>(childLvl - 1); ++iter)
	{
		elem_t* e = *iter;
		aaNumChildren[e] = mg.num_children<elem_t>(e);
	}

	for(int lvl = childLvl - 2; lvl >= baseLvl; --lvl){
	//	copy from v-slaves to vmasters
		if(glm.has_layout<elem_t>(INT_V_SLAVE))
			m_intfcCom.send_data(glm.get_layout<elem_t>(INT_V_SLAVE).layout_on_level(lvl + 1),
								 compolCopy);
		if(glm.has_layout<elem_t>(INT_V_MASTER))
			m_intfcCom.receive_data(glm.get_layout<elem_t>(INT_V_MASTER).layout_on_level(lvl + 1),
									compolCopy);
		m_intfcCom.communicate();

	//	accumulate child counts in parent elements on lvl
		for(ElemIter iter = mg.begin<elem_t>(lvl); iter != mg.end<elem_t>(lvl); ++iter)
		{
			elem_t* e = *iter;
			aaNumChildren[e] = 0;
			size_t numChildren = mg.num_children<elem_t>(e);
			for(size_t i = 0; i < numChildren; ++i)
				aaNumChildren[e] += aaNumChildren[mg.get_child<elem_t>(e, i)];
		}
	}

//	if child counts are required in vmasters on the base-level, copy them now...
	if(copyToVMastersOnBaseLvl){
		if(glm.has_layout<elem_t>(INT_V_SLAVE))
			m_intfcCom.send_data(glm.get_layout<elem_t>(INT_V_SLAVE).layout_on_level(baseLvl),
								 compolCopy);
		if(glm.has_layout<elem_t>(INT_V_MASTER))
			m_intfcCom.receive_data(glm.get_layout<elem_t>(INT_V_MASTER).layout_on_level(baseLvl),
									compolCopy);
		m_intfcCom.communicate();
	}
	UG_DLOG(LIB_GRID, 1, "Partitioner_Parmetis-stop gather_child_counts_from_level\n");
}


template<int dim>
template<class TIter>
void Partitioner_Parmetis<dim>::
fill_multi_constraint_vertex_weights(vector<idx_t>& vwgtsOut,
									 int baseLvl, int maxLvl, AInt aInt,
									 bool fillVMastersOnBaseLvl,
									 TIter baseElemsBegin, TIter baseElemsEnd,
									 int numBaseElements)
{
	GDIST_PROFILE_FUNC();
	MultiGrid& mg = *m_mg;
	Grid::AttachmentAccessor<elem_t, AInt> aaNumChildren(mg, aInt);

	int ncons = maxLvl - baseLvl + 1;
	vwgtsOut.resize(numBaseElements * ncons);

//	write 1 into the base-level constraints
	int wgtInd = 0;
	const int wgtOffset = maxLvl - baseLvl + 1;

	for(TIter iter = baseElemsBegin; iter != baseElemsEnd; ++iter, wgtInd += wgtOffset)
	{
		vwgtsOut[wgtInd] = 1;
	}

	for(int ci = 1; ci < ncons; ++ci){
		gather_child_counts_from_level(baseLvl, baseLvl + ci, aInt, fillVMastersOnBaseLvl);
		wgtInd = ci;
		for(TIter iter = baseElemsBegin; iter != baseElemsEnd; ++iter, wgtInd += wgtOffset)
		{
			vwgtsOut[wgtInd] = aaNumChildren[*iter];
		}
	}
}


template<int dim>
void Partitioner_Parmetis<dim>::
partition(size_t baseLvl, size_t elementThreshold)
{
	GDIST_PROFILE_FUNC();
	UG_DLOG(LIB_GRID, 1, "Partitioner_Parmetis-start partition\n");

//	UG_LOG("Partitioning - current process hierarchy:\n");
//	UG_LOG(m_processHierarchy->to_string() << endl);
//	if(m_nextProcessHierarchy.valid()){
//		UG_LOG("Partitioning - next process hierarchy:\n");
//		UG_LOG(m_nextProcessHierarchy->to_string() << endl);
//	}

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

	GridLayoutMap& glm = mg.distributed_grid_manager()->grid_layout_map();
	UG_ASSERT(!(glm.has_layout<elem_t>(INT_H_MASTER) || glm.has_layout<elem_t>(INT_H_SLAVE)),
			  "Elements for which a partitioning is generated may not be contained in"
			  " H-Interfaces!");

//	assign all elements below baseLvl to the local process
	for(int i = 0; i < (int)baseLvl; ++i)
		m_sh.assign_subset(mg.begin<elem_t>(i), mg.end<elem_t>(i), localProc);

//	UG_LOG("New Partitioning. Local elements (including ghosts): " << mg.num<elem_t>() << endl);

	const ProcessHierarchy* procH;
	if(m_nextProcessHierarchy.valid())
		procH = m_nextProcessHierarchy.get();
	else
		procH = m_processHierarchy.get();

//	iterate through all hierarchy levels and perform rebalancing for all
//	hierarchy-sections which contain levels higher than baseLvl
	pcl::ProcessCommunicator globCom;
	for(size_t hlevel = 0; hlevel < procH->num_hierarchy_levels(); ++ hlevel)
	{
//		UG_LOG("h-level: " << hlevel << endl);
	//	make sure that certain processes don't get ahead of others...
		if(hlevel > 0)
			globCom.barrier();

		int minLvl = procH->grid_base_level(hlevel);
		if(minLvl > (int)mg.top_level())
			break;

		int maxLvl = mg.top_level();
		if(hlevel + 1 < procH->num_hierarchy_levels()){
			maxLvl = min<int>(maxLvl,
						(int)procH->grid_base_level(hlevel + 1) - 1);
		}

		if(minLvl < (int)baseLvl)
			minLvl = (int)baseLvl;

		if(maxLvl < minLvl)
			break;


		int numProcsOnLvl = procH->num_global_procs_involved(hlevel);

		if(numProcsOnLvl <= 1){
			int targetProcId = 0;
			for(int i = minLvl; i <= maxLvl; ++i)
				m_sh.assign_subset(mg.begin<elem_t>(i), mg.end<elem_t>(i), targetProcId);
			continue;
		}


//	//	since the old and new hierarchy may be different, we'll get the old hlevel
//	//	using the current grid level
//		size_t gridLvl = procH->grid_base_level(hlevel);
//		size_t oldHLvl = m_processHierarchy.hierarchy_level_from_grid_level(gridLvl);
//		pcl::ProcessCommunicator hlvlCom = m_processHierarchy.global_proc_com(oldHLvl);
//		if(hlvlCom.empty()){
//			int oldMinLvl = m_processHierarchy.grid_base_level(oldHLvl);
//			if(oldMinLvl < (int)baseLvl)
//				oldMinLvl = (int)baseLvl;
//			int oldMaxLvl = mg.top_level();
//			if(oldHLvl + 1 < m_processHierarchy.num_hierarchy_levels()){
//				oldMaxLvl = min<int>(oldMaxLvl,
//							(int)m_processHierarchy.grid_base_level(oldHLvl + 1) - 1);
//			}
//
//		//	the local level is not contained in the current hlvl of the proc-hierarchy.
//		//	make sure that it doesn't contain any elements, since we would most likely
//		//	get problems during h- or v-interface communication later on.
//			UG_COND_THROW(mg.num<VertexBase>(oldMinLvl) > 0,
//					  "Process " << localProc
//					  << " shouldn't contain vertices on this level: " << oldMinLvl);
//
//		//	since we communicate with elements in maxLvl below, this assertion has to hold.
//			UG_COND_THROW(mg.num<elem_t>(oldMaxLvl) > 0,
//					  "Process " << localProc
//					  << " shouldn't contain elements on this level: " << oldMaxLvl);
//			continue;
//		}

	//	create a communicator for all processes which participate in load-balancing
	//	on this h-level. If the old and new process hierarchies match, this is the same
	//	as the old-hierarchies global-proc-com. If they don't match, we have to communicate...

		pcl::ProcessCommunicator hlvlCom;
		if(m_nextProcessHierarchy.valid()){
		//	check if this process contains vertices between minLvl and maxLvl
			bool containsVrts = false;
			for(int i = minLvl; i <= maxLvl; ++i){
				if(mg.num<VertexBase>(i) != 0){
					containsVrts = true;
					continue;
				}
			}

			pcl::ProcessCommunicator globalCom;
			hlvlCom = globalCom.create_sub_communicator(containsVrts);

			if(!containsVrts)
				continue;
		}
		else{
			hlvlCom = m_processHierarchy->global_proc_com(hlevel);
			if(hlvlCom.empty()){
			//	the local level is not contained in the current hlvl of the proc-hierarchy.
			//	make sure that it doesn't contain any elements, since we would most likely
			//	get problems during h- or v-interface communication later on.
				UG_COND_THROW(mg.num<VertexBase>(minLvl) > 0,
						  "Process " << localProc
						  << " shouldn't contain vertices on this level: " << minLvl);

			//	since we communicate with elements in maxLvl below, this assertion has to hold.
				UG_COND_THROW(mg.num<elem_t>(maxLvl) > 0,
						  "Process " << localProc
						  << " shouldn't contain elements on this level: " << maxLvl);
				continue;
			}
		}

	//	adjust the number of target processes so that each receives at least
	//	#elementThreshold elements.
		int numGhosts = 0;
		if(glm.has_layout<elem_t>(INT_V_MASTER)){
			numGhosts = glm.get_layout<elem_t>(INT_V_MASTER).
							layout_on_level(minLvl).num_interface_elements();
		}

		int numLocalElems = (int)mg.num<elem_t>(minLvl) - numGhosts;
		int numTargetProcs = numProcsOnLvl;
		if(elementThreshold > 0){
			int numGlobalElems = hlvlCom.allreduce(numLocalElems, PCL_RO_SUM);
			numTargetProcs = clip<int>(numGlobalElems / elementThreshold, 1, numProcsOnLvl);
		}

//		UG_LOG("  numTargetProcs: " << numTargetProcs << endl);
		if(numTargetProcs == 1){
		//	all elements have to be sent to process 0
			for(int i = minLvl; i <= maxLvl; ++i)
				m_sh.assign_subset(mg.begin<elem_t>(i), mg.end<elem_t>(i), 0);
			continue;
		}

	//	we have to find out how many of the level's processes already contain a grid.
		int gotGrid = 0;
		if(numLocalElems > 0)
			gotGrid = 1;

		int numProcsWithGrid = hlvlCom.allreduce(gotGrid, PCL_RO_SUM);

//		UG_LOG("  numProcsWithGrid: " << numProcsWithGrid << endl);

		if(numProcsWithGrid == 0)
			break;

		if(numProcsWithGrid == 1){
			//if(gotGrid)
			partition_level_metis(minLvl, maxLvl, numTargetProcs);
		}
		else{
			partition_level_parmetis(minLvl, maxLvl, numTargetProcs, hlvlCom, pdg);
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
					m_intfcCom.send_data(glm.get_layout<elem_t>(INT_V_MASTER).layout_on_level(lvl+1),
										 compolSHCopy);
				}
				if(glm.has_layout<elem_t>(INT_V_SLAVE)){
					m_intfcCom.receive_data(glm.get_layout<elem_t>(INT_V_SLAVE).layout_on_level(lvl+1),
											compolSHCopy);
				}
				m_intfcCom.communicate();
			}
		}
	}

	if(m_nextProcessHierarchy.valid()){
		*m_processHierarchy = *m_nextProcessHierarchy;
		m_nextProcessHierarchy = SPProcessHierarchy(NULL);
	}

	UG_DLOG(LIB_GRID, 1, "Partitioner_Parmetis-stop partition\n");
}


template<int dim>
void Partitioner_Parmetis<dim>::
partition_level_metis(int baseLvl, int maxLvl, int numTargetProcs)
{
	GDIST_PROFILE_FUNC();
	UG_DLOG(LIB_GRID, 1, "Partitioner_Parmetis-start partition_level_metis\n");
	typedef typename Grid::traits<elem_t>::iterator ElemIter;
	assert(m_mg);
	MultiGrid& mg = *m_mg;
	int lvl = baseLvl;

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
	idx_t edgeCut;
	idx_t numParts = (idx_t)numTargetProcs;
	vector<idx_t> partitionMap(nVrts);

//todo: consider specified balance and connection weights!
	idx_t nConstraints = maxLvl - baseLvl + 1;
	vector<idx_t> vrtWgtMap;
	fill_multi_constraint_vertex_weights(vrtWgtMap, baseLvl, maxLvl, m_aNumChildren, false,
										 mg.begin<elem_t>(baseLvl), mg.end<elem_t>(baseLvl),
										 mg.num<elem_t>(baseLvl));

//	create a weight map for the vertices based on the number of children+1
//	for each graph-vertex. This is not necessary, if we're already on the top level
//	accumulate_child_counts(baseLvl, maxLvl, m_aNumChildren);
//	if(lvl < (int)mg.top_level()){
//		vrtWgtMap.reserve(nVrts);
//		for(ElemIter iter = mg.begin<elem_t>(lvl); iter != mg.end<elem_t>(lvl); ++iter){
//			vrtSizeMap.push_back(m_childWeight * m_aaNumChildren[*iter] + 1);
//		}
//
//		assert((int)vrtSizeMap.size() == nVrts);
//		pVrtWgtMap = &vrtWgtMap.front();
//	}

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
					adjwgts[ie] = m_siblingWeight;
				}
			}
		}

		pAdjWgts = &adjwgts.front();
	}

	if(!adjacencyMap.empty()){
		UG_DLOG(LIB_GRID, 1, "CALLING METIS\n");
		GDIST_PROFILE(METIS);
		int metisRet =	METIS_PartGraphKway(&nVrts, &nConstraints,
											&adjacencyMapStructure.front(),
											&adjacencyMap.front(),
											&vrtWgtMap.front(), NULL, pAdjWgts,
											&numParts, NULL, NULL, options,
											&edgeCut, &partitionMap.front());
		GDIST_PROFILE_END();
		UG_DLOG(LIB_GRID, 1, "METIS DONE\n");

		if(metisRet != METIS_OK){
			UG_THROW("METIS FAILED while partitioning the grid on level " << lvl);
		}
	}

//	assign partition-subsets from graph-colors
	int counter = 0;
	for(ElemIter iter = mg.begin<elem_t>(lvl); iter != mg.end<elem_t>(lvl); ++iter){
		m_sh.assign_subset(*iter, partitionMap[counter++]);
	}

//	copy subset indices from vertical slaves to vertical masters,
//	since partitioning was only performed on vslaves
//	this is only important if e.g. process 0 has only ghosts and process 1 has all slaves.
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
partition_level_parmetis(int baseLvl, int maxLvl, int numTargetProcs,
						 const pcl::ProcessCommunicator& procComAll,
						 ParallelDualGraph<elem_t, idx_t>& pdg)
{
	GDIST_PROFILE_FUNC();

	PCL_DEBUG_BARRIER(procComAll);

	UG_DLOG(LIB_GRID, 1, "Partitioner_Parmetis-start partition_level_parmetis\n");
	typedef typename Grid::traits<elem_t>::iterator ElemIter;
	assert(m_mg);
	MultiGrid& mg = *m_mg;
	int lvl = baseLvl;

//	generate the parallel graph. H-Interface communication involved.
	pdg.set_grid(m_mg);
	pdg.generate_graph(lvl, procComAll);

	idx_t partOptions[3];
	partOptions[0] = 0;//default values
//	partOptions[0] = 1;
//	partOptions[1] = 127;
//	partOptions[2] = 1;

	//idx_t refineOptions[4] = {1, 0, 0, PARMETIS_PSR_UNCOUPLED};
	idx_t nConstraints = maxLvl - baseLvl + 1;
	idx_t refineOptions[4]; refineOptions[0] = 0;
	idx_t nVrts = pdg.num_graph_vertices();
	idx_t edgeCut;
	idx_t wgtFlag = 3;//vertices and edges are weighted
	idx_t numFlag = 0;
	idx_t numParts = (idx_t)numTargetProcs;
	vector<idx_t> partitionMap(nVrts);
	vector<real_t> tpwgts(numParts * nConstraints, 1. / (number)numParts);
	vector<real_t> ubvec(nConstraints, 1.05);
	real_t comVsRedistRation = m_comVsRedistRatio;

//todo: consider specified balance and connection weights!
	vector<idx_t> vrtWgtMap;
	fill_multi_constraint_vertex_weights(vrtWgtMap, baseLvl, maxLvl, m_aNumChildren, false,
										 pdg.elements_begin(), pdg.elements_end(),
										 pdg.num_graph_vertices());

//	(THE CODE BELOW WAS REPLACED BY THE CODE ABOVE...)
//	create a weight map for the vertices based on the number of children+1
//	accumulate_child_counts(baseLvl, maxLvl, m_aNumChildren);
//	vrtWgtMap.reserve(nVrts);
//	for(ElemIter iter = mg.begin<elem_t>(lvl); iter != mg.end<elem_t>(lvl); ++iter){
//		if(pdg.was_considered(*iter))
//			vrtWgtMap.push_back(m_childWeight * m_aaNumChildren[*iter] + 1);
//	}

	UG_DLOG(LIB_GRID, 2, "  parallel dual graph #vrts: " << (int)pdg.num_graph_vertices()
						<< ", #edges: " << (int)pdg.num_graph_edges() << "\n");

	pcl::ProcessCommunicator procCom = pdg.process_communicator();

	if(!procCom.empty()){
		idx_t* pAdjWgts = NULL;
		vector<idx_t> adjwgts;
		adjwgts.resize(pdg.num_graph_edges(), 1);
		if(base_class::clustered_siblings_enabled()){
		//	we have to weighten edges between siblings higher than between non-siblings
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
						adjwgts[ie] = m_siblingWeight;
					}
				}
			}
		}
		pAdjWgts = &adjwgts.front();

//		UG_LOG("parallel offset map: ");
//		for(size_t i = 0; i <= procCom.size(); ++i){
//			UG_LOG(" " << pdg.parallel_offset_map()[i]);
//		}
//		UG_LOG(endl);
//		UG_LOG("vertex edge offsets: ");
//		for(size_t i = 0; i <= pdg.num_graph_vertices(); ++i){
//			UG_LOG(" " << pdg.adjacency_map_structure()[i]);
//		}
//		UG_LOG(endl);
//		UG_LOG("connected vertices: ");
//		for(size_t i = 0; i < pdg.num_graph_edges(); ++i){
//			UG_LOG(" " << pdg.adjacency_map()[i]);
//		}
//		UG_LOG(endl);

		GDIST_PROFILE(PARMETIS);
		MPI_Comm mpiCom = procCom.get_mpi_communicator();
		if((int)procCom.size() != numTargetProcs){
			//UG_DLOG(LIB_GRID, 1, "Calling Parmetis_V3_PartKWay...");
//			UG_LOG("Calling Parmetis_V3_PartKWay...\n");
			int metisRet =	ParMETIS_V3_PartKway(pdg.parallel_offset_map(),
												pdg.adjacency_map_structure(),
												pdg.adjacency_map(),
												&vrtWgtMap.front(), pAdjWgts, &wgtFlag,
												&numFlag, &nConstraints,
												&numParts, &tpwgts.front(), &ubvec.front(), partOptions,
												&edgeCut, &partitionMap.front(),
												&mpiCom);
			UG_DLOG(LIB_GRID, 1, "done\n");

			if(metisRet != METIS_OK){
				UG_THROW("ParMETIS_V3_PartKway failed on process " << pcl::GetProcRank()
						 << " while partitioning level " << lvl);
			}
		}
		else{
		//	defines the redistribution cost for each vertex
		//	Since vrtWgtMap holds the number of children for each constraint on each level,
		//	we can simply sum them up to retrieve the redistribution cost for that element
		//todo: Better distributions may possibly be obtained if vrtSizeMap would contain 1
		//		for all elements... Check that out!
			vector<idx_t> vrtSizeMap;
			vrtSizeMap.reserve(pdg.num_graph_vertices());
			for(int ivrt = 0; ivrt < (int)vrtWgtMap.size(); ivrt += nConstraints){
				idx_t size = 0;
				for(int ic = 0; ic < (int)nConstraints; ++ic){
					size += vrtWgtMap[ivrt + ic];
				}
				vrtSizeMap.push_back(size);
			}


			//UG_DLOG(LIB_GRID, 1, "Calling ParMETIS_V3_AdaptiveRepart...");
//			UG_LOG("Calling ParMETIS_V3_AdaptiveRepart...\n");
			int metisRet =	ParMETIS_V3_AdaptiveRepart(pdg.parallel_offset_map(),
												pdg.adjacency_map_structure(),
												pdg.adjacency_map(),
												&vrtWgtMap.front(), &vrtSizeMap.front(), pAdjWgts, &wgtFlag,
												&numFlag, &nConstraints,
												&numParts, &tpwgts.front(), &ubvec.front(),
												&comVsRedistRation, refineOptions,
												&edgeCut, &partitionMap.front(),
												&mpiCom);
			UG_DLOG(LIB_GRID, 1, "done\n");
			if(metisRet != METIS_OK){
				UG_THROW("Parmetis_V3_RefineKWay failed on process " << pcl::GetProcRank()
						 << " while partitioning level " << lvl);
			}
		}
		GDIST_PROFILE_END();

	//	assign partition-subsets from graph-colors
		int counter = 0;
		for(ElemIter iter = mg.begin<elem_t>(lvl); iter != mg.end<elem_t>(lvl); ++iter){
			if(pdg.was_considered(*iter))
				m_sh.assign_subset(*iter, partitionMap[counter++]);
		}
	}

//	UG_LOG("PARMETIS DONE\n");

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
		GDIST_PROFILE(cluster_siblings);
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

template<int dim>
void Partitioner_Parmetis<dim>::
set_child_weight(int w)
{
	m_childWeight = w;
}

template<int dim>
void Partitioner_Parmetis<dim>::
set_sibling_weight(int w)
{
	m_siblingWeight = w;
}

template<int dim>
void Partitioner_Parmetis<dim>::
set_itr_factor(float itr)
{
	m_comVsRedistRatio = itr;
}

template class Partitioner_Parmetis<1>;
template class Partitioner_Parmetis<2>;
template class Partitioner_Parmetis<3>;

}// end of namespace
