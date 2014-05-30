// created by Sebastian Reiter
// s.b.reiter@gmail.com
// Nov 4, 2013

#include "partitioner_dynamic_bisection.h"
#include "load_balancer_util.h"
#include "distributed_grid.h"
#include "lib_grid/parallelization/util/compol_copy_attachment.h"
#include "lib_grid/parallelization/util/compol_subset.h"
#include "lib_grid/parallelization/parallelization_util.h"
#include "lib_grid/algorithms/attachment_util.h"
#include "lib_grid/algorithms/subset_util.h"

using namespace std;

namespace ug{

template <class TElem, int dim>
Partitioner_DynamicBisection<TElem, dim>::
Partitioner_DynamicBisection() :
	m_mg(NULL),
	m_staticPartitioning(false),
	m_tolerance(0.99),
	m_splitImproveIterations(10),
	m_highestRedistLevel(-1)
{
	m_processHierarchy = SPProcessHierarchy(new ProcessHierarchy);
	m_processHierarchy->add_hierarchy_level(0, 1);

	m_balanceWeights = make_sp(new IBalanceWeights());
}

template <class TElem, int dim>
Partitioner_DynamicBisection<TElem, dim>::
~Partitioner_DynamicBisection()
{
}

template <class TElem, int dim>
void Partitioner_DynamicBisection<TElem, dim>::
set_grid(MultiGrid* mg, Attachment<MathVector<dim> > aPos)
{
	m_mg = mg;
	if(m_sh.valid())
		m_sh->assign_grid(m_mg);
	m_aPos = aPos;
	m_aaPos.access(*m_mg, m_aPos);
}

template <class TElem, int dim>
void Partitioner_DynamicBisection<TElem, dim>::
set_subset_handler(SmartPtr<SubsetHandler> sh)
{
	m_sh = sh;
	if(m_mg)
		m_sh->assign_grid(m_mg);
}

template <class TElem, int dim>
void Partitioner_DynamicBisection<TElem, dim>::
set_next_process_hierarchy(SPProcessHierarchy procHierarchy)
{
	m_nextProcessHierarchy = procHierarchy;
}

template <class TElem, int dim>
void Partitioner_DynamicBisection<TElem, dim>::
set_balance_weights(SPBalanceWeights balanceWeights)
{
	m_balanceWeights = balanceWeights;
}
//
//template <class TElem, int dim>
//void Partitioner_DynamicBisection<TElem, dim>::
//set_connection_weights(SmartPtr<ConnectionWeights<dim> >)
//{
//}


template <class TElem, int dim>
ConstSPProcessHierarchy Partitioner_DynamicBisection<TElem, dim>::
current_process_hierarchy() const
{
	return m_processHierarchy;
}


template <class TElem, int dim>
ConstSPProcessHierarchy Partitioner_DynamicBisection<TElem, dim>::
next_process_hierarchy() const
{
	return m_nextProcessHierarchy;
}


template <class TElem, int dim>
bool Partitioner_DynamicBisection<TElem, dim>::
supports_balance_weights() const
{
	return false;
}

template <class TElem, int dim>
bool Partitioner_DynamicBisection<TElem, dim>::
supports_connection_weights() const
{
	return false;
}

template <class TElem, int dim>
SubsetHandler& Partitioner_DynamicBisection<TElem, dim>::
get_partitions()
{
	if(m_sh.invalid()){
		if(m_mg)
			m_sh = make_sp(new SubsetHandler(*m_mg));
		else
			m_sh = make_sp(new SubsetHandler());
	}
	return *m_sh;
}

template <class TElem, int dim>
const std::vector<int>* Partitioner_DynamicBisection<TElem, dim>::
get_process_map() const
{
	if(static_partitioning_enabled())
		return &m_procMap;
	else
		return NULL;
}

template <class TElem, int dim>
bool Partitioner_DynamicBisection<TElem, dim>::
partition(size_t baseLvl, size_t elementThreshold)
{
	GDIST_PROFILE_FUNC();
	typedef typename Grid::traits<elem_t>::iterator ElemIter;

	UG_COND_THROW(m_mg == NULL,
			"No grid was specified for Partitioner_DynamicBisection. "
			"partitioning can't be executed without a specified grid.");

	if(m_balanceWeights.invalid())
		m_balanceWeights = make_sp(new IBalanceWeights());

	MultiGrid& mg = *m_mg;
	if(m_sh.invalid())
		m_sh = make_sp(new SubsetHandler(mg));
	SubsetHandler& sh = *m_sh;
	sh.clear();
	//	int localProc = pcl::ProcRank();

	ANumber aWeight;
	mg.attach_to<elem_t>(aWeight);

//	assign all elements below baseLvl to the local process
	for(int i = 0; i < (int)baseLvl; ++i)
		sh.assign_subset(mg.begin<elem_t>(i), mg.end<elem_t>(i), 0);

	const ProcessHierarchy* procH;
	if(m_nextProcessHierarchy.valid())
		procH = m_nextProcessHierarchy.get();
	else
		procH = m_processHierarchy.get();

	m_procMap.clear();

	m_problemsOccurred = false;

//	iteprocHierarchy levels and perform rebalancing for all
//	hierarchy-sections which contain levels higher than baseLvl
	int oldHighestRedistLvl = m_highestRedistLevel;	// only used if static_partitioning is enabled
	for(size_t hlevel = 0; hlevel < procH->num_hierarchy_levels(); ++ hlevel)
	{
		int numProcs = procH->num_global_procs_involved(hlevel);

		int minLvl = procH->grid_base_level(hlevel);
		int maxLvl = (int)mg.top_level();

		if(m_balanceWeights->has_level_offsets()){
			++maxLvl;

			if(mg.top_level() < procH->grid_base_level(hlevel)){
			//	if the previous process-hierarchy had the same number of processes,
			//	we will silently ignore this hierarchy level. Only if it had
			//	a different amount of processes, m_problemsOccurred will be set
			//	to true, to indicate that another redistribution is necessary
			//	in order to distribute the grid over all involved processes.
				if((hlevel == 0) ||
					((int)procH->num_global_procs_involved(hlevel - 1) != numProcs))
				{
					UG_LOG("Partitioner_DynamicBisection: Ignoring hierarchy level "
						<< hlevel << " since it doesn't contain any elements yet\n");
					m_problemsOccurred = true;
				}
				continue;
			}
		}

		if(hlevel + 1 < procH->num_hierarchy_levels()){
		//	TODO:	currently elements which shall be considered in a level above
		//			and which do lie directly below a new hierarchy level would have to
		//			be partitioned separately
			if(m_balanceWeights->has_level_offsets() &&
				(mg.top_level() < procH->grid_base_level(hlevel+1)))
			{
				maxLvl = min<int>(maxLvl,
							(int)procH->grid_base_level(hlevel + 1));
			}
			else{
				maxLvl = min<int>(maxLvl,
							(int)procH->grid_base_level(hlevel + 1) - 1);
			}

		//	this code always starts partitioning one level above the highest
		//	hierarchy level to make sure that marked elements of the
		//	highest hierarchy level are considered
			// maxLvl = min<int>(maxLvl,
			// 			(int)procH->grid_base_level(hlevel + 1));
		}

		if(minLvl < (int)baseLvl)
			minLvl = (int)baseLvl;

		if(maxLvl < minLvl)
			continue;

		int maxValidLvl = min<int>(maxLvl, mg.top_level());

		int numPartitions = numProcs;
		if(static_partitioning_enabled())
			numPartitions = (int)procH->cluster_procs(hlevel).size();

		if((numProcs <= 1) || (numPartitions == 1)){
			for(int i = minLvl; i <= maxValidLvl; ++i)
				sh.assign_subset(mg.begin<elem_t>(i), mg.end<elem_t>(i), 0);
			continue;
		}

		if(static_partitioning_enabled() && ((int)hlevel <= m_highestRedistLevel)){
			for(int i = minLvl; i <= maxValidLvl; ++i)
				sh.assign_subset( mg.begin<elem_t>(i), mg.end<elem_t>(i), 0);
			continue;
		}

	//	if clustered siblings are enabled, we'll perform partitioning on the level
	//	below minLvl (if such a level exists). However, only the partition-map
	//	of minLvl and levels above will be adjusted.
		int partitionLvl = minLvl;
		pcl::ProcessCommunicator com(pcl::PCD_LOCAL);

		if((minLvl > 0) && base_class::clustered_siblings_enabled()){
			partitionLvl = minLvl - 1;
			size_t partitionHLvl = m_processHierarchy->hierarchy_level_from_grid_level(partitionLvl);
			if(!static_partitioning_enabled())
				com = m_processHierarchy->global_proc_com(partitionHLvl);
		}
		else if(!static_partitioning_enabled()){
			com = procH->global_proc_com(hlevel);
		}

		perform_bisection_new(numPartitions, minLvl, maxLvl, partitionLvl, aWeight, com);

		for(int i = minLvl; i < maxValidLvl; ++i){
			copy_partitions_to_children(sh, i);
		}

		if(static_partitioning_enabled()){
			m_procMap = procH->cluster_procs(hlevel);
			m_highestRedistLevel = hlevel;
		}
	}

//	make sure that everybody knows about the highestRedistLevel!
	pcl::ProcessCommunicator globCom;
	m_highestRedistLevel = globCom.allreduce(m_highestRedistLevel, PCL_RO_MAX);

	if(m_nextProcessHierarchy.valid()){
		*m_processHierarchy = *m_nextProcessHierarchy;
		m_nextProcessHierarchy = SPProcessHierarchy(NULL);
	}

	mg.detach_from<elem_t>(aWeight);

	if(static_partitioning_enabled()){
		if(m_procMap.empty() && (sh.num_subsets() > 0)){
			if(sh.num_subsets() != 1){
				UG_THROW("Something went wrong during partitioning. At this point"
						" either exactly one subset or a filled process map should exist.");
			}
			m_procMap.push_back(pcl::ProcRank());
		}
	}

//	debugging
//	static int execCounter = 0;
//	stringstream ss;
//	ss << "partition-map-" << execCounter << "-p" << pcl::ProcRank() << ".ugx";
//	AssignSubsetColors(sh);
//	SaveGridHierarchyTransformed(mg, sh, ss.str().c_str(), 20);
//	++execCounter;

	if(static_partitioning_enabled())
		return m_highestRedistLevel != oldHighestRedistLvl;
	else
		return true;
}


template <class TElem, int dim>
void Partitioner_DynamicBisection<TElem, dim>::
copy_partitions_to_children(ISubsetHandler& partitionSH, int lvl)
{
	GDIST_PROFILE_FUNC();
	typedef typename Grid::traits<elem_t>::iterator ElemIter;
	MultiGrid& mg = *m_mg;

//	assign partitions to all children in this hierarchy level
	for(ElemIter iter = mg.begin<elem_t>(lvl); iter != mg.end<elem_t>(lvl); ++iter)
	{
		size_t numChildren = mg.num_children<elem_t>(*iter);
		int si = partitionSH.get_subset_index(*iter);
		for(size_t i = 0; i < numChildren; ++i)
			partitionSH.assign_subset(mg.get_child<elem_t>(*iter, i), si);
	}

	if(mg.is_parallel()){
		GridLayoutMap& glm = mg.distributed_grid_manager()->grid_layout_map();
	//	communicate partitions from v-masters to v-slaves, since v-slaves
	//	havn't got no parents on their procs.
		ComPol_Subset<layout_t>	compolSHCopy(partitionSH, true);
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


template <class TElem, int dim>
void Partitioner_DynamicBisection<TElem, dim>::
gather_weights_from_level(int baseLvl, int childLvl, ANumber aWeight,
							   bool copyToVMastersOnBaseLvl)
{
	GDIST_PROFILE_FUNC();
	UG_DLOG(LIB_GRID, 1, "Partitioner_DynamicBisection-start gather_weights_from_level\n");
	typedef typename Grid::traits<elem_t>::iterator ElemIter;

	assert(m_mg);
	assert(m_mg->is_parallel());
	assert(m_mg->has_attachment<elem_t>(aWeight));

	IBalanceWeights& bw = *m_balanceWeights;
	MultiGrid& mg = *m_mg;
	Grid::AttachmentAccessor<elem_t, ANumber> aaWeight(mg, aWeight);

//	childLvl may be at most one level above the highest level of the grid hierarchy.
//	This is allowed to consider top-level elements which are to be considered one level
//	above their actual level...
	if(childLvl < baseLvl || childLvl > (int)mg.num_levels()){
		SetAttachmentValues(aaWeight, mg.begin<elem_t>(baseLvl), mg.end<elem_t>(baseLvl), 0);
		return;
	}
	else{
		DistributedGridManager& dgm = *mg.distributed_grid_manager();
		GridLayoutMap& glm = mg.distributed_grid_manager()->grid_layout_map();
		ComPol_CopyAttachment<layout_t, ANumber> compolCopy(mg, aWeight);

		if(childLvl < (int)mg.num_levels()){
			for(ElemIter iter = mg.begin<elem_t>(childLvl);
				iter != mg.end<elem_t>(childLvl); ++iter)
			{
				elem_t* e = *iter;
				aaWeight[e] = bw.get_weight(*iter);
			}
		}

		for(int lvl = childLvl - 1; lvl >= baseLvl; --lvl){
		//	copy from v-slaves to vmasters
			if(lvl + 1 < (int)mg.num_levels()){
				if(glm.has_layout<elem_t>(INT_V_SLAVE))
					m_intfcCom.send_data(glm.get_layout<elem_t>(INT_V_SLAVE).layout_on_level(lvl + 1),
										 compolCopy);
				if(glm.has_layout<elem_t>(INT_V_MASTER))
					m_intfcCom.receive_data(glm.get_layout<elem_t>(INT_V_MASTER).layout_on_level(lvl + 1),
											compolCopy);
				m_intfcCom.communicate();
			}

		//	accumulate child counts in parent elements on lvl
			if(bw.has_level_offsets()){
				for(ElemIter iter = mg.begin<elem_t>(lvl);
					iter != mg.end<elem_t>(lvl); ++iter)
				{
					elem_t* e = *iter;
					size_t numChildren = mg.num_children<elem_t>(e);
					if((lvl == childLvl - 1) && (numChildren == 0) && (!dgm.is_ghost(e))
						&& (bw.consider_in_level_above(e)))
					{
						aaWeight[e] = bw.get_refined_weight(e);
					}
					else{
					//	gather weights from children
						aaWeight[e] = 0;
						for(size_t i = 0; i < numChildren; ++i)
							aaWeight[e] += aaWeight[mg.get_child<elem_t>(e, i)];
					}
				}
			}
			else{
				for(ElemIter iter = mg.begin<elem_t>(lvl);
					iter != mg.end<elem_t>(lvl); ++iter)
				{
					elem_t* e = *iter;
					aaWeight[e] = 0;
					size_t numChildren = mg.num_children<elem_t>(e);
					for(size_t i = 0; i < numChildren; ++i)
						aaWeight[e] += aaWeight[mg.get_child<elem_t>(e, i)];
				}
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
	}

	UG_DLOG(LIB_GRID, 1, "Partitioner_DynamicBisection-stop gather_weights_from_level\n");
}


template <class TElem, int dim>
int Partitioner_DynamicBisection<TElem, dim>::
classify_elem(elem_t* e, int splitDim, number splitValue)
{
	typename elem_t::ConstVertexArray vrts = e->vertices();
	const size_t numVrts = e->num_vertices();

//	count the number of elements which are completely on the left or on the
//	right of the splitting plane
	int location = UNCLASSIFIED;
	for(size_t i = 0; i < numVrts; ++i){
		number vrtVal = m_aaPos[vrts[i]][splitDim];
		if(vrtVal < splitValue)
			location |= LEFT;
		else if(vrtVal > splitValue)
			location |= RIGHT;
		else
			location |= CUTTING;
	}

	return location;
}


template <class TElem, int dim>
void Partitioner_DynamicBisection<TElem, dim>::
perform_bisection_new(int numTargetProcs, int minLvl, int maxLvl, int partitionLvl,
				  ANumber aWeight, pcl::ProcessCommunicator com)
{
	// UG_LOG("Performing bisection for " << numTargetProcs << " procs on level "
	// 		<< partitionLvl << " for levels " << minLvl << " to " << maxLvl << "\n");

	SubsetHandler& sh = *m_sh;

	GDIST_PROFILE_FUNC();
	typedef typename MultiGrid::traits<elem_t>::iterator iter_t;

//	clear subset-assignment on partitionLvl
	MultiGrid& mg = *m_mg;
	DistributedGridManager* pdgm = mg.distributed_grid_manager();

	Grid::AttachmentAccessor<elem_t, ANumber> aaWeight(mg, aWeight);
	SetAttachmentValues(aaWeight, mg.begin<elem_t>(partitionLvl),
						mg.end<elem_t>(partitionLvl), 0);

//	arrays on which we'll perform partitioning
//	vector<elem_t*>	elems;
	m_entries.reserve(mg.num<elem_t>(partitionLvl));
	vector<TreeNode> treeNodes(1);

	vector<int> origSubsetIndices;
	if(partitionLvl < minLvl){
		origSubsetIndices.reserve(mg.num<elem_t>(partitionLvl));
		for(iter_t eiter = mg.begin<elem_t>(partitionLvl);
			eiter != mg.end<elem_t>(partitionLvl); ++eiter)
		{
			origSubsetIndices.push_back(sh.get_subset_index(*eiter));
		}
	}

//	invalidate target partitions of all elements in partitionLvl
	sh.assign_subset(mg.begin<elem_t>(partitionLvl),
					   mg.end<elem_t>(partitionLvl), -1);

//	iterate over all levels and gather child-counts in the partitionLvl
//TODO:	if partitionLvl == minLvl, we currently ignore level-offsets of some elements.
	for(int i_lvl = maxLvl; i_lvl >= minLvl; --i_lvl){
		gather_weights_from_level(partitionLvl, i_lvl, aWeight, false);

	//	now collect elements on partitionLvl, which have children in i_lvl but
	//	have not yet been partitioned.
		m_entries.clear();
		TreeNode& root = treeNodes.front();
		root.elems.set_entry_list(&m_entries);
		ElemList& elems = root.elems;
		elems.clear();

		number maxChildWeight = 0;// per element on partitionLvl
		for(iter_t eiter = mg.begin<elem_t>(partitionLvl);
			eiter != mg.end<elem_t>(partitionLvl); ++eiter)
		{
			elem_t* elem = *eiter;
			if((aaWeight[elem] > 0) && (sh.get_subset_index(elem) == -1)
				&& ((!pdgm) || (!pdgm->is_ghost(elem))))
			{
				m_entries.push_back(Entry(elem));
				elems.add(m_entries.size() - 1);
				maxChildWeight = max<number>(maxChildWeight, aaWeight[elem]);
			}
		}

		if(maxChildWeight == 0)
			maxChildWeight = 1;

		if(!com.empty()){
			maxChildWeight = com.allreduce(maxChildWeight, PCL_RO_MAX);

			root.firstProc = 0;
			root.numTargetProcs = numTargetProcs;

			control_bisection(sh, treeNodes, aWeight, maxChildWeight, com);
		}
	}

	if(partitionLvl < minLvl){
		UG_ASSERT(partitionLvl == minLvl - 1,
				  "partitionLvl and minLvl should be neighbors");

	//	copy subset indices from partition-level to minLvl
		for(int i = partitionLvl; i < minLvl; ++i){
			copy_partitions_to_children(sh, i);
		}

	//	reset partitions in the specified partition-level
		size_t counter = 0;
		for(iter_t eiter = mg.begin<elem_t>(partitionLvl);
			eiter != mg.end<elem_t>(partitionLvl); ++eiter, ++counter)
		{
			sh.assign_subset(*eiter, origSubsetIndices[counter]);
		}
	}
	else if(pdgm){
	//	copy subset indices from vertical slaves to vertical masters,
	//	since partitioning was only performed on vslaves
		GridLayoutMap& glm = pdgm->grid_layout_map();
		ComPol_Subset<layout_t>	compolSHCopy(sh, true);

		if(glm.has_layout<elem_t>(INT_V_SLAVE))
			m_intfcCom.send_data(glm.get_layout<elem_t>(INT_V_SLAVE).layout_on_level(partitionLvl),
								 compolSHCopy);
		if(glm.has_layout<elem_t>(INT_V_MASTER))
			m_intfcCom.receive_data(glm.get_layout<elem_t>(INT_V_MASTER).layout_on_level(partitionLvl),
									compolSHCopy);
		m_intfcCom.communicate();
	}
}


template <class TElem, int dim>
void Partitioner_DynamicBisection<TElem, dim>::
control_bisection(ISubsetHandler& partitionSH, std::vector<TreeNode>& treeNodes,
				  ANumber aWeight, number maxChildWeight, pcl::ProcessCommunicator& com)
{
//	for(size_t i = 0; i < treeNodes.size(); ++i){
//		UG_LOG("tree-node " << i << ": targetProc " << treeNodes[i].numTargetProcs
//				<< ", firstProc " << treeNodes[i].firstProc << "\n");
//	}

//	UG_LOG("\nCONTROL-BISECTION-STARTS\n");
//	at the beginning of this function, each node in treeNodes has to have
//	valid 'elem', 'firstProc', 'numTargetProcs' and members.
//	numTargetProcs has to be at least 2.
	GDIST_PROFILE_FUNC();

//	calculate the number of child-nodes required and the split-ratio for each parent node
	vector<TreeNode> childNodes;
	childNodes.reserve(treeNodes.size() * 2);
	for(size_t iNode = 0; iNode < treeNodes.size(); ++iNode){
		TreeNode& tn = treeNodes[iNode];
		tn.bisectionComplete = false;

		if(tn.numTargetProcs < 2){
			tn.bisectionComplete = true;
			tn.firstChildNode = s_invalidIndex;
			tn.ratioLeft = 0.5;
		//	the tree node is a leaf node and we can thus assign all of its elements
		//	to the final partition
			size_t numElems = 0;
			for(size_t i = tn.elems.first(); i != s_invalidIndex; i = tn.elems.next(i), ++numElems)
				partitionSH.assign_subset(tn.elems.elem(i), tn.firstProc);

//			UG_LOG("assigned " << numElems << " elements to proc " << tn.firstProc << "\n");

			continue;
		}

		int numTargetProcsLeft = tn.numTargetProcs / 2;
		int numTargetProcsRight = tn.numTargetProcs - numTargetProcsLeft;
		tn.ratioLeft = (number)numTargetProcsLeft / (number)tn.numTargetProcs;

		tn.firstChildNode = childNodes.size();
		childNodes.resize(childNodes.size() + 2);
		TreeNode& childNodeLeft = childNodes[tn.firstChildNode];
		TreeNode& childNodeRight = childNodes[tn.firstChildNode + 1];

		childNodeLeft.elems.set_entry_list(&m_entries);
		childNodeLeft.firstProc = tn.firstProc;
		childNodeLeft.numTargetProcs = numTargetProcsLeft;

		childNodeRight.elems.set_entry_list(&m_entries);
		childNodeRight.firstProc = tn.firstProc + numTargetProcsLeft;
		childNodeRight.numTargetProcs = numTargetProcsRight;

	}

//	if there are no child nodes, we're done.
	if(childNodes.empty())
		return;

//	perform the actual bisection
	bisect_elements(childNodes, treeNodes, aWeight, maxChildWeight, com, 0);

//	perform the recursion
	control_bisection(partitionSH, childNodes, aWeight, maxChildWeight, com);

}


template <class TElem, int dim>
void Partitioner_DynamicBisection<TElem, dim>::
calculate_global_dimensions(vector<TreeNode>& treeNodes,
							number maxChildWeight, ANumber aWeight,
							pcl::ProcessCommunicator& com)
{
	GDIST_PROFILE_FUNC();
	MultiGrid& mg = *m_mg;
	Grid::AttachmentAccessor<elem_t, ANumber> aaWeight(mg, aWeight);

	vector<double> weightBuf(treeNodes.size());
	vector<double> centerBuf(treeNodes.size() * dim);
	vector<double> boxMinBuf(treeNodes.size() * dim);
	vector<double> boxMaxBuf(treeNodes.size() * dim);

	for(size_t iNode = 0; iNode < treeNodes.size(); ++iNode){
		TreeNode& tn = treeNodes[iNode];
		tn.totalWeight = 0;
		VecSet(tn.center, 0);
		VecSet(tn.boxMin, numeric_limits<double>::max());
		VecSet(tn.boxMax, -numeric_limits<double>::max());

		if(tn.firstChildNode != s_invalidIndex){
			ElemList& elems = tn.elems;
			for(size_t iElem = elems.first(); iElem != s_invalidIndex;
				iElem = elems.next(iElem))
			{
				elem_t* e = elems.elem(iElem);
				vector_t p = CalculateCenter(e, m_aaPos);
				number w = aaWeight[e] / maxChildWeight;
				VecScaleAdd(tn.center, 1, tn.center, w, p);
				tn.totalWeight += w;

				for(int i_dim = 0; i_dim < dim; ++i_dim){
					if(p[i_dim] < tn.boxMin[i_dim])
						tn.boxMin[i_dim] = p[i_dim];
					if(p[i_dim] > tn.boxMax[i_dim])
						tn.boxMax[i_dim] = p[i_dim];
				}
			}
		}

		weightBuf[iNode] = tn.totalWeight;
		for(int d = 0; d < dim; ++d){
			int ind = iNode * dim + d;
			centerBuf[ind] = tn.center[d];
			boxMinBuf[ind] = tn.boxMin[d];
			boxMaxBuf[ind] = tn.boxMax[d];
		}
	}


	vector<double> gWeightBuf, gCenterBuf, gBoxMinBuf, gBoxMaxBuf;
	com.allreduce(weightBuf, gWeightBuf, PCL_RO_SUM);
	com.allreduce(centerBuf, gCenterBuf, PCL_RO_SUM);
	com.allreduce(boxMinBuf, gBoxMinBuf, PCL_RO_MIN);
	com.allreduce(boxMaxBuf, gBoxMaxBuf, PCL_RO_MAX);


	for(size_t iNode = 0; iNode < treeNodes.size(); ++iNode){
		TreeNode& tn = treeNodes[iNode];
		for(int d = 0; d < dim; ++d){
			int ind = iNode * dim + d;
			tn.center[d] = gCenterBuf[ind] / gWeightBuf[iNode];
			tn.boxMin[d] = gBoxMinBuf[ind];
			tn.boxMax[d] = gBoxMaxBuf[ind];
		}
	}
}

template <class TElem, int dim>
void Partitioner_DynamicBisection<TElem, dim>::
improve_split_values(vector<TreeNode>& treeNodes,
					 size_t maxIterations, ANumber aWeight,
					 pcl::ProcessCommunicator& com)
{
	vector<bool> gotSplitValue(treeNodes.size(), false);
	for(size_t iNode = 0; iNode < treeNodes.size(); ++iNode){
		if(treeNodes[iNode].firstChildNode == s_invalidIndex)
			gotSplitValue[iNode] = true;
	}

	Grid::AttachmentAccessor<elem_t, ANumber> aaWeight(*m_mg, aWeight);
	vector<double> weights, gWeights;
	weights.reserve(treeNodes.size() * NUM_CONSTANTS);

	for(size_t iteration = 0; iteration < maxIterations; ++iteration){
		weights.assign(treeNodes.size() * NUM_CONSTANTS, 0);

		for(size_t iNode = 0; iNode < treeNodes.size(); ++iNode){
			if(gotSplitValue[iNode])
				continue;

			TreeNode& tn = treeNodes[iNode];
//			UG_LOG(">node " << iNode << " sv: " << tn.splitValue
//					<< ", svmin: " << tn.minSplitValue << ", svmax: " << tn.maxSplitValue << "\n");

			size_t firstWgt = iNode * NUM_CONSTANTS;
			ElemList& elems = tn.elems;

			for(size_t i = elems.first(); i != s_invalidIndex; i = elems.next(i)){
				elem_t* e = elems.elem(i);
				int location = classify_elem(e, tn.splitDim, tn.splitValue);
				weights[firstWgt + location] += aaWeight[e];
				weights[firstWgt + TOTAL] += aaWeight[e];
				if(location == CUTTING){
					if(CalculateCenter(e, m_aaPos)[tn.splitDim] < tn.splitValue)
						weights[firstWgt + CUTTING_CENTER_LEFT] += aaWeight[e];
					else
						weights[firstWgt + CUTTING_CENTER_RIGHT] += aaWeight[e];
				}
			}
		}

		com.allreduce(weights, gWeights, PCL_RO_SUM);

	//		UG_LOG("weights unclassified: " << gWeights[UNCLASSIFIED] << endl);
	//		UG_LOG("weights left: " << gWeights[LEFT] << endl);
	//		UG_LOG("weights right: " << gWeights[RIGHT] << endl);
	//		UG_LOG("weights cutting: " << gWeights[CUTTING] << endl);
	//		UG_LOG("weights cutting-center-left: " << gWeights[CUTTING_CENTER_LEFT] << endl);
	//		UG_LOG("weights cutting-center-right: " << gWeights[CUTTING_CENTER_RIGHT] << endl);
	//		UG_LOG("weights total: " << gWeights[TOTAL] << endl);

	//	check whether both sides are below the splitRatio
		for(size_t iNode = 0; iNode < treeNodes.size(); ++iNode){
//			UG_LOG("NODE " << iNode << "\n");
			if(gotSplitValue[iNode])
				continue;

//			UG_LOG("improving split value of node " << iNode << endl);
			TreeNode& tn = treeNodes[iNode];
//			UG_LOG("Node ratio left: " << tn.ratioLeft << endl);

			size_t firstWgt = iNode * NUM_CONSTANTS;
			if(gWeights[firstWgt + TOTAL] > 0){
				bool leftOk = (gWeights[firstWgt + LEFT] / gWeights[firstWgt + TOTAL] <= tn.ratioLeft);
				bool rightOk = (gWeights[firstWgt + RIGHT] / gWeights[firstWgt + TOTAL] <= (1. - tn.ratioLeft));

				if(!(leftOk && rightOk)){
				//	check if the ratio would be good if we also considered those with CENTER_LEFT and CENTER_RIGHT
					number ratioLeft = (gWeights[firstWgt + LEFT] + gWeights[firstWgt + CUTTING_CENTER_LEFT])
										/ (gWeights[firstWgt + TOTAL] * tn.ratioLeft);
					number ratioRight = (gWeights[firstWgt + RIGHT] + gWeights[firstWgt + CUTTING_CENTER_RIGHT])
										/ (gWeights[firstWgt + TOTAL] * (1. - tn.ratioLeft));

					if((ratioLeft > m_tolerance) && (ratioRight > m_tolerance)){
//						UG_LOG(" ratio-for-center-split: " << ratio << endl);
//						UG_LOG(" center split\n");
					//	the split-value is fine!
						gotSplitValue[iNode] = true;
						continue;
					}
				}

				if(!leftOk){
//					UG_LOG("left with ratio: " << (number)gWeights[firstWgt + LEFT] / (number)gWeights[firstWgt + TOTAL] << endl);
				//	move the split-value to the left
					tn.maxSplitValue = tn.splitValue;
					tn.splitValue = (tn.minSplitValue + tn.splitValue) / 2.;
				}
				else if(!rightOk){
//					UG_LOG("right with ratio: " << (number)gWeights[firstWgt + RIGHT] / (number)gWeights[firstWgt + TOTAL] << endl);
				//	move the split-value to the right
					tn.minSplitValue = tn.splitValue;
					tn.splitValue = (tn.splitValue + tn.maxSplitValue) / 2.;
				}
				else{
//					UG_LOG("Got good split value: " << tn.splitValue << "\n");
				//	the split value seems to be fine
					gotSplitValue[iNode] = true;
				}
			}
		}
	}
}

template <class TElem, int dim>
void Partitioner_DynamicBisection<TElem, dim>::
bisect_elements(vector<TreeNode>& childNodesOut, vector<TreeNode>& parentNodes,
				ANumber aWeight, number maxChildWeight, pcl::ProcessCommunicator& com,
				int cutRecursion)
{
	GDIST_PROFILE_FUNC();
	MultiGrid& mg = *m_mg;

//	calculate bounding box and center for the current tree-nodes
	calculate_global_dimensions(parentNodes, maxChildWeight, aWeight, com);

//	we will now calculate the split-dim and split-value for each tree-node
	for(size_t iNode = 0; iNode < parentNodes.size(); ++iNode){
		TreeNode& tn = parentNodes[iNode];
		if(!tn.bisectionComplete){
			tn.splitDim = 0;
			for(int i = 1; i < dim; ++i){
				if((tn.boxMax[i] - tn.boxMin[i]) > (tn.boxMax[tn.splitDim] - tn.boxMin[tn.splitDim]))
					tn.splitDim = i;
			}

		//	this is an initial guess
			tn.minSplitValue = tn.boxMin[tn.splitDim];
			tn.maxSplitValue = tn.boxMax[tn.splitDim];
			tn.splitValue = (1. - 2. * tn.ratioLeft) * tn.minSplitValue
							+ 2. * tn.ratioLeft * tn.center[tn.splitDim];

//			UG_LOG("node " << iNode << ":\n");
//			UG_LOG("  center: " << tn.center << ", boxMin: " << tn.boxMin << ", boxMax: " << tn.boxMax << endl);
//			UG_LOG("  splitDim: " << tn.splitDim << ", splitValue: " << tn.splitValue << endl);
		}
	}

	improve_split_values(parentNodes, m_splitImproveIterations, aWeight, com);


	Grid::AttachmentAccessor<elem_t, ANumber> aaWeight(mg, aWeight);

	vector<ElemList> elemsCutVec(parentNodes.size());
	for(size_t i = 0; i < elemsCutVec.size(); ++i){
		elemsCutVec[i].set_entry_list(parentNodes[i].elems.entries());
	}

	if(cutRecursion < dim - 1){
		const size_t numWgtsPerNode = 5;
		vector<double> weights(parentNodes.size() * numWgtsPerNode, 0);

		for(size_t iNode = 0; iNode < parentNodes.size(); ++iNode){
			TreeNode& tn = parentNodes[iNode];
			if(!tn.bisectionComplete){
				size_t firstWgtInd = iNode * numWgtsPerNode;
				ElemList& elems = tn.elems;
				ElemList& elemsCut = elemsCutVec[iNode];
				ElemList& elemsLeftOut = childNodesOut[tn.firstChildNode].elems;
				ElemList& elemsRightOut = childNodesOut[tn.firstChildNode + 1].elems;
				for(size_t i = elems.first(); i != s_invalidIndex;){
					elem_t* e = elems.elem(i);
					size_t iNext = elems.next(i);

					int loc = classify_elem(e, tn.splitDim, tn.splitValue);
					switch(loc){
						case LEFT:
							elemsLeftOut.add(i);
							weights[firstWgtInd] += aaWeight[e];
							break;
						case RIGHT:
							elemsRightOut.add(i);
							weights[firstWgtInd + 1] += aaWeight[e];
							break;
						case CUTTING:
							elemsCut.add(i);
							weights[firstWgtInd + 2] += aaWeight[e];
						//	check whether the center is left or right, since we can
						//	eventually avoid additional bisection
							if(CalculateCenter(e, m_aaPos)[tn.splitDim] < tn.splitValue)
								weights[firstWgtInd + 3] += aaWeight[e];
							else
								weights[firstWgtInd + 4] += aaWeight[e];
							break;
						default:
							UG_THROW("INVALID CLASSIFICATION DURING BISECTION\n");
							break;
					}

					i = iNext;
				}
			}

		//	since all elements from the elem list have been assigned to other arrays,
		//	we have to clear elems, since it would be invalid anyways
			tn.elems.clear();
		}

		vector<double> gWeights;
		com.allreduce(weights, gWeights, PCL_RO_SUM);

		for(size_t iNode = 0; iNode < parentNodes.size(); ++iNode){
			TreeNode& tn = parentNodes[iNode];
			if(!tn.bisectionComplete){
//				UG_LOG("Performing bisection with ratio: " << tn.ratioLeft << endl);

				size_t firstWgtInd = iNode * numWgtsPerNode;
				ElemList& elemsCut = elemsCutVec[iNode];
				ElemList& elemsLeftOut = childNodesOut[tn.firstChildNode].elems;
				ElemList& elemsRightOut = childNodesOut[tn.firstChildNode + 1].elems;

			//	calculate the number of elements which have to go to the left side
				double gWTotal = gWeights[firstWgtInd] + gWeights[firstWgtInd + 1] + gWeights[firstWgtInd + 2];
				double gMissingLeft = tn.ratioLeft * gWTotal - gWeights[firstWgtInd];
				double gMissingRight = (1. - tn.ratioLeft) * gWTotal - gWeights[firstWgtInd + 1];

		//		UG_LOG("weights left:  " << gWeights[0] << endl);
		//		UG_LOG("weights right: " << gWeights[1] << endl);
		//		UG_LOG("weights cut:   " << gWeights[2] << endl);

			//	bisect the cutting elements
				if(gMissingLeft <= 0){
//					UG_LOG("Adding all to right\n");
				//	add all cut-elements to the right side
					for(size_t i = elemsCut.first(); i != s_invalidIndex;){
						size_t iNext = elemsCut.next(i);
						elemsRightOut.add(i);
						i = iNext;
					}
					elemsCut.clear();
					tn.bisectionComplete = true;
				}
				else if(gMissingRight <= 0){
//					UG_LOG("Adding all to left\n");
				//	add all cut-elements to the left side
					for(size_t i = elemsCut.first(); i != s_invalidIndex;){
						size_t iNext = elemsCut.next(i);
						elemsLeftOut.add(i);
						i = iNext;
					}

					elemsCut.clear();
					tn.bisectionComplete = true;
				}
				else{
				//	if the ratios would be fine if one simply bisected by midpoints, we'll do that
//					number wLeft = gWeights[0] + gWeights[3];
//					number wRight = gWeights[1] + gWeights[4];
//					number ratio;
//					if(wLeft < wRight)	ratio = wLeft / wRight;
//					else				ratio = wRight / wLeft;
//
//					if(ratio > 0.99){
					number ratioLeft = (gWeights[0] + gWeights[3]) / (gWTotal * tn.ratioLeft);
					number ratioRight = (gWeights[1] + gWeights[4]) / (gWTotal * (1. - tn.ratioLeft));

					if((ratioLeft > m_tolerance) && (ratioRight > m_tolerance)){
//						UG_LOG("direct cut bisection on node " << iNode
//								<< " with ratioLeft = " << ratioLeft
//								<< " and ratioRight = " << ratioRight << "\n");

						for(size_t i = elemsCut.first(); i != s_invalidIndex;){
							elem_t* e = elemsCut.elem(i);
							size_t iNext = elemsCut.next(i);
							if(CalculateCenter(e, m_aaPos)[tn.splitDim] < tn.splitValue)
								elemsLeftOut.add(i);
							else
								elemsRightOut.add(i);
							i = iNext;
						}
						elemsCut.clear();
						tn.bisectionComplete = true;
					}
				}
			//	prepare for recursion
				if(!tn.bisectionComplete){
//					UG_LOG("Preparing for recursion in node " << iNode << endl);
//					UG_LOG("elemsLeftOut.size(): " << elemsLeftOut.size() <<
//							", elemsRightOut.size(): " << elemsRightOut.size() << "\n");
					number newRatioLeft = (number)gMissingLeft / (number)(gMissingLeft + gMissingRight);
					tn.ratioLeft = newRatioLeft;
					tn.elems = elemsCut;
				}
			}
		}

		bisect_elements(childNodesOut, parentNodes, aWeight,
						maxChildWeight, com, cutRecursion + 1);
	}
	else{
	//	perform a simple bisection
		for(size_t iNode = 0; iNode < parentNodes.size(); ++iNode){
			TreeNode& tn = parentNodes[iNode];
			if(!tn.bisectionComplete){
//				UG_LOG("performing simple bisection on node " << iNode
//						<< " with ratio " << tn.ratioLeft << "\n");
				ElemList& elems = tn.elems;
				ElemList& elemsLeftOut = childNodesOut[tn.firstChildNode].elems;
				ElemList& elemsRightOut = childNodesOut[tn.firstChildNode + 1].elems;
				for(size_t i = elems.first(); i != s_invalidIndex;){
					size_t iNext = elems.next(i);
					elem_t* e = elems.elem(i);
					if(CalculateCenter(e, m_aaPos)[tn.splitDim] <= tn.splitValue)
						elemsLeftOut.add(i);
					else
						elemsRightOut.add(i);
					i = iNext;
				}
				elems.clear();
			}
		}
	}
}

template <class TElem, int dim>
void Partitioner_DynamicBisection<TElem, dim>::
enable_static_partitioning(bool enable)
{
	m_staticPartitioning = true;
}

template <class TElem, int dim>
bool Partitioner_DynamicBisection<TElem, dim>::
static_partitioning_enabled() const
{
	return m_staticPartitioning;
}


template class Partitioner_DynamicBisection<Edge, 1>;
template class Partitioner_DynamicBisection<Edge, 2>;
template class Partitioner_DynamicBisection<Face, 2>;
template class Partitioner_DynamicBisection<Edge, 3>;
template class Partitioner_DynamicBisection<Face, 3>;
template class Partitioner_DynamicBisection<Volume, 3>;

}// end of namespace
