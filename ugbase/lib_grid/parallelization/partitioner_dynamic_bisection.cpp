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

template <int dim>
Partitioner_DynamicBisection<dim>::
Partitioner_DynamicBisection() :
	m_mg(NULL),
	m_staticPartitioning(false),
	m_highestRedistLevel(-1)
{
	m_processHierarchy = SPProcessHierarchy(new ProcessHierarchy);
	m_processHierarchy->add_hierarchy_level(0, 1);
}

template<int dim>
Partitioner_DynamicBisection<dim>::
~Partitioner_DynamicBisection()
{
}

template<int dim>
void Partitioner_DynamicBisection<dim>::
set_grid(MultiGrid* mg, Attachment<MathVector<dim> > aPos)
{
	m_mg = mg;
	m_sh.assign_grid(m_mg);
	m_aPos = aPos;
	m_aaPos.access(*m_mg, m_aPos);
}

template<int dim>
void Partitioner_DynamicBisection<dim>::
set_next_process_hierarchy(SPProcessHierarchy procHierarchy)
{
	m_nextProcessHierarchy = procHierarchy;
}

template<int dim>
void Partitioner_DynamicBisection<dim>::
set_balance_weights(SmartPtr<BalanceWeights<dim> > balanceWeights)
{
	m_balanceWeights = balanceWeights;
}

template<int dim>
void Partitioner_DynamicBisection<dim>::
set_connection_weights(SmartPtr<ConnectionWeights<dim> >)
{
}


template<int dim>
ConstSPProcessHierarchy Partitioner_DynamicBisection<dim>::
current_process_hierarchy() const
{
	return m_processHierarchy;
}


template<int dim>
ConstSPProcessHierarchy Partitioner_DynamicBisection<dim>::
next_process_hierarchy() const
{
	return m_nextProcessHierarchy;
}


template<int dim>
bool Partitioner_DynamicBisection<dim>::
supports_balance_weights() const
{
	return false;
}

template<int dim>
bool Partitioner_DynamicBisection<dim>::
supports_connection_weights() const
{
	return false;
}

template<int dim>
number Partitioner_DynamicBisection<dim>::
estimate_distribution_quality(std::vector<number>* pLvlQualitiesOut)
{
//todo	Consider connection weights in the final quality!
	typedef typename Grid::traits<elem_t>::iterator ElemIter;
	using std::min;

	MultiGrid& mg = *m_mg;
	DistributedGridManager& distGridMgr = *mg.distributed_grid_manager();

	number minQuality = 1;

	if(pLvlQualitiesOut)
		pLvlQualitiesOut->clear();

//	calculate the quality estimate.
//todo The quality of a level could be weighted by the total amount of elements
//		in each level.
	const ProcessHierarchy* procH;
	if(m_nextProcessHierarchy.valid())
		procH = m_nextProcessHierarchy.get();
	else
		procH = m_processHierarchy.get();

	for(size_t lvl = 0; lvl < mg.num_levels(); ++lvl){
		size_t hlvl = procH->hierarchy_level_from_grid_level(lvl);
		int numProcs = procH->num_global_procs_involved(hlvl);
		if(numProcs <= 1){
			if(pLvlQualitiesOut)
				pLvlQualitiesOut->push_back(1.0);
			continue;
		}

		pcl::ProcessCommunicator procComAll = procH->global_proc_com(hlvl);
		if(!procComAll.empty()){
			int localWeight = 0;
			for(ElemIter iter = mg.begin<elem_t>(lvl);
				iter != mg.end<elem_t>(lvl); ++iter)
			{
				if(!distGridMgr.is_ghost(*iter))
					localWeight += 1;
			}

			int maxWeight = procComAll.allreduce(localWeight, PCL_RO_MAX);
			int minWeight = procComAll.allreduce(localWeight, PCL_RO_MIN);
			//int totalWeight = procComAll.allreduce(localWeight, PCL_RO_SUM);
			//number averageWeight = totalWeight / procComAll.size();

			if(maxWeight > 0){
				number quality = (number)minWeight / (number)maxWeight;
				minQuality = min(minQuality, quality);
				if(pLvlQualitiesOut)
					pLvlQualitiesOut->push_back(quality);
			}
			else if(pLvlQualitiesOut)
				pLvlQualitiesOut->push_back(-1);
		}
		else if(pLvlQualitiesOut)
			pLvlQualitiesOut->push_back(-1);
	}

	pcl::ProcessCommunicator comGlobal;
	return comGlobal.allreduce(minQuality, PCL_RO_MIN);
}

template<int dim>
bool Partitioner_DynamicBisection<dim>::
partition(size_t baseLvl, size_t elementThreshold)
{
	GDIST_PROFILE_FUNC();
	typedef typename Grid::traits<elem_t>::iterator ElemIter;

	assert(m_mg);
	MultiGrid& mg = *m_mg;
	m_sh.clear();
//	int localProc = pcl::GetProcRank();

	ANumber aWeight;
	mg.attach_to<elem_t>(aWeight);

//	assign all elements below baseLvl to the local process
	for(int i = 0; i < (int)baseLvl; ++i)
		m_sh.assign_subset(mg.begin<elem_t>(i), mg.end<elem_t>(i), 0);

	const ProcessHierarchy* procH;
	if(m_nextProcessHierarchy.valid())
		procH = m_nextProcessHierarchy.get();
	else
		procH = m_processHierarchy.get();

//	iteprocHierarchy levels and perform rebalancing for all
//	hierarchy-sections which contain levels higher than baseLvl
	int oldHighestRedistLvl = m_highestRedistLevel;	// only used if static_partitioning is enabled
	for(size_t hlevel = 0; hlevel < procH->num_hierarchy_levels(); ++ hlevel)
	{
		int minLvl = procH->grid_base_level(hlevel);
		int maxLvl = (int)mg.top_level();
		if(hlevel + 1 < procH->num_hierarchy_levels()){
			maxLvl = min<int>(maxLvl,
						(int)procH->grid_base_level(hlevel + 1) - 1);
		}

		if(minLvl < (int)baseLvl)
			minLvl = (int)baseLvl;

		if(maxLvl < minLvl)
			continue;

		int numProcs = procH->num_global_procs_involved(hlevel);

		if(numProcs <= 1){
			for(int i = minLvl; i <= maxLvl; ++i)
				m_sh.assign_subset(mg.begin<elem_t>(i), mg.end<elem_t>(i), 0);
			continue;
		}


		if(static_partitioning_enabled() && ((int)hlevel <= m_highestRedistLevel)){
			for(int i = minLvl; i <= maxLvl; ++i)
				m_sh.assign_subset( mg.begin<elem_t>(i), mg.end<elem_t>(i), 0);
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

		int numPartitions = numProcs;
		if(static_partitioning_enabled())
			numPartitions = (int)procH->cluster_procs(hlevel).size();

		perform_bisection(numPartitions, minLvl, maxLvl, partitionLvl, aWeight, com);

		for(int i = minLvl; i < maxLvl; ++i){
			copy_partitions_to_children(m_sh, i);
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

//	debugging
//	static int execCounter = 0;
//	stringstream ss;
//	ss << "partition-map-" << execCounter << "-p" << pcl::GetProcRank() << ".ugx";
//	AssignSubsetColors(m_sh);
//	SaveGridHierarchyTransformed(mg, m_sh, ss.str().c_str(), 20);
//	++execCounter;

	if(static_partitioning_enabled())
		return m_highestRedistLevel != oldHighestRedistLvl;
	else
		return true;
}


template <int dim>
void Partitioner_DynamicBisection<dim>::
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


template <int dim>
void Partitioner_DynamicBisection<dim>::
perform_bisection(int numTargetProcs, int minLvl, int maxLvl, int partitionLvl,
				  ANumber aWeight, pcl::ProcessCommunicator com)
{
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
	ElemList elems(m_entries);

	vector<int> origSubsetIndices;
	if(partitionLvl < minLvl){
		origSubsetIndices.reserve(mg.num<elem_t>(partitionLvl));
		for(iter_t eiter = mg.begin<elem_t>(partitionLvl);
			eiter != mg.end<elem_t>(partitionLvl); ++eiter)
		{
			origSubsetIndices.push_back(m_sh.get_subset_index(*eiter));
		}
	}

//	invalidate target partitions of all elements in partitionLvl
	m_sh.assign_subset(mg.begin<elem_t>(partitionLvl),
					   mg.end<elem_t>(partitionLvl), -1);

//	iterate over all levels and gather child-counts in the partitionLvl
	for(int i_lvl = maxLvl; i_lvl >= minLvl; --i_lvl){
		gather_weights_from_level(partitionLvl, i_lvl, aWeight, false);

	//	now collect elements on partitionLvl, which have children in i_lvl but
	//	have not yet been partitioned.
		m_entries.clear();
		elems.clear();
		number maxChildWeight = 0;// per element on partitionLvl
		for(iter_t eiter = mg.begin<elem_t>(partitionLvl);
			eiter != mg.end<elem_t>(partitionLvl); ++eiter)
		{
			elem_t* elem = *eiter;
			if((aaWeight[elem] > 0) && (m_sh.get_subset_index(elem) == -1)
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
			control_bisection(m_sh, elems, maxChildWeight, numTargetProcs,
							  0, aWeight, com);
		}
	}

	if(partitionLvl < minLvl){
		UG_ASSERT(partitionLvl = minLvl - 1,
				  "partitionLvl and minLvl should be neighbors");

	//	copy subset indices from partition-level to minLvl
		for(int i = partitionLvl; i < minLvl; ++i){
			copy_partitions_to_children(m_sh, i);
		}

	//	reset partitions in the specified partition-level
		size_t counter = 0;
		for(iter_t eiter = mg.begin<elem_t>(partitionLvl);
			eiter != mg.end<elem_t>(partitionLvl); ++eiter, ++counter)
		{
			m_sh.assign_subset(*eiter, origSubsetIndices[counter]);
		}
	}
	else if(pdgm){
	//	copy subset indices from vertical slaves to vertical masters,
	//	since partitioning was only performed on vslaves
		GridLayoutMap& glm = pdgm->grid_layout_map();
		ComPol_Subset<layout_t>	compolSHCopy(m_sh, true);

		if(glm.has_layout<elem_t>(INT_V_SLAVE))
			m_intfcCom.send_data(glm.get_layout<elem_t>(INT_V_SLAVE).layout_on_level(partitionLvl),
								 compolSHCopy);
		if(glm.has_layout<elem_t>(INT_V_MASTER))
			m_intfcCom.receive_data(glm.get_layout<elem_t>(INT_V_MASTER).layout_on_level(partitionLvl),
									compolSHCopy);
		m_intfcCom.communicate();
	}
}


template <int dim>
void Partitioner_DynamicBisection<dim>::
control_bisection(ISubsetHandler& partitionSH, ElemList& elems,
				  number maxChildWeight, int numTargetProcs, int firstProc,
				  ANumber aWeight, pcl::ProcessCommunicator& com)
{
	GDIST_PROFILE_FUNC();

	int numTargetProcsLeft = numTargetProcs / 2;
	int numTargetProcsRight = numTargetProcs - numTargetProcsLeft;
	number ratioLeft = (number)numTargetProcsLeft / (number)numTargetProcs;

	ElemList elemsLeft(elems.entries()), elemsRight(elems.entries());

	size_t numElemsIn = elems.size();
	bisect_elements(elemsLeft, elemsRight, elems, ratioLeft, aWeight,
					maxChildWeight, com, 0);

	UG_COND_THROW(numElemsIn != elemsLeft.size() + elemsRight.size(),
				  "Not all elements have been assigned to one of the two sides: "
				  << "in: " << numElemsIn << ", left: " << elemsLeft.size()
				  << ", right: " << elemsRight.size());

	if(numTargetProcs == 2){
		for(size_t i = elemsLeft.first(); i != s_invalidIndex; i = elemsLeft.next(i))
			partitionSH.assign_subset(elemsLeft.elem(i), firstProc);
		for(size_t i = elemsRight.first(); i != s_invalidIndex; i = elemsRight.next(i))
			partitionSH.assign_subset(elemsRight.elem(i), firstProc + 1);
	}
	else if(numTargetProcs == 3){
		for(size_t i = elemsLeft.first(); i != s_invalidIndex; i = elemsLeft.next(i))
			partitionSH.assign_subset(elemsLeft.elem(i), firstProc);
		control_bisection(partitionSH, elemsRight, maxChildWeight, 2, firstProc + 1,
						  aWeight, com);
	}
	else if(numTargetProcs > 3){
		control_bisection(partitionSH, elemsLeft, maxChildWeight, numTargetProcsLeft,
						  firstProc, aWeight, com);
		control_bisection(partitionSH, elemsRight, maxChildWeight, numTargetProcsRight,
						  firstProc + numTargetProcsLeft, aWeight, com);
	}
	else{
		UG_THROW("At least 2 target processes have to be specified for a bisection");
	}
}


template<int dim>
void Partitioner_DynamicBisection<dim>::
bisect_elements(ElemList& elemsLeftOut, ElemList& elemsRightOut,
		   	    ElemList& elems, number ratioLeft, ANumber aWeight,
		   	    number maxChildWeight, pcl::ProcessCommunicator& com, int cutRecursion)
{
//	UG_LOG("Performing bisection with ratio: " << ratioLeft << endl);
	GDIST_PROFILE_FUNC();
	MultiGrid& mg = *m_mg;

	Grid::AttachmentAccessor<elem_t, ANumber> aaWeight(mg, aWeight);

	vector_t gCenter, gBoxMin, gBoxMax;
	calculate_global_dimensions(gCenter, gBoxMin, gBoxMax, elems,
								maxChildWeight, aWeight, com);

//	split the array into 2 sub-arrays along the global-box's largest dimension
	int splitDim = 0;
	for(int i = 1; i < dim; ++i){
		if((gBoxMax[i] - gBoxMin[i]) > (gBoxMax[splitDim] - gBoxMin[splitDim]))
			splitDim = i;
	}

	//UG_LOG("BBox: " << gBoxMin << ", " << gBoxMax << ", Center: " << gCenter << endl);
	number splitValue = find_split_value(elems, splitDim, ratioLeft, gCenter[splitDim],
										 gBoxMin[splitDim], gBoxMax[splitDim],
										 10, aWeight, com);

//	UG_LOG("splitDim: " << splitDim << ", splitValue: " << splitValue << endl);

	if(cutRecursion < dim - 1){
		ElemList elemsCut(elems.entries());

		double weights[5] = {0, 0, 0, 0, 0};

		for(size_t i = elems.first(); i != s_invalidIndex;){
			elem_t* e = elems.elem(i);
			size_t iNext = elems.next(i);

			int loc = classify_elem(e, splitDim, splitValue);
			switch(loc){
				case LEFT:
					elemsLeftOut.add(i);
					weights[0] += aaWeight[e];
					break;
				case RIGHT:
					elemsRightOut.add(i);
					weights[1] += aaWeight[e];
					break;
				case CUTTING:
					elemsCut.add(i);
					weights[2] += aaWeight[e];
				//	check whether the center is left or right, since we can
				//	eventually avoid additional bisection
					if(CalculateCenter(e, m_aaPos)[splitDim] < splitValue)
						weights[3] += aaWeight[e];
					else
						weights[4] += aaWeight[e];
					break;
				default:
					UG_THROW("INVALID CLASSIFICATION DURING BISECTION\n");
					break;
			}

			i = iNext;
		}

	//	since all elements from the elems list have been assigned to other arrays,
	//	we have to clear elems, since it would be invalid anyways
		elems.clear();

		double gWeights[5];
		com.allreduce(weights, gWeights, 5, PCL_RO_SUM);

		//	calculate the number of elements which have to go to the left side
		double gWTotal = gWeights[0] + gWeights[1] + gWeights[2];
		double gMissingLeft = ratioLeft * gWTotal - gWeights[0];
		double gMissingRight = (1. - ratioLeft) * gWTotal - gWeights[1];

//		UG_LOG("weights left:  " << gWeights[0] << endl);
//		UG_LOG("weights right: " << gWeights[1] << endl);
//		UG_LOG("weights cut:   " << gWeights[2] << endl);

	//	bisect the cutting elements
		if(gMissingLeft <= 0){
//			UG_LOG("Adding all to right\n");
		//	add all cut-elements to the right side
			for(size_t i = elemsCut.first(); i != s_invalidIndex;){
				size_t iNext = elemsCut.next(i);
				elemsRightOut.add(i);
				i = iNext;
			}
			elemsCut.clear();
			return;
		}

		if(gMissingRight <= 0){
//			UG_LOG("Adding all to left\n");
		//	add all cut-elements to the left side
			for(size_t i = elemsCut.first(); i != s_invalidIndex;){
				size_t iNext = elemsCut.next(i);
				elemsLeftOut.add(i);
				i = iNext;
			}

			elemsCut.clear();
			return;
		}

	//	if the ratios would be fine if one simply bisected by midpoints, we'll do that
		number wLeft = gWeights[0] + gWeights[3];
		number wRight = gWeights[1] + gWeights[4];
		number ratio;
		if(wLeft < wRight)	ratio = wLeft / wRight;
		else				ratio = wRight / wLeft;

		if(ratio > 0.99){
//			UG_LOG("direct cut bisection\n");
			for(size_t i = elemsCut.first(); i != s_invalidIndex;){
				elem_t* e = elemsCut.elem(i);
				size_t iNext = elemsCut.next(i);
				if(CalculateCenter(e, m_aaPos)[splitDim] < splitValue)
					elemsLeftOut.add(i);
				else
					elemsRightOut.add(i);
				i = iNext;
			}
			elemsCut.clear();
			return;
		}

	//	perform a recursion
		number newRatioLeft = (number)gMissingLeft / (number)(gMissingLeft + gMissingRight);
		bisect_elements(elemsLeftOut, elemsRightOut, elemsCut, newRatioLeft, aWeight,
						maxChildWeight, com, cutRecursion + 1);
	}
	else{
	//	perform a simple bisection
		for(size_t i = elems.first(); i != s_invalidIndex;){
			size_t iNext = elems.next(i);
			elem_t* e = elems.elem(i);
			if(CalculateCenter(e, m_aaPos)[splitDim] <= splitValue)
				elemsLeftOut.add(i);
			else
				elemsRightOut.add(i);
			i = iNext;
		}
		elems.clear();
	}
}


template <int dim>
void Partitioner_DynamicBisection<dim>::
calculate_global_dimensions(vector_t& centerOut, vector_t& boxMinOut,
							vector_t& boxMaxOut, const ElemList& elems,
							number maxChildWeight, ANumber aWeight,
							pcl::ProcessCommunicator& com)
{
	GDIST_PROFILE_FUNC();
	MultiGrid& mg = *m_mg;
	Grid::AttachmentAccessor<elem_t, ANumber> aaWeight(mg, aWeight);

	number totalWeight = 0;
	vector_t center;
	VecSet(center, 0);

//	calculate weighted center
	vector_t minCorner, maxCorner;
	VecSet(minCorner, 0);
	VecSet(maxCorner, 0);

	if(!elems.empty())
		minCorner = maxCorner = CalculateCenter(elems.elem(elems.first()), m_aaPos);

	for(size_t i = elems.first(); i != s_invalidIndex; i = elems.next(i)){
		elem_t* e = elems.elem(i);
		vector_t p = CalculateCenter(e, m_aaPos);
		number w = aaWeight[e] / maxChildWeight;
		VecScaleAdd(center, 1, center, w, p);
		totalWeight += w;

		for(int i_dim = 0; i_dim < dim; ++i_dim){
			if(p[i_dim] < minCorner[i_dim])
				minCorner[i_dim] = p[i_dim];
			else if(p[i_dim] > maxCorner[i_dim])
				maxCorner[i_dim] = p[i_dim];
		}
	}

//	UG_LOG("num elems: " << elems.size() << endl);
//	UG_LOG("local BBox: " << minCorner << ", " << maxCorner << endl);
//	if(totalWeight > 0 && elems.size() > 0){
//		vector_t tmpCenter;
//		VecScale(tmpCenter, center, 1./ totalWeight);
//		UG_LOG("local center: " << tmpCenter << endl);
//	}
//	else{
//		UG_LOG("---\n");
//	}

//	communicate data to proc 0
	static const int ofsNumElems = 0;
	static const int ofsCenter = 1;
	static const int ofsBoxMin = ofsCenter + dim;
	static const int ofsBoxMax = ofsBoxMin + dim;
	static const int ofsWeight = ofsBoxMax + dim;
	static const int comDataSize = ofsWeight + 1;

	double localComData[comDataSize];
	localComData[ofsNumElems] = elems.size();
	for(int i = 0; i < dim; ++i){
		localComData[i + ofsCenter] = center[i];
		localComData[i + ofsBoxMin] = minCorner[i];
		localComData[i + ofsBoxMax] = maxCorner[i];
	}
	localComData[ofsWeight] = totalWeight;

	vector<double> recvBuf(com.size() * comDataSize);
	com.gather(localComData, comDataSize, PCL_DT_DOUBLE, &recvBuf.front(),
			   comDataSize, PCL_DT_DOUBLE, 0);

//	calculate global center and global bounding box on proc 0
	vector_t gCenter;
	VecSet(gCenter, 0);
	vector_t gMinCorner = minCorner;
	vector_t gMaxCorner = maxCorner;
	bool validBox = false;

	number totalNumElems = 0;
	totalWeight = 0;
	//UG_LOG("com.get_local_proc_id(): " << com.get_local_proc_id() << endl);
	if(com.get_local_proc_id() == 0){
		for(size_t i = 0; i < recvBuf.size(); i += comDataSize){
			//UG_LOG("from proc " << i/comDataSize << ":\n");
			//UG_LOG("  #elems: " << recvBuf[i + ofsNumElems] << endl);
		//	ignore contributions of empty processes
			if(recvBuf[i + ofsNumElems] == 0)
				continue;

			totalNumElems += recvBuf[i + ofsNumElems];
			vector_t c;
			for(int i_dim = 0; i_dim < dim; ++i_dim)
				c[i_dim] = recvBuf[i + ofsCenter + i_dim];

			if(validBox){
				for(int i_dim = 0; i_dim < dim; ++i_dim){
					if(recvBuf[i + ofsBoxMin + i_dim] < gMinCorner[i_dim])
						gMinCorner[i_dim] = recvBuf[i + ofsBoxMin + i_dim];
					else if(recvBuf[i + ofsBoxMax + i_dim] > gMaxCorner[i_dim])
						gMaxCorner[i_dim] = recvBuf[i + ofsBoxMax + i_dim];
				}
			}
			else{
				for(int i_dim = 0; i_dim < dim; ++i_dim){
					gMinCorner[i_dim] = recvBuf[i + ofsBoxMin + i_dim];
					gMaxCorner[i_dim] = recvBuf[i + ofsBoxMax + i_dim];
				}
				validBox = true;
			}
			number w = recvBuf[i + ofsWeight];
			totalWeight += w;
			VecAdd(gCenter, gCenter, c);
		}
	}

	if(totalWeight > 0)
		VecScale(gCenter, gCenter, 1./totalWeight);

//	UG_LOG("global center pre-com: " << gCenter << endl);
//	broadcast the calculated center and bounding box to all involved processes.
	vector_t broadcastData[3];
	broadcastData[0] = gCenter;
	broadcastData[1] = gMinCorner;
	broadcastData[2] = gMaxCorner;

	com.broadcast((char*)broadcastData, sizeof(vector_t) * 3, 0);

	centerOut = broadcastData[0];
	boxMinOut = broadcastData[1];
	boxMaxOut = broadcastData[2];
//	UG_LOG("global center post-com: " << centerOut << endl);

}

template <int dim>
void Partitioner_DynamicBisection<dim>::
gather_weights_from_level(int baseLvl, int childLvl, ANumber aWeight,
							   bool copyToVMastersOnBaseLvl)
{
	GDIST_PROFILE_FUNC();
	UG_DLOG(LIB_GRID, 1, "Partitioner_DynamicBisection-start gather_weights_from_level\n");
	typedef typename Grid::traits<elem_t>::iterator ElemIter;

	assert(m_mg);
	assert(m_mg->is_parallel());
	assert(m_mg->has_attachment<elem_t>(aWeight));

	MultiGrid& mg = *m_mg;
	Grid::AttachmentAccessor<elem_t, ANumber> aaWeight(mg, aWeight);


	if(childLvl < baseLvl){
		SetAttachmentValues(aaWeight, mg.begin<elem_t>(baseLvl), mg.end<elem_t>(baseLvl), 0);
		return;
	}
	else if(childLvl == baseLvl){
		if(m_balanceWeights.valid()){
			BalanceWeights<dim>& bw = *m_balanceWeights;
			for(ElemIter iter = mg.begin<elem_t>(baseLvl);
				iter != mg.end<elem_t>(baseLvl); ++iter)
			{
				elem_t* e = *iter;
				aaWeight[e] = bw.get_weight(*iter);
			}
		}
		else
			SetAttachmentValues(aaWeight, mg.begin<elem_t>(baseLvl), mg.end<elem_t>(baseLvl), 1);
		return;
	}
	else{
		GridLayoutMap& glm = mg.distributed_grid_manager()->grid_layout_map();
		ComPol_CopyAttachment<layout_t, ANumber> compolCopy(mg, aWeight);

	//	initialize weights.
	//	we don't fill the weights on childLvl but on the level below.
		if(m_balanceWeights.valid()){
		//	we have to write a balance-weight into all elements of childLvl
			BalanceWeights<dim>& bw = *m_balanceWeights;
			for(ElemIter iter = mg.begin<elem_t>(childLvl - 1);
				iter != mg.end<elem_t>(childLvl - 1); ++iter)
			{
				elem_t* parent = *iter;
				aaWeight[parent] = 0;
				size_t numChildren = mg.num_children<elem_t>(parent);
				for(size_t ichild = 0; ichild < numChildren; ++ichild)
				{
					aaWeight[parent] += bw.get_weight(mg.get_child<elem_t>(parent, ichild));
				}
			}
		}
		else{
			for(ElemIter iter = mg.begin<elem_t>(childLvl - 1);
				iter != mg.end<elem_t>(childLvl - 1); ++iter)
			{
				elem_t* e = *iter;
				aaWeight[e] = mg.num_children<elem_t>(e);
			}
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
				aaWeight[e] = 0;
				size_t numChildren = mg.num_children<elem_t>(e);
				for(size_t i = 0; i < numChildren; ++i)
					aaWeight[e] += aaWeight[mg.get_child<elem_t>(e, i)];
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


template<int dim>
SubsetHandler& Partitioner_DynamicBisection<dim>::
get_partitions()
{
	return m_sh;
}

template<int dim>
const std::vector<int>* Partitioner_DynamicBisection<dim>::
get_process_map() const
{
	if(static_partitioning_enabled())
		return &m_procMap;
	else
		return NULL;
}


template<int dim>
int Partitioner_DynamicBisection<dim>::
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


template<int dim>
number Partitioner_DynamicBisection<dim>::
find_split_value(const ElemList& elems, int splitDim,
				 number splitRatio, number initialGuess,
				 number minValue, number maxValue,
				 size_t maxIterations, ANumber aWeight,
				 pcl::ProcessCommunicator& com)
{
	//UG_LOG("splitRatio: " << splitRatio << ", initialGuess: " << initialGuess << endl);
	if(splitRatio < SMALL)
		return minValue + SMALL;
	if(splitRatio > 1. - SMALL)
		return maxValue - SMALL;

	Grid::AttachmentAccessor<elem_t, ANumber> aaWeight(*m_mg, aWeight);

	number splitValue = initialGuess;
	for(size_t iteration = 0; iteration < maxIterations; ++iteration){
		int numElems[NUM_CONSTANTS];
		for(int i = 0; i < NUM_CONSTANTS; ++i)
			numElems[i] = 0;

		for(size_t i = elems.first(); i != s_invalidIndex; i = elems.next(i)){
			elem_t* e = elems.elem(i);
			int location = classify_elem(e, splitDim, splitValue);
			numElems[location] += aaWeight[e];
			numElems[TOTAL] += aaWeight[e];
		}

		int gNumElems[NUM_CONSTANTS];
		com.allreduce(numElems, gNumElems, NUM_CONSTANTS, PCL_RO_SUM);

	//	check whether both sides are below the splitRatio
		if(gNumElems[TOTAL] > 0){
			bool leftOk = ((number)gNumElems[LEFT] / (number)gNumElems[TOTAL] <= splitRatio);
			bool rightOk = ((number)gNumElems[RIGHT] / (number)gNumElems[TOTAL] <= (1. - splitRatio));

			if(!leftOk){
				//UG_LOG("left with ratio: " << (number)gNumElems[LEFT] / (number)gNumElems[TOTAL] << endl);
			//	move the split-value to the left
				maxValue = splitValue;
				splitValue = (minValue + splitValue) / 2.;
			}
			else if(!rightOk){
				//UG_LOG("right with ratio: " << (number)gNumElems[RIGHT] / (number)gNumElems[TOTAL] << endl);
			//	move the split-value to the right
				minValue = splitValue;
				splitValue = (splitValue + maxValue) / 2.;
			}
			else{
				//UG_LOG("Found split-value after " << iteration + 1 << " iterations\n");
				return splitValue;
			}
		}
	}

	//UG_LOG("No perfect split-value found!\n");
	return splitValue;
}

template<int dim>
void Partitioner_DynamicBisection<dim>::
enable_static_partitioning(bool enable)
{
	m_staticPartitioning = true;
}

template<int dim>
bool Partitioner_DynamicBisection<dim>::
static_partitioning_enabled() const
{
	return m_staticPartitioning;
}


template class Partitioner_DynamicBisection<1>;
template class Partitioner_DynamicBisection<2>;
template class Partitioner_DynamicBisection<3>;

}// end of namespace
