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
	m_mg(NULL)
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
}

template<int dim>
void Partitioner_DynamicBisection<dim>::
set_next_process_hierarchy(SPProcessHierarchy procHierarchy)
{
	m_nextProcessHierarchy = procHierarchy;
}

template<int dim>
void Partitioner_DynamicBisection<dim>::
set_balance_weights(SmartPtr<BalanceWeights<dim> >)
{
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
void Partitioner_DynamicBisection<dim>::
partition(size_t baseLvl, size_t elementThreshold)
{
	UG_LOG("Dynamic Bisection Partitioning begins...\n");

	typedef typename Grid::traits<elem_t>::iterator ElemIter;

	assert(m_mg);
	MultiGrid& mg = *m_mg;
	m_sh.clear();
	int localProc = pcl::GetProcRank();

	AInt aInt;
	mg.attach_to<elem_t>(aInt);

//	assign all elements below baseLvl to the local process
	for(int i = 0; i < (int)baseLvl; ++i)
		m_sh.assign_subset(mg.begin<elem_t>(i), mg.end<elem_t>(i), localProc);

	const ProcessHierarchy* procH;
	if(m_nextProcessHierarchy.valid())
		procH = m_nextProcessHierarchy.get();
	else
		procH = m_processHierarchy.get();

//	iteprocHierarchy levels and perform rebalancing for all
//	hierarchy-sections which contain levels higher than baseLvl
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

	//	if clustered siblings are enabled, we'll perform partitioning on the level
	//	below minLvl (if such a level exists). However, only the partition-map
	//	of minLvl and levels above will be adjusted.
		int partitionLvl = minLvl;
		pcl::ProcessCommunicator com = procH->global_proc_com(hlevel);
		if((minLvl > 0) && base_class::clustered_siblings_enabled()){
			partitionLvl = minLvl - 1;
			size_t partitionHLvl = m_processHierarchy->hierarchy_level_from_grid_level(partitionLvl);
			com = m_processHierarchy->global_proc_com(partitionHLvl);
		}

		perform_bisection(numProcs, minLvl, maxLvl, partitionLvl, aInt, com);

		for(int i = minLvl; i < maxLvl; ++i){
			copy_partitions_to_children(m_sh, i);
		}
	}

//	make sure that everybody knows about the highestRedistLevel!
	if(m_nextProcessHierarchy.valid()){
		*m_processHierarchy = *m_nextProcessHierarchy;
		m_nextProcessHierarchy = SPProcessHierarchy(NULL);
	}

	mg.detach_from<elem_t>(aInt);
	UG_LOG("Dynamic Bisection Partitioning ends...\n");

//	debugging
	static int execCounter = 0;
	stringstream ss;
	ss << "partition-map-" << execCounter << "-p" << pcl::GetProcRank() << ".ugx";
	AssignSubsetColors(m_sh);
	SaveGridHierarchyTransformed(mg, m_sh, ss.str().c_str(), 0.1);
	++execCounter;
}


template <int dim>
void Partitioner_DynamicBisection<dim>::
copy_partitions_to_children(ISubsetHandler& partitionSH, int lvl)
{
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
				  AInt aChildCount, pcl::ProcessCommunicator com)
{
	typedef typename MultiGrid::traits<elem_t>::iterator iter_t;

//	clear subset-assignment on partitionLvl
	MultiGrid& mg = *m_mg;
	DistributedGridManager* pdgm = mg.distributed_grid_manager();

	Grid::AttachmentAccessor<elem_t, AInt> aaChildCount(mg, aChildCount);
	SetAttachmentValues(aaChildCount, mg.begin<elem_t>(partitionLvl),
						mg.end<elem_t>(partitionLvl), 0);

//	arrays on which we'll perform partitioning
	vector<elem_t*>	elems;
	elems.reserve(mg.num<elem_t>(partitionLvl));

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
		gather_child_counts_from_level(partitionLvl, i_lvl, aChildCount, false);

	//	now collect elements on partitionLvl, which have children in i_lvl but
	//	have not yet been partitioned.
		elems.clear();
		int maxNumChildren = 0;// per element on partitionLvl
		for(iter_t eiter = mg.begin<elem_t>(partitionLvl);
			eiter != mg.end<elem_t>(partitionLvl); ++eiter)
		{
			elem_t* elem = *eiter;
			if((aaChildCount[elem] > 0) && (m_sh.get_subset_index(elem) == -1)
				&& ((!pdgm) || (!pdgm->is_ghost(elem))))
			{
				elems.push_back(elem);
				maxNumChildren = max(maxNumChildren, aaChildCount[elem]);
			}
		}

		if(maxNumChildren == 0)
			maxNumChildren = 1;

		if(!com.empty()){
			maxNumChildren = com.allreduce(maxNumChildren, PCL_RO_MAX);
			bisect_elements(m_sh, elems, maxNumChildren, numTargetProcs,
							0, aChildCount, com);
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
bisect_elements(ISubsetHandler& partitionSH, const vector<elem_t*>& elems,
				int maxNumChildren, int numTargetProcs, int firstProc,
				AInt aChildCount, pcl::ProcessCommunicator& com)
{
	MultiGrid& mg = *m_mg;

	Grid::AttachmentAccessor<elem_t, AInt> aaChildCount(mg, aChildCount);
	Grid::VertexAttachmentAccessor<apos_t> aaPos(mg, m_aPos);

	vector_t gCenter, gBoxMin, gBoxMax;
	calculate_global_dimensions(gCenter, gBoxMin, gBoxMax, elems,
								maxNumChildren, aChildCount, com);

//todo	Adjust center if numTargetProcs can't be divided by 2
//		this can e.g. be done by performing some iterations of the following algo:
//		- After an initial splitValue and minSplit and maxSplit values have been chosen:
//		- count elements left and right of splitValue.
//		- adjust splitValue, minSplit and maxSplit values using interval-bisection
//		- repeat with step 2

//	split the array into 2 sub-arrays along the global-box's largest dimension
	int splitDim = 0;
	for(int i = 1; i < dim; ++i){
		if((gBoxMax[i] - gBoxMin[i]) > (gBoxMax[splitDim] - gBoxMin[splitDim]))
			splitDim = i;
	}

	number splitValue = gCenter[splitDim];

	if(numTargetProcs == 2){
		for(size_t i = 0; i < elems.size(); ++i){
			elem_t* e = elems[i];
			if(CalculateCenter(e, aaPos)[splitDim] <= splitValue){
				//partitionSH.assign_subset(e, procMap[firstProc]);
				partitionSH.assign_subset(e, firstProc);
			}
			else{
				//partitionSH.assign_subset(e, procMap[firstProc + 1]);
				partitionSH.assign_subset(e, firstProc + 1);
			}
		}
	}
	else if(numTargetProcs == 3){
		vector<elem_t*>	elemsRight;
		elemsRight.reserve(elems.size() * 0.55);
		for(size_t i = 0; i < elems.size(); ++i){
			elem_t* e = elems[i];
			if(CalculateCenter(e, aaPos)[splitDim] <= splitValue){
				//partitionSH.assign_subset(e, procMap[firstProc]);
				partitionSH.assign_subset(e, firstProc);
			}
			else{
				elemsRight.push_back(e);
			}
		}
		bisect_elements(partitionSH, elemsRight, maxNumChildren, 2, firstProc + 1,
						aChildCount, com);
	}
	else if(numTargetProcs > 3){
		vector<elem_t*>	elemsLeft, elemsRight;
		elemsLeft.reserve(elems.size() * 0.55);
		elemsRight.reserve(elems.size() * 0.55);

		for(size_t i = 0; i < elems.size(); ++i){
			elem_t* e = elems[i];
			if(CalculateCenter(e, aaPos)[splitDim] <= splitValue){
				elemsLeft.push_back(e);
			}
			else{
				elemsRight.push_back(e);
			}
		}

		int numTargetProcsLeft = numTargetProcs / 2;
		int numTargetProcsRight = numTargetProcs - numTargetProcsLeft;

		bisect_elements(partitionSH, elemsLeft, maxNumChildren, numTargetProcsLeft,
						firstProc, aChildCount, com);
		bisect_elements(partitionSH, elemsRight, maxNumChildren, numTargetProcsRight,
						firstProc + numTargetProcsLeft, aChildCount, com);
	}
	else{
		UG_THROW("At least 2 target processes have to be specified for a bisection");
	}
}


template <int dim>
void Partitioner_DynamicBisection<dim>::
calculate_global_dimensions(vector_t& centerOut, vector_t& boxMinOut,
							vector_t& boxMaxOut, const vector<elem_t*>& elems,
							int maxNumChildren, AInt aChildCount,
							pcl::ProcessCommunicator& com)
{
	MultiGrid& mg = *m_mg;
	Grid::AttachmentAccessor<elem_t, AInt> aaChildCount(mg, aChildCount);
	Grid::VertexAttachmentAccessor<apos_t> aaPos(mg, m_aPos);

	number totalWeight = 0;
	vector_t center;
	VecSet(center, 0);

//	calculate weighted center
	vector_t minCorner, maxCorner;
	VecSet(minCorner, 0);
	VecSet(maxCorner, 0);

	if(!elems.empty())
		minCorner = maxCorner = CalculateCenter(elems.front(), aaPos);

	for(size_t i = 0; i < elems.size(); ++i){
		elem_t* e = elems[i];
		vector_t p = CalculateCenter(e, aaPos);
		number w = (number)aaChildCount[e] / maxNumChildren;
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
//	UG_LOG("local center: ");
//	if(totalWeight > 0 && elems.size() > 0){
//		vector_t tmpCenter;
//		VecScale(tmpCenter, center, 1./(totalWeight * (number)elems.size()));
//		UG_LOG(tmpCenter << endl);
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
	if(com.get_local_proc_id() == 0){
		for(size_t i = 0; i < recvBuf.size(); i += comDataSize){
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

//	broadcast the calculated center and bounding box to all involved processes.
	vector_t broadcastData[3];
	broadcastData[0] = gCenter;
	broadcastData[1] = gMinCorner;
	broadcastData[2] = gMaxCorner;

	com.broadcast((char*)broadcastData, sizeof(vector_t) * 3, 0);

	centerOut = broadcastData[0];
	boxMinOut = broadcastData[1];
	boxMaxOut = broadcastData[2];
}

template <int dim>
void Partitioner_DynamicBisection<dim>::
gather_child_counts_from_level(int baseLvl, int childLvl, AInt aInt,
							   bool copyToVMastersOnBaseLvl)
{
	GDIST_PROFILE_FUNC();
	UG_DLOG(LIB_GRID, 1, "Partitioner_DynamicBisection-start gather_child_counts_from_level\n");
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
	UG_DLOG(LIB_GRID, 1, "Partitioner_DynamicBisection-stop gather_child_counts_from_level\n");
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
	return NULL;
}

template class Partitioner_DynamicBisection<1>;
template class Partitioner_DynamicBisection<2>;
template class Partitioner_DynamicBisection<3>;

}// end of namespace
