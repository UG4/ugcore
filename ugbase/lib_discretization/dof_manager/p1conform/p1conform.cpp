#include "./p1conform.h"
#include <vector>
#include <queue>
#include <algorithm>


namespace ug{

///////////////////////////////////////
// P1StorageManager
///////////////////////////////////////

void
P1StorageManager::
set_subset_handler(ISubsetHandler& sh)
{
	if(m_pSH != NULL) clear();

	m_pSH = &sh;
	m_pSH->enable_subset_attachments(true);
}

void
P1StorageManager::
clear_subset_handler()
{
	if(m_pSH != NULL) clear();
	m_pSH->enable_subset_attachments(false);
	m_pSH = NULL;
}

void
P1StorageManager::
clear()
{
	if(m_pSH == NULL) return;

	for(size_t si = 0; si < m_vSubsetInfo.size(); ++si)
	{
		m_pSH->detach_from<VertexBase>(m_vSubsetInfo[si].aDoF, si);
		m_vSubsetInfo[si].aaDoFVRT.invalidate();
	}
	m_vSubsetInfo.clear();
}

void
P1StorageManager::
update_attachments()
{
	if(m_pSH == NULL) return;

	size_t num_subsets =  m_pSH->num_subsets();

	// Create level dof distributors
	for(size_t si = m_vSubsetInfo.size(); si < num_subsets; ++si)
	{
		m_vSubsetInfo.push_back(SubsetInfo());
		m_pSH->attach_to<VertexBase>(m_vSubsetInfo[si].aDoF, si);
		m_vSubsetInfo[si].aaDoFVRT.access(*m_pSH, m_vSubsetInfo[si].aDoF, si);
	}
}


///////////////////////////////////////
// P1ConformFunctionPattern
///////////////////////////////////////

bool
P1ConformFunctionPattern::
add_discrete_function(const char* name, LocalShapeFunctionSetID id, int dim)
{
	// for a P1 dof manager only Lagrange P1 function space is permitted
	if(id != LocalShapeFunctionSetID(LocalShapeFunctionSetID::LAGRANGE, 1))
	{
		UG_LOG("P1ConformDoFDistributor: Only LSFS_LAGRANGEP1 functions are supported.\n");
		return false;
	}

	return FunctionPattern::add_discrete_function(name, id, dim);
}

bool
P1ConformFunctionPattern::
add_discrete_function(const char* name, LocalShapeFunctionSetID id, const SubsetGroup& SubsetIndices, int dim)
{
	// for a P1 dof manager only Lagrange P1 function space is permitted
	if(id != LocalShapeFunctionSetID(LocalShapeFunctionSetID::LAGRANGE, 1))
	{
		UG_LOG("P1ConformDoFDistributor: Only LSFS_LAGRANGEP1 functions are supported.\n");
		return false;
	}

	return FunctionPattern::add_discrete_function(name, id, SubsetIndices, dim);
}


///////////////////////////////////////
// P1ConformDoFDistribution
///////////////////////////////////////

size_t
P1ConformDoFDistribution::
num_indices(ReferenceObjectID refID, int si, const FunctionGroup& funcGroup) const
{
	const ReferenceElement& refElem = ReferenceElementFactory::get_reference_element(refID);

	size_t numFct = 0;
	for(size_t fct = 0; fct < funcGroup.num_fct(); ++fct)
	{
		if(is_def_in_subset(funcGroup[fct], si))
			numFct++;
	}

	return numFct * refElem.num_obj(0);
}

bool
P1ConformDoFDistribution::
prepare_indices(ReferenceObjectID refID, int si, LocalIndices& ind, bool withHanging) const
{
	const ReferenceElement& refElem = ReferenceElementFactory::get_reference_element(refID);

	if(!withHanging)
	{
		ind.clear();
		size_t numInd = 0;
		size_t numFct = 0;
		for(size_t fct = 0; fct < ind.num_fct(); ++fct)
		{
			if(!is_def_in_subset(ind.fct_id(fct), si)) continue;
			for(size_t dof = 0; dof < refElem.num_obj(0); ++dof)
			{
				LocalIndices::multi_index_type dof_ind;
				dof_ind[0] = dof + numFct * refElem.num_obj(0);
				dof_ind[1] = 0;
				ind.add_dof(fct, dof_ind);
			}
			numInd += refElem.num_obj(0);
			numFct++;
		}
		ind.set_num_indices(numInd);
	}
	return true;
}


size_t
P1ConformDoFDistribution::
num_inner_indices(ReferenceObjectID refID, int si, const FunctionGroup& funcGroup) const
{
	if(refID != ROID_VERTEX) return 0;
	else return num_indices(refID, si, funcGroup);
}

bool
P1ConformDoFDistribution::
prepare_inner_indices(ReferenceObjectID refID, int si, LocalIndices& ind) const
{
	ind.clear();
	if(refID != ROID_VERTEX) return true;
	else return prepare_indices(refID, si, ind);
}

///////////// creation /////////////////

bool
P1ConformDoFDistribution::
distribute_dofs()
{
	if(m_pStorageManager == NULL)
	{
		UG_LOG("In 'P1ConformDoFDistribution::distribute_dofs:"
				"Storage Manager not set. Aborting.\n");
		return false;
	}

	if(m_pFunctionPattern == NULL)
	{
		UG_LOG("In 'P1ConformDoFDistribution::distribute_dofs:"
				"Function Pattern not set. Aborting.\n");
		return false;
	}

	// Attach indices
	m_pStorageManager->update_attachments();

	// iterators
	geometry_traits<VertexBase>::iterator iter, iterBegin, iterEnd;

	// loop subsets
	size_t i = 0;
	size_t i_per_subset = 0;

	// reset number of dofs
	m_vNumDoFs.clear(); m_vNumDoFs.resize(num_subsets(), 0);
	m_vvOffsets.clear(); m_vvOffsets.resize(num_subsets());

	// loop subsets
	for(int si = 0; si < num_subsets(); ++si)
	{
		// create offsets
		size_t count = 0;
		m_vvOffsets[si].resize(m_pFunctionPattern->num_fct());
		for(size_t fct = 0; fct < m_pFunctionPattern->num_fct(); ++fct)
		{
			if(!is_def_in_subset(fct, si)) m_vvOffsets[si][fct] = (size_t) -1;
			else m_vvOffsets[si][fct] = count++;
		}

		// skip if no dofs to be distributed
		if(!(m_pFunctionPattern->num_fct(si)>0)) continue;

		iterBegin = m_goc.begin<VertexBase>(si);
		iterEnd =  m_goc.end<VertexBase>(si);

		// remember number of functions
		size_t num_fct =  m_pFunctionPattern->num_fct(si);

		// loop Vertices
		i_per_subset = 0;
		for(iter = iterBegin; iter != iterEnd; ++iter)
		{
			// get vertex
			VertexBase* vrt = *iter;

			// write index
			m_pStorageManager->m_vSubsetInfo[si].aaDoFVRT[vrt] = i;
			i += num_fct;
			i_per_subset += num_fct;
		}
		m_vNumDoFs[si] = i_per_subset;
	}
	m_numDoFs = i;

//	order
	if(!order_cuthill_mckee())
	{
		UG_LOG("In 'P1ConformDoFDistribution::distribute_dofs':"
				" Error while ordering dofs.\n");
		return false;
	}

	return true;
}

///////////////////////////////////////
// GroupedP1ConformDoFDistribution
///////////////////////////////////////


size_t
GroupedP1ConformDoFDistribution::
num_indices(ReferenceObjectID refID, int si, const FunctionGroup& funcGroup) const
{
	const ReferenceElement& refElem = ReferenceElementFactory::get_reference_element(refID);

	for(size_t fct = 0; fct < funcGroup.num_fct(); ++fct)
	{
		if(is_def_in_subset(funcGroup[fct], si))
			return refElem.num_obj(0);
	}

	return 0;
}


size_t
GroupedP1ConformDoFDistribution::
num_inner_indices(ReferenceObjectID refID, int si, const FunctionGroup& funcGroup) const
{
	if(refID != ROID_VERTEX) return 0;
	else return num_indices(refID, si, funcGroup);
}

bool
GroupedP1ConformDoFDistribution::
prepare_indices(ReferenceObjectID refID, int si, LocalIndices& ind, bool withHanging) const
{
	const ReferenceElement& refElem = ReferenceElementFactory::get_reference_element(refID);

	ind.clear();
	size_t numInd = 0;
	for(size_t fct = 0; fct < ind.num_fct(); ++fct)
	{
		if(!is_def_in_subset(ind.fct_id(fct), si)) continue;

		for(size_t dof = 0; dof < refElem.num_obj(0); ++dof)
		{
			LocalIndices::multi_index_type dof_ind;
			dof_ind[0] = dof;
			dof_ind[1] = ind.fct_id(fct);
			ind.add_dof(fct, dof_ind);
		}
		numInd = refElem.num_obj(0);
	}
	ind.set_num_indices(numInd);
	return true;
}

bool
GroupedP1ConformDoFDistribution::
prepare_inner_indices(ReferenceObjectID refID, int si, LocalIndices& ind) const
{
	ind.clear();
	if(refID != ROID_VERTEX) return 0;
	else return prepare_indices(refID, si, ind);
}

///////// Creation /////////////

bool
GroupedP1ConformDoFDistribution::
distribute_dofs()
{
	if(m_pStorageManager == NULL)
	{
		UG_LOG("In 'P1ConformDoFDistribution::distribute_dofs:"
				"Storage Manager not set. Aborting.\n");
		return false;
	}

	if(m_pFunctionPattern == NULL)
	{
		UG_LOG("In 'P1ConformDoFDistribution::distribute_dofs:"
				"Function Pattern not set. Aborting.\n");
		return false;
	}

	// Attach indices
	m_pStorageManager->update_attachments();

	// iterators
	geometry_traits<VertexBase>::iterator iter, iterBegin, iterEnd;

	// loop subsets
	size_t i = 0;
	size_t i_per_subset = 0;

	// reset number of dofs
	m_vNumDoFs.clear(); m_vNumDoFs.resize(num_subsets(), 0);
	m_vvOffsets.clear(); m_vvOffsets.resize(num_subsets());

	// loop subsets
	for(int si = 0; si < num_subsets(); ++si)
	{
		// create offsets
		size_t count = 0;
		m_vvOffsets[si].resize(m_pFunctionPattern->num_fct());
		for(size_t fct = 0; fct < m_pFunctionPattern->num_fct(); ++fct)
		{
			if(!is_def_in_subset(fct, si)) m_vvOffsets[si][fct] = (size_t) -1;
			else m_vvOffsets[si][fct] = count++;
		}

		// skip if no dofs to be distributed
		if(!(m_pFunctionPattern->num_fct(si)>0)) continue;

		iterBegin = m_goc.begin<VertexBase>(si);
		iterEnd =  m_goc.end<VertexBase>(si);

		// loop Vertices
		i_per_subset = 0;
		for(iter = iterBegin; iter != iterEnd; ++iter)
		{
			// get vertex
			VertexBase* vrt = *iter;

			// write index
			m_pStorageManager->m_vSubsetInfo[si].aaDoFVRT[vrt] = i;
			i ++;
			i_per_subset ++;
		}
		m_vNumDoFs[si] = i_per_subset;
	}
	m_numDoFs = i;

	return true;
}

//////////////////////////////
//////////////////////////////
// Cuthill - McKee
//////////////////////////////
//////////////////////////////

struct VertexInfo
{
	VertexInfo() : pVertex(NULL), si(-1) {}
	VertexBase* pVertex;
	int si;
	size_t oldNr;
	std::vector<size_t> vAdjacentVertex;
	bool handled;

	size_t degree() const {return vAdjacentVertex.size();}
};

struct SortClass {
		SortClass(const std::vector<VertexInfo>& vInfo)
		 	: m_vInfo(vInfo) {}
		const std::vector<VertexInfo>& m_vInfo;

		bool operator() (size_t i,size_t j) { return (m_vInfo[i].degree()<m_vInfo[j].degree());}
};

struct NewInfo{
		NewInfo(VertexBase* vrt_, int si_)
		: vrt(vrt_), si(si_) {}
	VertexBase* vrt;
	int si;
};



bool
P1ConformDoFDistribution::
order_cuthill_mckee(bool bReverse)
{
	UG_ASSERT(m_pStorageManager != NULL, "No Storage Manager");
	UG_ASSERT(m_pISubsetHandler != NULL, "No Subset Handler");

	if(num_subsets() == 0)
	{
		UG_LOG("Cuthill_McKee: No subsets. Done.\n");
		return true;
	}

//	check that in all subsets same number of functions
	size_t num_fct = m_pFunctionPattern->num_fct(0);
	for(int si = 0; si < num_subsets(); ++si)
	{
		if(num_fct != m_pFunctionPattern->num_fct(si))
		{
			UG_LOG("Cuthill_McKee: Currently only implemented iff same"
					" number of functions in all subsets.\n");
			return true;
		}
	}

//	check that numbering is correct
	if(m_numDoFs % num_fct != 0)
	{
		UG_LOG("Cuthill_McKee: Cannot divide number of dofs / num fct. Done.\n");
		return false;
	}

//	Adjacent Edges
	std::vector<EdgeBase*> vEdges;

// 	Grid
	Grid* grid = m_pStorageManager->m_pSH->get_assigned_grid();

// 	Iterators
	geometry_traits<VertexBase>::iterator iter, iterBegin, iterEnd;

//	create list of current numbering
	std::vector<VertexInfo> vOldIndex(m_numDoFs / num_fct);

//	Loop vertices
	size_t ind = 0;
	for(int si = 0; si < num_subsets(); ++si)
	{
		iterBegin = m_goc.begin<VertexBase>(si);
		iterEnd =  m_goc.end<VertexBase>(si);

		for(iter = iterBegin; iter != iterEnd; ++iter)
		{
		// 	Get vertex
			VertexBase* vrt = *iter;

		//	Set infos
			vOldIndex[ind].pVertex = vrt;
			vOldIndex[ind].si = si;
			vOldIndex[ind].oldNr = m_pStorageManager->m_vSubsetInfo[si].aaDoFVRT[vrt] / num_fct;
			vOldIndex[ind].vAdjacentVertex.clear();
			vOldIndex[ind].handled = false;

		//	Get Edges
			CollectEdges(vEdges, *grid, vrt);

		//	Get connected indices
			for(size_t ed = 0; ed < vEdges.size(); ++ed)
			{
				for(size_t i = 0; i < vEdges[ed]->num_vertices(); ++i)
				{
				//	Get Vertices of adjacent edges
					VertexBase* vrt1 = vEdges[ed]->vertex(i);

				//	skip own vertex
					if(vrt1 == vrt) continue;

				//	get index
					const int si1 = m_pISubsetHandler->get_subset_index(vrt1);

				//	skip iff in no subset
				//  This can happen, when no subset is set from the beginning or
				//  even when a surface grid is considered with vertical
				//	copy nodes.
					if(si1 < 0) continue;

				//	get adjacent index
					const size_t adjInd = m_pStorageManager->m_vSubsetInfo[si1].aaDoFVRT[vrt1] / num_fct;

				//	Add vertex to list of vertices
					vOldIndex[ind].vAdjacentVertex.push_back(adjInd);
				}
			}
			ind ++;
		}
	}
	if(ind*num_fct != m_numDoFs)
	{
		UG_LOG("ERROR in order_cuthill_mckee: "
				"Sorting #ind = "<< ind * num_fct<<", but we have " << m_numDoFs << " Indices.\n");
		return false;
	}

//	Sort adjacent vertices by degree
	SortClass mySortClass(vOldIndex);
	for(size_t i = 0; i < vOldIndex.size(); ++i)
	{
		std::sort(	vOldIndex[i].vAdjacentVertex.begin(),
					vOldIndex[i].vAdjacentVertex.end(), mySortClass);
	}

// 	Create list of mapping
	std::vector<NewInfo> vNewIndex; vNewIndex.reserve(m_numDoFs / num_fct);

	while(true)
	{
	//	find first unhandled vertex
		size_t start = vOldIndex.size();
		for(size_t i = 0; i < vOldIndex.size(); ++i)
		{
			if(vOldIndex[i].handled == false)
			{
				start = i; break;
			}
		}

	//	check if one unhandled vertex left
		if(start == vOldIndex.size())
			break;

	//	Find node with smallest degree
		for(size_t i = 0; i < vOldIndex.size(); ++i)
		{
			if(vOldIndex[i].handled == false &&
				vOldIndex[i].degree() < vOldIndex[start].degree())
				start = i;
		}

	//	Add start vertex to mapping
		vNewIndex.push_back(NewInfo(vOldIndex[start].pVertex, vOldIndex[start].si));
		vOldIndex[start].handled = true;

	//	Create queue of adjacent vertices
		std::queue<size_t> qAdjacent;
		for(size_t i = 0; i < vOldIndex[start].degree(); ++i)
		{
			qAdjacent.push(vOldIndex[start].vAdjacentVertex[i]);
		}

	//	add adjacent vertices to mapping
		while(!qAdjacent.empty())
		{
		//	get next index
			const size_t front = qAdjacent.front();

		//	if not handled
			if(vOldIndex[front].handled == false)
			{
			//	Add to mapping
				vNewIndex.push_back(NewInfo(vOldIndex[front].pVertex, vOldIndex[front].si));
				vOldIndex[front].handled = true;

			//	add adjacent to queue
				for(size_t i = 0; i < vOldIndex[front].degree(); ++i)
				{
					const size_t ind = vOldIndex[front].vAdjacentVertex[i];
					if(vOldIndex[ind].handled == false)
						qAdjacent.push(ind);
				}
			}

		//	pop index
			qAdjacent.pop();
		}
	}

//	Check that all indices have been handled
	if(vNewIndex.size() != vOldIndex.size())
	{
		UG_LOG("ERROR in order_cuthill_mckee:"
				" Number of new Indices (" << vNewIndex.size() * num_fct << ") does not match number "
				" of old Indices (" << vOldIndex.size() * num_fct << ").\n");
		return false;
	}

//	substitute old new indices
	for(size_t i = 0; i < vNewIndex.size(); ++i)
	{
		VertexBase* vrt = vNewIndex[i].vrt;
		int si = vNewIndex[i].si;

		if(bReverse)
			m_pStorageManager->m_vSubsetInfo[si].aaDoFVRT[vrt] = (vNewIndex.size()-1-i) * num_fct;
		else
			m_pStorageManager->m_vSubsetInfo[si].aaDoFVRT[vrt] = i * num_fct;
	}

//	we're done
	return true;
}

bool
GroupedP1ConformDoFDistribution::
order_cuthill_mckee(bool bReverse)
{
	UG_ASSERT(m_pStorageManager != NULL, "No Storage Manager");
	UG_ASSERT(m_pISubsetHandler != NULL, "No Subset Handler");

//	Adjacend Edges
	std::vector<EdgeBase*> vEdges;

// 	Grid
	Grid* grid = m_pStorageManager->m_pSH->get_assigned_grid();

// 	Iterators
	geometry_traits<VertexBase>::iterator iter, iterBegin, iterEnd;

//	create list of current numbering
	std::vector<VertexInfo> vOldIndex(m_numDoFs);

//	Loop vertices
	size_t ind = 0;
	for(int si = 0; si < num_subsets(); ++si)
	{
		iterBegin = m_goc.begin<VertexBase>(si);
		iterEnd =  m_goc.end<VertexBase>(si);

		for(iter = iterBegin; iter != iterEnd; ++iter)
		{
		// 	Get vertex
			VertexBase* vrt = *iter;

		//	Set infos
			vOldIndex[ind].pVertex = vrt;
			vOldIndex[ind].si = si;
			vOldIndex[ind].oldNr = m_pStorageManager->m_vSubsetInfo[si].aaDoFVRT[vrt];
			vOldIndex[ind].vAdjacentVertex.clear();
			vOldIndex[ind].handled = false;

		//	Get Edges
			CollectEdges(vEdges, *grid, vrt);

		//	Get connected indices
			for(size_t ed = 0; ed < vEdges.size(); ++ed)
			{
				for(size_t i = 0; i < vEdges[ed]->num_vertices(); ++i)
				{
				//	Get Vertices of adjacent edges
					VertexBase* vrt1 = vEdges[ed]->vertex(i);

				//	skip own vertex
					if(vrt1 == vrt) continue;

				//	get index
					const int si1 = m_pISubsetHandler->get_subset_index(vrt1);

				//	skip iff in no subset
				//  This can happen, when no subset is set from the beginning or
				//  even when a surface grid is considered with vertical
				//	copy nodes.
					if(si1 < 0) continue;

				//	get adjacent index
					const size_t adjInd = m_pStorageManager->m_vSubsetInfo[si1].aaDoFVRT[vrt1];

				//	Add vertex to list of vertices
					vOldIndex[ind].vAdjacentVertex.push_back(adjInd);
				}
			}
		}
		ind++;
	}
//	Check
	if(ind != m_numDoFs)
	{
		UG_LOG("ERROR in order_cuthill_mckee: "
				"Sorting #ind = "<< ind <<", but we have " << m_numDoFs << " Indices.\n");
		return false;
	}

//	Sort adjacent vertices by degree
	SortClass mySortClass(vOldIndex);
	for(size_t i = 0; i < vOldIndex.size(); ++i)
	{
		std::sort(	vOldIndex[i].vAdjacentVertex.begin(),
					vOldIndex[i].vAdjacentVertex.end(), mySortClass);
	}

// 	Create list of mapping
	std::vector<NewInfo> vNewIndex; vNewIndex.reserve(m_numDoFs);

	while(true)
	{
	//	find first unhandled vertex
		size_t start = vOldIndex.size();
		for(size_t i = 0; i < vOldIndex.size(); ++i)
		{
			if(vOldIndex[i].handled == false)
			{
				start = i; break;
			}
		}

	//	check if one unhandled vertex left
		if(start == vOldIndex.size())
			break;

	//	Find node with smallest degree
		for(size_t i = 0; i < vOldIndex.size(); ++i)
		{
			if(vOldIndex[i].handled == false &&
				vOldIndex[i].degree() < vOldIndex[start].degree())
				start = i;
		}

	//	Add start vertex to mapping
		vNewIndex.push_back(NewInfo(vOldIndex[start].pVertex, vOldIndex[start].si));
		vOldIndex[start].handled = true;

	//	Create queue of adjacent vertices
		std::queue<size_t> qAdjacent;
		for(size_t i = 0; i < vOldIndex[start].degree(); ++i)
		{
			qAdjacent.push(vOldIndex[start].vAdjacentVertex[i]);
		}

	//	add adjacent vertices to mapping
		while(!qAdjacent.empty())
		{
		//	get next index
			const size_t front = qAdjacent.front();

		//	if not handled
			if(vOldIndex[front].handled == false)
			{
			//	Add to mapping
				vNewIndex.push_back(NewInfo(vOldIndex[front].pVertex, vOldIndex[front].si));
				vOldIndex[front].handled = true;

			//	add adjacent to queue
				for(size_t i = 0; i < vOldIndex[front].degree(); ++i)
				{
					const size_t ind = vOldIndex[front].vAdjacentVertex[i];
					if(vOldIndex[ind].handled == false)
						qAdjacent.push(ind);
				}
			}

		//	pop index
			qAdjacent.pop();
		}
	}

//	Check that all indices have been handled
	if(vNewIndex.size() != vOldIndex.size())
	{
		UG_LOG("ERROR in order_cuthill_mckee:"
				" Number of new Indices (" << vNewIndex.size() << ") does not match number "
				" of old Indices (" << vOldIndex.size() << ").\n");
		return false;
	}

//	substitute old new indices
	for(size_t i = 0; i < vNewIndex.size(); ++i)
	{
		VertexBase* vrt = vNewIndex[i].vrt;
		int si = vNewIndex[i].si;

		if(bReverse)
			m_pStorageManager->m_vSubsetInfo[si].aaDoFVRT[vrt] = (vNewIndex.size()-1-i);
		else
			m_pStorageManager->m_vSubsetInfo[si].aaDoFVRT[vrt] = i;
	}

//	we're done
	return true;
}




} // end namespace ug

