#include "./p1conform.h"
#include <vector>
#include <queue>
#include <algorithm>

namespace ug{

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// P1StorageManager
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

void
P1StorageManager::
set_subset_handler(ISubsetHandler& sh)
{
//	do nothing if same subset handler
	if(m_pSH == &sh) return;

//	clear first if already subset handler set
	clear();

//	set SubsetHandler and grid
	m_pSH = &sh;
	m_pGrid = m_pSH->get_assigned_grid();
}

void
P1StorageManager::
clear()
{
//	if no subsethandler given, nothing is attached
	if(m_pSH == NULL) return;

//	detach DoFs
	m_pGrid->detach_from<VertexBase>(m_aDoF);
	m_aaDoFVRT.invalidate();

//	reset SubsetHandler
	m_pSH = NULL;
	m_pGrid = NULL;
}

bool
P1StorageManager::
update_attachments()
{
//	check, that everything has been set
	if(m_pSH == NULL)
	{
		UG_LOG("ERROR in 'P1StorageManager::update_attachments':"
				" Updating indices, but no SubsetHandler set.\n");
		return false;
	}
	if(m_pGrid == NULL)
	{
		UG_LOG("ERROR in 'P1StorageManager::update_attachments':"
				" Updating indices, but no Grid in SubsetHandler set.\n");
		return false;
	}

//	attach DoFs to vertices
	m_pGrid->attach_to<VertexBase>(m_aDoF);

//	access the
	m_aaDoFVRT.access(*m_pGrid, m_aDoF);

//	done
	return true;
}



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// P1ConformDoFDistribution
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

size_t
P1ConformDoFDistribution::
get_free_index(size_t si)
{
//	The idea is as follows:
//	- 	If a free index is left, the a free index is returned. This
//		changes the number of (used) dofs, but the index set remains
//		the same. (one hole less)
// 	-	If no free index is left (i.e. no holes in index set and therefore
//		m_numDoFs == m_sizeIndexSet), the index set is increased and
//		the newly created index is returned. This changes the size of
//		the index set and the number of dofs.

//	start with default index to be returned
	size_t freeIndex = m_sizeIndexSet;

//	check if free index available
	if(!m_vFreeIndex.empty())
	{
	//	return free index instead and pop index from free index list
		freeIndex = m_vFreeIndex.back(); m_vFreeIndex.pop_back();
	}
	else
	{
	//	if using new index, increase size of index set
		m_sizeIndexSet += num_fct(si);
	}

//	adjust counters
	m_numDoFs += num_fct(si);
	m_vNumDoFs[si] += num_fct(si);

//	return new index
	return freeIndex;
}

void
P1ConformDoFDistribution::
push_free_index(size_t freeIndex, size_t si)
{
//	remember index
	m_vFreeIndex.push_back(freeIndex);

//	decrease number of distributed indices
	m_numDoFs -= num_fct(si);
	m_vNumDoFs[si] -= num_fct(si);
}

void
P1ConformDoFDistribution::
create_offsets()
{
//	clear offsets
	m_vvOffsets.clear();

//	resize for all subsets
	m_vvOffsets.resize(num_subsets());

// 	loop subsets
	for(int si = 0; si < num_subsets(); ++si)
	{
	// 	counter
		size_t count = 0;

	//	resize for each function
		m_vvOffsets[si].resize(num_fct());

	//	loop functions
		for(size_t fct = 0; fct < num_fct(); ++fct)
		{
		//	set offset for each function defined in the subset
			if(!is_def_in_subset(fct, si)) m_vvOffsets[si][fct] = (size_t) -1;
			else m_vvOffsets[si][fct] = count++;
		}
	}
}

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

	if(m_pFuncPattern == NULL)
	{
		UG_LOG("In 'P1ConformDoFDistribution::distribute_dofs:"
				"Function Pattern not set. Aborting.\n");
		return false;
	}

// 	Attach indices
	if(!m_pStorageManager->update_attachments())
		return false;

// 	iterators
	geometry_traits<VertexBase>::iterator iter, iterBegin, iterEnd;

// 	create offsets
	create_offsets();

// 	reset counter for all dofs
	m_numDoFs = 0;

// 	reset number of dofs
	m_vNumDoFs.clear(); m_vNumDoFs.resize(num_subsets(), 0);

// 	loop subsets
	for(int si = 0; si < num_subsets(); ++si)
	{
	// 	skip if no dofs to be distributed
		if(!(num_fct(si)>0)) continue;

		iterBegin = this->begin<VertexBase>(si);
		iterEnd =  this->end<VertexBase>(si);

	// 	remember number of functions
		const size_t numFct = num_fct(si);

	// 	loop Vertices
		m_vNumDoFs[si] = 0;
		for(iter = iterBegin; iter != iterEnd; ++iter)
		{
		// 	get vertex
			VertexBase* vrt = *iter;

		//	skip shadows
			if(m_pSurfaceView != NULL)
				if(m_pSurfaceView->is_shadow(vrt))
					continue;

		// 	write index
			first_index(vrt, si) = m_numDoFs;

		//	increase number of DoFs
			m_numDoFs += numFct;

		//	increase number of dofs on subset
			m_vNumDoFs[si] += numFct;
		}
	}

//	the size of the index set is the number of DoFs
	m_sizeIndexSet = m_numDoFs;

//	post process shadows.
	if(m_pSurfaceView != NULL)
	{
		for(int si = 0; si < num_subsets(); ++si)
		{

			iterBegin = this->begin<VertexBase>(si);
			iterEnd =  this->end<VertexBase>(si);
			for(iter = iterBegin; iter != iterEnd; ++iter)
			{
			//	get vertex
				VertexBase* vrt = *iter;

			//	skip non-shadows
				if(!m_pSurfaceView->is_shadow(vrt)) continue;

			//	get child
				VertexBase* vrtChild = m_pSurfaceView->get_child(vrt);

			//	get indices of child
				const size_t indexChild = first_index(vrtChild, si);

			//	set index of shadow to index of child
				first_index(vrt, si) = indexChild;
			}
		}
	}

//	we're done
	return true;
}

bool
P1ConformDoFDistribution::
permute_indices(std::vector<size_t>& vIndNew)
{
//	check, that storage is initialized
	if(m_pStorageManager == NULL)
	{
		UG_LOG("ERROR in 'P1ConformDoFDistribution::permute_indices':"
				" No Storage Manager");
		return false;
	}
	if(m_pISubsetHandler == NULL)
	{
		UG_LOG("ERROR in 'P1ConformDoFDistribution::permute_indices':"
				" No Subset Handler");
		return false;
	}

//	check, that passed index fields have the same size
	if(this->num_dofs() != vIndNew.size())
	{
		UG_LOG("ERROR in 'P1ConformDoFDistribution::permute_indices': New index set"
				" must have same cardinality for swap indices.\n");
		return false;
	}

// 	loop subsets
	for(int si = 0; si < num_subsets(); ++si)
	{
	//	get iterators
		geometry_traits<VertexBase>::iterator iter, iterBegin, iterEnd;
		iterBegin = this->begin<VertexBase>(si);
		iterEnd =  this->end<VertexBase>(si);

	// 	loop Vertices
		for(iter = iterBegin; iter != iterEnd; ++iter)
		{
		// 	get vertex
			VertexBase* vrt = *iter;

		// 	get current (old) index
			const size_t oldIndex = first_index(vrt, si);

		//	replace old index by new one
			first_index(vrt, si) = vIndNew[oldIndex];

		//	set correct permutation in permuting vector
			for(size_t fct = 0; fct < num_fct(si); ++fct)
				vIndNew[oldIndex+fct] = vIndNew[oldIndex] + fct;
		}
	}

//	we're done
	return true;
}

bool
P1ConformDoFDistribution::
get_connections(std::vector<std::vector<size_t> >& vvConnection)
{
	UG_ASSERT(m_pStorageManager != NULL, "No Storage Manager");
	UG_ASSERT(m_pISubsetHandler != NULL, "No Subset Handler");

//	clear neighbours
	vvConnection.clear(); vvConnection.resize(m_numDoFs);

//	if no subset given, we're done
	if(num_subsets() == 0) return true;

//	check that in all subsets same number of functions and at least one
	size_t numFct = num_fct(0);
	if(numFct == 0) return true;
	for(int si = 0; si < num_subsets(); ++si)
	{
		if(numFct != num_fct(si))
		{
			UG_LOG("ERROR in 'P1ConformDoFDistribution::get_algebraic_neighbours':"
					" Currently only implemented iff same number of functions"
					" in all subsets.\n");
			return true;
		}
	}

//	Adjacent Edges
	std::vector<EdgeBase*> vEdges;

// 	Grid
	Grid* grid = m_pStorageManager->get_assigned_grid();

// 	Iterators
	geometry_traits<VertexBase>::iterator iter, iterBegin, iterEnd;

//	Loop vertices
	for(int si = 0; si < num_subsets(); ++si)
	{
		iterBegin = this->begin<VertexBase>(si);
		iterEnd =  this->end<VertexBase>(si);

		for(iter = iterBegin; iter != iterEnd; ++iter)
		{
		// 	Get vertex
			VertexBase* vrt = *iter;

		//	skip shadows
			if(m_pSurfaceView != NULL)
				if(m_pSurfaceView->is_shadow(vrt))
					continue;

		//	get index
			const size_t index = first_index(vrt, si);

		//	always connection with itself
			vvConnection[index].push_back(index);

		//	Get Edges
			CollectEdges(vEdges, *grid, vrt);

		//	Get connections via shadow
			if(m_pSurfaceView != NULL)
			{
			//	skip if not shadowing
				if(m_pSurfaceView->shadows(vrt))
				{
				//	get parent
					VertexBase* vrtParent =
						dynamic_cast<VertexBase*>(m_pSurfaceView->get_parent(vrt));

				//	Get Edges
					if(vrtParent != NULL)
						CollectEdges(vEdges, *grid, vrtParent, false);
				}
			}

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
					const size_t adjInd = first_index(vrt1, si1);

				//	search for index in adjacend indices
					std::vector<size_t>::iterator it;
					it = find(vvConnection[index].begin(), vvConnection[index].end(),
					          adjInd);

				//	Add vertex to list of vertices
					if(it == vvConnection[index].end())
						vvConnection[index].push_back(adjInd);
				}
			}
		}
	}

//	we're done
	return true;
}

void
P1ConformDoFDistribution::
grid_obj_added(VertexBase* vrt)
{
	UG_ASSERT(m_pISubsetHandler != NULL, "No Subset Handler.");

//	for newly created vertices, we have to add the integers

//	get subset index
	const int si = m_pISubsetHandler->get_subset_index(vrt);

//	Only for the surface dof manager:
//	If the added vertex is a shadow, use index of child (shadowing) vertex
	if(m_pSurfaceView != NULL && m_pSurfaceView->is_shadow(vrt))
	{
	//	get child
		VertexBase* vrtChild = m_pSurfaceView->get_child(vrt);

	//	get indices of child
		const size_t indexChild = first_index(vrtChild, si);

	//	set index of shadow to index of child
		first_index(vrt, si) = indexChild;
	}
//	normally, we can set a new free index
	else
	{
		const size_t index = get_free_index(si);

	// 	write next free index
		first_index(vrt, si) = index;
	}
}

void
P1ConformDoFDistribution::
grid_obj_to_be_removed(VertexBase* vrt)
{
//	get subset index
	const int si = m_pISubsetHandler->get_subset_index(vrt);

// 	remember free index
	push_free_index(first_index(vrt, si), si);
}

void
P1ConformDoFDistribution::
grid_obj_replaced(VertexBase* vrtNew, VertexBase* vrtOld)
{
//	get subset index
	const int si = m_pISubsetHandler->get_subset_index(vrtOld);

	UG_ASSERT(m_pISubsetHandler->get_subset_index(vrtNew) == si,
	          "New vertex does not have same subset as replaced on.");

//	copy index
	first_index(vrtNew, si) = first_index(vrtOld, si);
}


bool
P1ConformDoFDistribution::
defragment()
{
//	we loop all indices and those with highest degree are replaced by a free
//	index. This operation is performed in O(number of Indices)
//	All indices with an index >= m_numDoFs are replaced.

//	check, if holes exist. If not, we're done
	if(m_vFreeIndex.empty())
	{
	//	the index set might have been increased. Thus resize grid functions
		num_indices_changed(m_numDoFs);

	//	we're done
		return true;
	}

//	pairs of replaced indices
	std::vector<std::pair<size_t, size_t> > m_vReplaced;

	for(int si = 0; si < num_subsets(); ++si)
	{
	// 	skip if no DoFs to be distributed
		if(!(num_fct(si)>0)) continue;

		geometry_traits<VertexBase>::iterator iter, iterBegin, iterEnd;
		iterBegin = this->begin<VertexBase>(si);
		iterEnd =  this->end<VertexBase>(si);

	// 	loop Vertices
		for(iter = iterBegin; iter != iterEnd; ++iter)
		{
		// 	get vertex
			VertexBase* vrt = *iter;

		//	get old (current) index
			const size_t oldIndex = first_index(vrt, si);

		// 	check if index must be replaced by lower one
			if(first_index(vrt, si) < m_numDoFs) continue;

		//	get free index
			const size_t newIndex = get_free_index(si);

		//	remember replacement
			m_vReplaced.push_back(std::pair<size_t,size_t>(oldIndex, newIndex));

		//	overwrite index
			first_index(vrt, si) = newIndex;

		//	adjust counters, since this was an replacement
			m_numDoFs -= num_fct(si);
			m_sizeIndexSet -= num_fct(si);
			m_vNumDoFs[si] -= num_fct(si);
		}
	}

	UG_LOG(" swapping done.\n");

//	check that all holes have been removed
	if(m_numDoFs != m_sizeIndexSet)
	{
		UG_LOG("ERROR in 'GroupedP1ConformDoFDistribution::compress': Still "
				" holes in index set after compression. Check implementation.\n");
		return false;
	}

//	copy values (from back into holes) for managed grid functions
	//\todo: HANDLE ALL INDICES
	if(!indices_swaped(m_vReplaced, true)) return false;

//	cut of unused tail of managed grid functions
	num_indices_changed(m_numDoFs);

//	we're done
	return true;
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// GroupedP1ConformDoFDistribution
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


size_t
GroupedP1ConformDoFDistribution::
get_free_index(size_t si)
{
//	The idea is as follows:
//	- 	If a free index is left, the a free index is returned. This
//		changes the number of (used) dofs, but the index set remains
//		the same. (one hole less)
// 	-	If no free index is left (i.e. no holes in index set and therefore
//		m_numDoFs == m_sizeIndexSet), the index set is increased and
//		the newly created index is returned. This changes the size of
//		the index set and the number of dofs.

//	strat with default index to be returned
	size_t freeIndex = m_sizeIndexSet;

//	check if free index available
	if(!m_vFreeIndex.empty())
	{
	//	return free index instead and pop index from free index list
		freeIndex = m_vFreeIndex.back(); m_vFreeIndex.pop_back();
	}
	else
	{
	//	if using new index, increase size of index set
		++m_sizeIndexSet;
	}

//	adjust counters
	++ m_numDoFs;
	++ (m_vNumDoFs[si]);

//	return new index
	return freeIndex;
}

void
GroupedP1ConformDoFDistribution::
push_free_index(size_t freeIndex, size_t si)
{
//	remember index
	m_vFreeIndex.push_back(freeIndex);

//	decrease number of distributed indices
	-- m_numDoFs;
	-- (m_vNumDoFs[si]);
}


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

	if(m_pFuncPattern == NULL)
	{
		UG_LOG("In 'P1ConformDoFDistribution::distribute_dofs:"
				"Function Pattern not set. Aborting.\n");
		return false;
	}

// 	Attach indices
	if(!m_pStorageManager->update_attachments())
		return false;

// 	iterators
	geometry_traits<VertexBase>::iterator iter, iterBegin, iterEnd;

//	counters
	size_t i = 0;

// 	reset number of DoFs
	m_vNumDoFs.clear(); m_vNumDoFs.resize(num_subsets(), 0);

// 	loop subsets
	for(int si = 0; si < num_subsets(); ++si)
	{
	// 	skip if no DoFs to be distributed
		if(!(num_fct(si)>0)) continue;

		iterBegin = this->begin<VertexBase>(si);
		iterEnd =  this->end<VertexBase>(si);

	// 	loop Vertices
		m_vNumDoFs[si] = 0;
		for(iter = iterBegin; iter != iterEnd; ++iter)
		{
		// 	get vertex
			VertexBase* vrt = *iter;

		// 	write index
			alg_index(vrt, si) = i ++;

		//	increase number of DoFs in subset
			m_vNumDoFs[si] ++;
		}
	}

//	remember total number of DoFs
	m_numDoFs = i;

//	the size of the index set is the number of DoFs
	m_sizeIndexSet = m_numDoFs;

//	we're done
	return true;
}

bool
GroupedP1ConformDoFDistribution::
permute_indices(std::vector<size_t>& vIndNew)
{
//	check, that storage is initialized
	if(m_pStorageManager == NULL)
	{
		UG_LOG("ERROR in 'GroupedP1ConformDoFDistribution::permute_indices':"
				" No Storage Manager");
		return false;
	}
	if(m_pISubsetHandler == NULL)
	{
		UG_LOG("ERROR in 'GroupedP1ConformDoFDistribution::permute_indices':"
				" No Subset Handler");
		return false;
	}

//	check, that passed index fields have the same size
	if(this->num_dofs() != vIndNew.size())
	{
		UG_LOG("ERROR in 'GroupedP1ConformDoFDistribution::permute_indices': "
				"New index set must have same cardinality for swap indices.\n");
		return false;
	}

// 	loop subsets
	for(int si = 0; si < num_subsets(); ++si)
	{
	//	get iterators
		geometry_traits<VertexBase>::iterator iter, iterBegin, iterEnd;
		iterBegin = this->begin<VertexBase>(si);
		iterEnd =  this->end<VertexBase>(si);

	// 	loop Vertices
		for(iter = iterBegin; iter != iterEnd; ++iter)
		{
		// 	get vertex
			VertexBase* vrt = *iter;

		// 	get current (old) index
			const size_t oldIndex = alg_index(vrt, si);

		//	replace old index by new one
			alg_index(vrt, si) = vIndNew[oldIndex];
		}
	}

//	we're done
	return true;
}

bool
GroupedP1ConformDoFDistribution::
get_connections(std::vector<std::vector<size_t> >& vvConnection)
{
	UG_ASSERT(m_pStorageManager != NULL, "No Storage Manager");
	UG_ASSERT(m_pISubsetHandler != NULL, "No Subset Handler");

//	clear neighbours
	vvConnection.clear(); vvConnection.resize(m_numDoFs);

//	if no subset given, we're done
	if(num_subsets() == 0) return true;

//	check that in all subsets same number of functions and at least one
	size_t numFct = num_fct(0);
	if(numFct == 0) return true;

//	Adjacent Edges
	std::vector<EdgeBase*> vEdges;

// 	Grid
	Grid* grid = m_pStorageManager->get_assigned_grid();

// 	Iterators
	geometry_traits<VertexBase>::iterator iter, iterBegin, iterEnd;

//	Loop vertices
	for(int si = 0; si < num_subsets(); ++si)
	{
		iterBegin = this->begin<VertexBase>(si);
		iterEnd =  this->end<VertexBase>(si);

		for(iter = iterBegin; iter != iterEnd; ++iter)
		{
		// 	Get vertex
			VertexBase* vrt = *iter;

		//	skip shadows
			if(m_pSurfaceView != NULL)
				if(m_pSurfaceView->is_shadow(vrt))
					continue;

		//	get index
			const size_t index = alg_index(vrt, si);

		//	always connection with itself
			vvConnection[index].push_back(index);

		//	Get Edges
			CollectEdges(vEdges, *grid, vrt);

		//	Get connections via shadow
			if(m_pSurfaceView != NULL)
			{
			//	skip if not shadowing
				if(m_pSurfaceView->shadows(vrt))
				{
				//	get parent
					VertexBase* vrtParent =
						dynamic_cast<VertexBase*>(m_pSurfaceView->get_parent(vrt));

				//	Get Edges
					if(vrtParent != NULL)
						CollectEdges(vEdges, *grid, vrtParent, false);
				}
			}

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
					const size_t adjInd = alg_index(vrt1, si1);

				//	search for index in adjacend indices
					std::vector<size_t>::iterator it;
					it = find(vvConnection[index].begin(), vvConnection[index].end(),
					          adjInd);

				//	Add vertex to list of vertices
					if(it == vvConnection[index].end())
						vvConnection[index].push_back(adjInd);
				}
			}
		}
	}

//	we're done
	return true;
}

void
GroupedP1ConformDoFDistribution::
grid_obj_added(VertexBase* vrt)
{
//	for newly created vertices, we have to add the integers

//	get subset index
	const int si = m_pISubsetHandler->get_subset_index(vrt);

//	Only for the surface dof manager:
//	If the added vertex is a shadow, use index of child (shadowing) vertex
	if(m_pSurfaceView != NULL && m_pSurfaceView->is_shadow(vrt))
	{
	//	get child
		VertexBase* vrtChild = m_pSurfaceView->get_child(vrt);

	//	get indices of child
		const size_t indexChild = alg_index(vrtChild, si);

	//	set index of shadow to index of child
		alg_index(vrt, si) = indexChild;
	}
//	normally, we can set a new free index
	else
	{
	// 	write next free index
		alg_index(vrt, si) = get_free_index(si);
	}
}

void
GroupedP1ConformDoFDistribution::
grid_obj_to_be_removed(VertexBase* vrt)
{
//	for vertices that will be erased, we have to remember the removed indices

//	get subset index
	const int si = m_pISubsetHandler->get_subset_index(vrt);

// 	remember free index
	push_free_index(alg_index(vrt, si), si);
}

void
GroupedP1ConformDoFDistribution::
grid_obj_replaced(VertexBase* vrtNew, VertexBase* vrtOld)
{
//	get subset index
	const int si = m_pISubsetHandler->get_subset_index(vrtOld);

	UG_ASSERT(m_pISubsetHandler->get_subset_index(vrtNew) == si,
	          "New vertex does not have same subset as replaced on.");

//	copy index
	alg_index(vrtNew, si) = alg_index(vrtOld, si);
}

bool
GroupedP1ConformDoFDistribution::
defragment()
{
//	we loop all indices and those with highest degree are replaced by a free
//	index. This operation is performed in O(number of Indices)
//	All indices with an index >= m_numDoFs are replaced.

//	check, if holes exist. If not, we're done
	if(m_vFreeIndex.empty())
	{
	//	the index set might have been increased. Thus resize grid functions
		num_indices_changed(m_numDoFs);

	//	we're done
		return true;
	}

//	pairs of replaced indices
	std::vector<std::pair<size_t, size_t> > m_vReplaced;

	for(int si = 0; si < num_subsets(); ++si)
	{
	// 	skip if no DoFs to be distributed
		if(!(num_fct(si)>0)) continue;

		geometry_traits<VertexBase>::iterator iter, iterBegin, iterEnd;
		iterBegin = this->begin<VertexBase>(si);
		iterEnd =  this->end<VertexBase>(si);

	// 	loop Vertices
		for(iter = iterBegin; iter != iterEnd; ++iter)
		{
		// 	get vertex
			VertexBase* vrt = *iter;

		//	get old (current) index
			const size_t oldIndex = alg_index(vrt, si);

		// 	check if index must be replaced by lower one
			if(alg_index(vrt, si) < m_numDoFs) continue;

		//	get free index
			const size_t newIndex = get_free_index(si);

		//	remember replacement
			m_vReplaced.push_back(std::pair<size_t,size_t>(oldIndex, newIndex));

		//	overwrite index
			alg_index(vrt, si) = newIndex;

		//	adjust counters, since this was an replacement
			--m_numDoFs;
			--m_sizeIndexSet;
			--(m_vNumDoFs[si]);
		}
	}

//	check that all holes have been removed
	if(m_numDoFs != m_sizeIndexSet)
	{
		UG_LOG("ERROR in 'GroupedP1ConformDoFDistribution::compress': Still "
				" holes in index set after compression. Check implementation.\n");
		return false;
	}

//	copy values (from back into holes) for managed grid functions
	if(!indices_swaped(m_vReplaced, true)) return false;

//	cut of unused tail of managed grid functions
	num_indices_changed(m_numDoFs);

//	we're done
	return true;
}

} // end namespace ug

