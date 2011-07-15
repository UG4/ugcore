#include "./p1conform.h"
#include <vector>
#include <queue>
#include <algorithm>

namespace ug{

///////////////////////////////////////////////////////////////////////////////
// P1StorageManager
///////////////////////////////////////////////////////////////////////////////

void P1StorageManager::set_subset_handler(ISubsetHandler& sh)
{
//	do nothing if same subset handler
	if(m_pSH == &sh) return;

//	clear first if already subset handler set
	clear();

//	set SubsetHandler and grid
	m_pSH = &sh;
	m_pGrid = m_pSH->get_assigned_grid();
}

void P1StorageManager::clear()
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

bool P1StorageManager::update_attachments()
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
// P1DoFDistribution
///////////////////////////////////////////////////////////////////////////////

template <bool bGrouped>
bool P1DoFDistribution<bGrouped>::has_dofs_on(GeometricBaseObject gbo) const
{
//	only in case of a Vertex, we have a DoF
	if(gbo == VERTEX) return true;
	else return false;
}

template <bool bGrouped>
bool P1DoFDistribution<bGrouped>::has_dofs_on(ReferenceObjectID roid) const
{
//	only in case of a Vertex, we have a DoF
	if(roid == ROID_VERTEX) return true;
	else return false;
}


template <bool bGrouped>
void P1DoFDistribution<bGrouped>::create_offsets()
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

template <bool bGrouped>
size_t P1DoFDistribution<bGrouped>::get_free_index(size_t si)
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
		if(!bGrouped) m_sizeIndexSet += num_fct(si);
		else ++m_sizeIndexSet;
	}

//	adjust counters
	if(!bGrouped)
	{
		m_numDoFs += num_fct(si);
		m_vNumDoFs[si] += num_fct(si);
	}
	else
	{
		++ m_numDoFs;
		++ (m_vNumDoFs[si]);
	}

//	return new index
	return freeIndex;
}

template <bool bGrouped>
void P1DoFDistribution<bGrouped>::push_free_index(size_t freeIndex, size_t si)
{
//	remember index
	m_vFreeIndex.push_back(freeIndex);

//	decrease number of distributed indices
	if(!bGrouped)
	{
		m_numDoFs -= num_fct(si);
		m_vNumDoFs[si] -= num_fct(si);
	}
	else
	{
		-- m_numDoFs;
		-- (m_vNumDoFs[si]);
	}
}

template <bool bGrouped>
bool P1DoFDistribution<bGrouped>::distribute_dofs()
{
//	storage manage required
	if(m_pStorageManager == NULL)
	{
		UG_LOG("In 'P1DoFDistribution::distribute_dofs:"
				"Storage Manager not set. Aborting.\n");
		return false;
	}

//	function pattern required
	if(this->m_pFuncPattern == NULL)
	{
		UG_LOG("In 'P1DoFDistribution::distribute_dofs:"
				"Function Pattern not set. Aborting.\n");
		return false;
	}

// 	Attach indices
	if(!m_pStorageManager->update_attachments()) return false;

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

	//	get iterators of vertices
		iterBegin = this->template begin<VertexBase>(si);
		iterEnd =  this->template end<VertexBase>(si);

	// 	remember number of functions
		const size_t numFct = num_fct(si);

	// 	loop Vertices
		m_vNumDoFs[si] = 0;
		for(iter = iterBegin; iter != iterEnd; ++iter)
		{
		// 	get vertex
			VertexBase* vrt = *iter;

		//	skip shadows
			if(this->m_pSurfaceView != NULL)
				if(this->m_pSurfaceView->is_shadow(vrt))
					continue;

		// 	write index
			first_index(vrt) = m_numDoFs;

		//	increase number of DoFs and DoFs on subset
			if(!bGrouped)
			{
				m_numDoFs += numFct;
				m_vNumDoFs[si] += numFct;
			}
			else
			{
				++ m_numDoFs;
				++ (m_vNumDoFs[si]);
			}
		}
	}

//	the size of the index set is the number of DoFs
	m_sizeIndexSet = m_numDoFs;

//	post process shadows.
	if(this->m_pSurfaceView != NULL)
	{
	//	loop subsets
		for(int si = 0; si < num_subsets(); ++si)
		{
		//	iterators
			iterBegin = this->template begin<VertexBase>(si);
			iterEnd =  this->template end<VertexBase>(si);

		//	loop vertices
			for(iter = iterBegin; iter != iterEnd; ++iter)
			{
			//	get vertex
				VertexBase* vrt = *iter;

			//	skip non-shadows
				if(!this->m_pSurfaceView->is_shadow(vrt)) continue;

			//	get child
				VertexBase* vrtChild = this->m_pSurfaceView->template get_shadow_child<VertexBase>(vrt);

			//	get indices of child
				const size_t indexChild = first_index(vrtChild);

			//	set index of shadow to index of child
				first_index(vrt) = indexChild;
			}
		}
	}

//	we're done
	return true;
}

template <bool bGrouped>
bool P1DoFDistribution<bGrouped>::
permute_indices(std::vector<size_t>& vIndNew)
{
//	storage manager required
	if(m_pStorageManager == NULL)
	{
		UG_LOG("ERROR in 'P1DoFDistribution::permute_indices':"
				" No Storage Manager");
		return false;
	}

//	subsethandler required
	if(m_pISubsetHandler == NULL)
	{
		UG_LOG("ERROR in 'P1DoFDistribution::permute_indices':"
				" No Subset Handler");
		return false;
	}

//	check, that passed index fields have the same size
	if(this->num_dofs() != vIndNew.size())
	{
		UG_LOG("ERROR in 'P1DoFDistribution::permute_indices': New index set"
				" must have same cardinality for swap indices.\n");
		return false;
	}

// 	loop subsets
	for(int si = 0; si < num_subsets(); ++si)
	{
	//	get iterators
		geometry_traits<VertexBase>::iterator iter, iterBegin, iterEnd;
		iterBegin = this->template begin<VertexBase>(si);
		iterEnd =  this->template end<VertexBase>(si);

	// 	loop Vertices
		for(iter = iterBegin; iter != iterEnd; ++iter)
		{
		// 	get vertex
			VertexBase* vrt = *iter;

		// 	get current (old) index
			const size_t oldIndex = first_index(vrt);

		//	replace old index by new one
			first_index(vrt) = vIndNew[oldIndex];

		//	set correct permutation in permuting vector
			if(!bGrouped)
				for(size_t fct = 0; fct < num_fct(si); ++fct)
					vIndNew[oldIndex+fct] = vIndNew[oldIndex] + fct;
		}
	}

//	we're done
	return true;
}

template <bool bGrouped>
bool P1DoFDistribution<bGrouped>::
get_connections(std::vector<std::vector<size_t> >& vvConnection)
{
//	storage manager required
	if(m_pStorageManager == NULL)
	{
		UG_LOG("ERROR in 'P1DoFDistribution::get_connections':"
				" No Storage Manager");
		return false;
	}

//	subset handler required
	if(m_pISubsetHandler == NULL)
	{
		UG_LOG("ERROR in 'P1DoFDistribution::get_connections':"
				" No Subset Handler");
		return false;
	}

//	clear neighbours
	vvConnection.clear(); vvConnection.resize(m_numDoFs);

//	if no subset given, we're done
	if(num_subsets() == 0) return true;

	const size_t numFct = num_fct(0);
//	check that in all subsets same number of functions and at least one
//	if this is not the case for ungrouped DoFs, we cannot allow reordering
	if(!bGrouped)
	{
		for(int si = 0; si < num_subsets(); ++si)
		{
			if(numFct != num_fct(si))
			{
				UG_LOG("ERROR in 'P1DoFDistribution::get_connections':"
						" Currently only implemented iff same number of functions"
						" in all subsets.\n");
				return true;
			}
		}
	//	if no functions given, return
		if(numFct == 0) return true;
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
	//	iterators
		iterBegin = this->template begin<VertexBase>(si);
		iterEnd =  this->template end<VertexBase>(si);

	//	loop vertices
		for(iter = iterBegin; iter != iterEnd; ++iter)
		{
		// 	Get vertex
			VertexBase* vrt = *iter;

		//	skip shadows
			if(this->m_pSurfaceView != NULL)
				if(this->m_pSurfaceView->is_shadow(vrt))
					continue;

		//	get index
			const size_t index = first_index(vrt);

		//	always connection with itself
			vvConnection[index].push_back(index);

		//	Get Edges
			CollectEdges(vEdges, *grid, vrt);

		//	Get connections via shadow
			if(this->m_pSurfaceView != NULL)
			{
			//	skip if not shadowing
				if(this->m_pSurfaceView->shadows(vrt))
				{
				//	get parent
					VertexBase* vrtParent =
						dynamic_cast<VertexBase*>(this->m_pSurfaceView->get_parent(vrt));

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
					const size_t adjInd = first_index(vrt1);

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

template <bool bGrouped>
void P1DoFDistribution<bGrouped>::grid_obj_added(VertexBase* vrt)
{
//	subset handler required
	if(m_pISubsetHandler == NULL)
	{
		UG_LOG("ERROR in 'P1DoFDistribution::get_connections':"
				" No Subset Handler");
		throw(UGFatalError("Subset Handler missing."));
	}

//	for newly created vertices, we have to add the integers

//	get subset index
	const int si = m_pISubsetHandler->get_subset_index(vrt);

//	Only for the surface dof manager:
//	If the added vertex is a shadow, use index of child (shadowing) vertex
	if(this->m_pSurfaceView != NULL && this->m_pSurfaceView->is_shadow(vrt))
	{
	//	get child
		VertexBase* vrtChild = this->m_pSurfaceView->template get_shadow_child<VertexBase>(vrt);

	//	get indices of child
		const size_t indexChild = first_index(vrtChild);

	//	set index of shadow to index of child
		first_index(vrt) = indexChild;
	}
//	normally, we can set a new free index
	else
	{
	// 	write next free index
		first_index(vrt) = get_free_index(si);
	}
}

template <bool bGrouped>
void P1DoFDistribution<bGrouped>::grid_obj_to_be_removed(VertexBase* vrt)
{
//	get subset index
	const int si = m_pISubsetHandler->get_subset_index(vrt);

// 	remember free index
	push_free_index(first_index(vrt), si);
}

template <bool bGrouped>
void P1DoFDistribution<bGrouped>::
grid_obj_replaced(VertexBase* vrtNew, VertexBase* vrtOld)
{
	UG_ASSERT(	m_pISubsetHandler->get_subset_index(vrtNew) ==
				m_pISubsetHandler->get_subset_index(vrtOld),
	          "New vertex does not have same subset as replaced on.");

//	copy index
	first_index(vrtNew) = first_index(vrtOld);
}


template <bool bGrouped>
bool P1DoFDistribution<bGrouped>::defragment()
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
		iterBegin = this->template begin<VertexBase>(si);
		iterEnd =  this->template end<VertexBase>(si);

	// 	loop Vertices
		for(iter = iterBegin; iter != iterEnd; ++iter)
		{
		// 	get vertex
			VertexBase* vrt = *iter;

		//	get old (current) index
			const size_t oldIndex = first_index(vrt);

		// 	check if index must be replaced by lower one
			if(first_index(vrt) < m_numDoFs) continue;

		//	get free index
			const size_t newIndex = get_free_index(si);

		//	remember replacement
			m_vReplaced.push_back(std::pair<size_t,size_t>(oldIndex, newIndex));

		//	overwrite index
			first_index(vrt) = newIndex;

		//	adjust counters, since this was an replacement
			if(!bGrouped)
			{
				m_numDoFs -= num_fct(si);
				m_sizeIndexSet -= num_fct(si);
				m_vNumDoFs[si] -= num_fct(si);
			}
			else
			{
				--m_numDoFs;
				--m_sizeIndexSet;
				--(m_vNumDoFs[si]);
			}
		}
	}

//	check that all holes have been removed
	if(m_numDoFs != m_sizeIndexSet)
	{
		UG_LOG("ERROR in 'P1DoFDistribution::defragment': Still holes in index "
				"set after compression. Check implementation.\n");
		return false;
	}

//	copy values (from back into holes) for managed grid functions
	if(!indices_swaped(m_vReplaced, true)) return false;

//	cut of unused tail of managed grid functions
	num_indices_changed(m_numDoFs);

//	we're done
	return true;
}

// explicit template instantiations
template class P1DoFDistribution<true>;
template class P1DoFDistribution<false>;

} // end namespace ug

