
//	standard included
#include <vector>
#include <queue>
#include <algorithm>

//	header
#include "conform.h"
#include "lib_discretization/local_finite_element/local_dof_set.h"

namespace ug{

///////////////////////////////////////////////////////////////////////////////
// ConformStorageManager
///////////////////////////////////////////////////////////////////////////////

void ConformStorageManager::set_subset_handler(ISubsetHandler& sh)
{
//	do nothing if same subset handler
	if(m_pSH == &sh) return;

//	clear first if already subset handler set
	clear();

//	set SubsetHandler and grid
	m_pSH = &sh;
	m_pGrid = m_pSH->get_assigned_grid();
}

void ConformStorageManager::clear()
{
//	if no subsethandler given, nothing is attached
	if(m_pSH == NULL) return;

//	detach DoFs
	m_pGrid->detach_from<VertexBase>(m_aIndex);
	m_pGrid->detach_from<EdgeBase>(m_aIndex);
	m_pGrid->detach_from<Face>(m_aIndex);
	m_pGrid->detach_from<Volume>(m_aIndex);
	m_aaIndexVRT.invalidate();
	m_aaIndexEDGE.invalidate();
	m_aaIndexFACE.invalidate();
	m_aaIndexVOL.invalidate();

//	reset SubsetHandler
	m_pSH = NULL;
	m_pGrid = NULL;
}

bool ConformStorageManager::update_attachments()
{
//	check, that everything has been set
	if(m_pSH == NULL)
	{
		UG_LOG("ERROR in 'ConformStorageManager::update_attachments':"
				" Updating indices, but no SubsetHandler set.\n");
		return false;
	}
	if(m_pGrid == NULL)
	{
		UG_LOG("ERROR in 'ConformStorageManager::update_attachments':"
				" Updating indices, but no Grid in SubsetHandler set.\n");
		return false;
	}

//	attach DoFs to vertices
	m_pGrid->attach_to<VertexBase>(m_aIndex);
	m_pGrid->attach_to<EdgeBase>(m_aIndex);
	m_pGrid->attach_to<Face>(m_aIndex);
	m_pGrid->attach_to<Volume>(m_aIndex);

//	access the
	m_aaIndexVRT.access(*m_pGrid, m_aIndex);
	m_aaIndexEDGE.access(*m_pGrid, m_aIndex);
	m_aaIndexFACE.access(*m_pGrid, m_aIndex);
	m_aaIndexVOL.access(*m_pGrid, m_aIndex);

//	done
	return true;
}

////////////////////////////////////////////////////////////////////////////////
//	DoFDistribution
////////////////////////////////////////////////////////////////////////////////

DoFDistribution::
DoFDistribution(GeometricObjectCollection goc,
                ISubsetHandler& sh, storage_manager_type& sm,
                FunctionPattern& fp)
: base_type(goc, fp), m_pISubsetHandler(&sh), m_pStorageManager(&sm),
  m_raaIndexVRT(sm.vertex_attachment_accessor()),
  m_raaIndexEDGE(sm.edge_attachment_accessor()),
  m_raaIndexFACE(sm.face_attachment_accessor()),
  m_raaIndexVOL(sm.volume_attachment_accessor()),
  m_numIndex(0), m_sizeIndexSet(0), m_bGrouped(false)
{
	m_vNumIndex.clear();
	m_vNumIndex.resize(num_subsets(), 0);

// 	Attach indices
	if(!m_pStorageManager->update_attachments())
		throw(UGFatalError("Attachment missing in DoF Storage Manager."));

// 	create offsets
	create_offsets();
}

DoFDistribution::
DoFDistribution(GeometricObjectCollection goc,
                ISubsetHandler& sh, storage_manager_type& sm,
                FunctionPattern& fp,
                const SurfaceView& surfView)
: base_type(goc, fp, surfView), m_pISubsetHandler(&sh), m_pStorageManager(&sm),
  m_raaIndexVRT(sm.vertex_attachment_accessor()),
  m_raaIndexEDGE(sm.edge_attachment_accessor()),
  m_raaIndexFACE(sm.face_attachment_accessor()),
  m_raaIndexVOL(sm.volume_attachment_accessor()),
  m_numIndex(0), m_sizeIndexSet(0), m_bGrouped(false)
{
	m_vNumIndex.clear();
	m_vNumIndex.resize(this->num_subsets(), 0);

// 	Attach indices
	if(!m_pStorageManager->update_attachments())
		throw(UGFatalError("Attachment missing in DoF Storage Manager."));

// 	create offsets
	create_offsets();
}

bool DoFDistribution::has_indices_on(GeometricBaseObject gbo) const
{
//	return if at least on dof on that type
	return (m_vMaxDoFsInDim[gbo] > 0);
}

bool DoFDistribution::has_indices_on(ReferenceObjectID roid) const
{
//	return if at least on dof on that type
	return (m_vMaxDoFsOnROID[roid] > 0);
}

void DoFDistribution::create_offsets(ReferenceObjectID roid)
{
// 	loop subsets
	for(int si = 0; si < num_subsets(); ++si)
	{
	// 	counter
		size_t count = 0;

	//	get dimension of subset
		int dim = dim_subset(si);

	//	loop functions
		for(size_t fct = 0; fct < num_fct(); ++fct)
		{
		//	reset to not defined (initially)
			m_vvvOffsets[roid][si][fct] = NOT_DEF_ON_SUBSET;

		//	if function is not defined, we leave the offset as invalid.
			if(!is_def_in_subset(fct, si))	continue;

		//	get local shape function id
			LFEID lfeID = local_finite_element_id(fct);

		//	get trial space
			const CommonLocalDoFSet& clds = LocalDoFSetProvider::get(lfeID, dim);

		//	get number of DoFs on the reference element need for the space
			const int numDoF = clds.num_dof(roid);

		//	check that numDoFs specified by this roid
			if(numDoF == -1) continue;

		//	set offset for each function defined in the subset
			m_vvvOffsets[roid][si][fct] = count;

		//	increase number of DoFs per type on this subset
			m_vvNumDoFsOnROID[roid][si] += numDoF;

		//	increase the count
			count += numDoF;
		}
	}
}

void DoFDistribution::create_offsets()
{
//	loop all reference element to resize the arrays
	for(int roid=ROID_VERTEX; roid < NUM_REFERENCE_OBJECTS; ++roid)
	{
	//	clear offsets
		m_vvvOffsets[roid].clear();
		m_vvNumDoFsOnROID[roid].clear();

	//	resize for all subsets
		m_vvvOffsets[roid].resize(num_subsets());
		m_vvNumDoFsOnROID[roid].resize(num_subsets(), 0);

	//	resize for each function
		for(int si = 0; si < num_subsets(); ++si)
			m_vvvOffsets[roid][si].resize(num_fct());
	}

//	loop all reference element, but not vertices (no disc there)
	for(int roid=ROID_VERTEX; roid < NUM_REFERENCE_OBJECTS; ++roid)
		create_offsets((ReferenceObjectID) roid);

//	reset dimension maximum
	for(size_t d=0; d <= 3; ++d) m_vMaxDoFsInDim[d] = 0;

//	get max number of dofs per roid
	for(int roid=ROID_VERTEX; roid < NUM_REFERENCE_OBJECTS; ++roid)
	{
	// 	lets find out the maximum number for a type on all subsets
		m_vMaxDoFsOnROID[roid] = 0;
		for(int si = 0; si < num_subsets(); ++si)
		{
			m_vMaxDoFsOnROID[roid] = std::max(m_vMaxDoFsOnROID[roid],
			                                  m_vvNumDoFsOnROID[roid][si]);
		}

	//	lets find out the maximum number of DoFs for objects in dimension
		const int d = ReferenceElementDimension((ReferenceObjectID)roid);
		m_vMaxDoFsInDim[d] = std::max(m_vMaxDoFsInDim[d],
		                              m_vMaxDoFsOnROID[roid]);
	}
}

size_t DoFDistribution::get_free_index(size_t si, ReferenceObjectID roid)
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
		if(!m_bGrouped) m_sizeIndexSet += m_vvNumDoFsOnROID[roid][si];
		else ++m_sizeIndexSet;
	}

//	adjust counters
	if(!m_bGrouped)
	{
		m_numIndex += m_vvNumDoFsOnROID[roid][si];
		m_vNumIndex[si] += m_vvNumDoFsOnROID[roid][si];
	}
	else
	{
		++ m_numIndex;
		++ (m_vNumIndex[si]);
	}

//	return new index
	return freeIndex;
}

void DoFDistribution::push_free_index(size_t freeIndex, size_t si, ReferenceObjectID roid)
{
//	remember index
	m_vFreeIndex.push_back(freeIndex);

//	decrease number of distributed indices
	if(!m_bGrouped)
	{
		m_numIndex -= m_vvNumDoFsOnROID[roid][si];
		m_vNumIndex[si] -= m_vvNumDoFsOnROID[roid][si];
	}
	else
	{
		-- m_numIndex;
		-- (m_vNumIndex[si]);
	}
}

template <typename TElem>
bool DoFDistribution::distribute_indices()
{
//	iterator type
	typedef typename geometry_traits<TElem>::iterator iterator;

//	get reference element type
	typedef typename reference_element_traits<TElem>::reference_element_type
			reference_element_type;

//	get reference object id
	ReferenceObjectID roid = reference_element_type::REFERENCE_OBJECT_ID;

// 	loop subsets
	for(int si = 0; si < num_subsets(); ++si)
	{
	// 	skip if no dofs to be distributed
		if(m_vvNumDoFsOnROID[roid][si] == 0) continue;

	//	get iterators of elems
		iterator iter = this->template begin<TElem>(si);
		iterator iterEnd =  this->template end<TElem>(si);

	// 	loop elems
		for(; iter != iterEnd; ++iter)
		{
		// 	get vertex
			TElem* elem = *iter;

		//	skip shadows
			if(this->m_pSurfaceView != NULL)
				if(this->m_pSurfaceView->is_shadow(elem))
					continue;

		// 	write index
			first_index(elem) = m_numIndex;

		//	increase number of DoFs and DoFs on subset
			if(!m_bGrouped)
			{
				m_numIndex += m_vvNumDoFsOnROID[roid][si];
				m_vNumIndex[si] += m_vvNumDoFsOnROID[roid][si];
			}
			else
			{
				++ m_numIndex;
				++ (m_vNumIndex[si]);
			}
		}
	}

//	post process shadows.
	if(this->m_pSurfaceView != NULL)
	{
	//	loop subsets
		for(int si = 0; si < num_subsets(); ++si)
		{
		//	get iterators of elems
			iterator iter = this->template begin<TElem>(si);
			iterator iterEnd =  this->template end<TElem>(si);

		//	loop elems
			for(; iter != iterEnd; ++iter)
			{
			//	get elem
				TElem* elem = *iter;

			//	skip non-shadows
				if(!this->m_pSurfaceView->is_shadow(elem)) continue;

			//	get child
				TElem* vrtChild = this->m_pSurfaceView->get_shadow_child(elem);

			//	get indices of child
				const size_t indexChild = first_index(vrtChild);

			//	set index of shadow to index of child
				first_index(elem) = indexChild;
			}
		}
	}

//	done
	return true;
}

bool DoFDistribution::distribute_indices()
{
//	storage manage required
	if(m_pStorageManager == NULL)
	{
		UG_LOG("In 'DoFDistribution::distribute_indices:"
				"Storage Manager not set. Aborting.\n");
		return false;
	}

//	function pattern required
	if(this->m_pFuncPattern == NULL)
	{
		UG_LOG("In 'DoFDistribution::distribute_indices:"
				"Function Pattern not set. Aborting.\n");
		return false;
	}

// 	Attach indices
	if(!m_pStorageManager->update_attachments())
	{
		UG_LOG("In 'DoFDistribution::distribute_indices:"
				"Cannot update index attachments. Aborting.\n");
		return false;
	}

// 	create offsets
	create_offsets();

// 	reset counter for all dofs
	m_numIndex = 0;

// 	reset number of dofs
	m_vNumIndex.clear(); m_vNumIndex.resize(num_subsets(), 0);

//	add dofs on all elems
	bool bSuccess = true;
	bSuccess &= distribute_indices<Vertex>();
	bSuccess &= distribute_indices<Edge>();
	// \todo: more....

//	the size of the index set is the number of DoFs
	m_sizeIndexSet = m_numIndex;

//	we're done
	return bSuccess;
}

template <typename TElem>
bool DoFDistribution::permute_indices(std::vector<size_t>& vIndNew)
{
//	get reference element type
	typedef typename reference_element_traits<TElem>::reference_element_type
			reference_element_type;

//	get reference object id
	ReferenceObjectID roid = reference_element_type::REFERENCE_OBJECT_ID;

// 	loop subsets
	for(int si = 0; si < num_subsets(); ++si)
	{
	//	get iterators
		typename geometry_traits<TElem>::iterator iter, iterBegin, iterEnd;
		iterBegin = this->template begin<TElem>(si);
		iterEnd =  this->template end<TElem>(si);

	// 	loop Vertices
		for(iter = iterBegin; iter != iterEnd; ++iter)
		{
		// 	get vertex
			TElem* elem = *iter;

		// 	get current (old) index
			const size_t oldIndex = first_index(elem);

		//	replace old index by new one
			UG_ASSERT(oldIndex < vIndNew.size(), "Index ("<<oldIndex<<") not "
			          	  	  	  	  "in vector of size "<<vIndNew.size());
			first_index(elem) = vIndNew[oldIndex];

		//	set correct permutation in permuting vector
			if(!m_bGrouped)
				for(size_t fct = 0; fct < m_vvNumDoFsOnROID[roid][si]; ++fct)
					vIndNew[oldIndex+fct] = vIndNew[oldIndex] + fct;
		}
	}

	return true;
}

bool DoFDistribution::permute_indices(std::vector<size_t>& vIndNew)
{
//	storage manager required
	if(m_pStorageManager == NULL)
	{
		UG_LOG("ERROR in 'DoFDistribution::permute_indices':"
				" No Storage Manager");
		return false;
	}

//	subsethandler required
	if(m_pISubsetHandler == NULL)
	{
		UG_LOG("ERROR in 'DoFDistribution::permute_indices':"
				" No Subset Handler");
		return false;
	}

//	check, that passed index fields have the same size
	if(this->num_indices() != vIndNew.size())
	{
		UG_LOG("ERROR in 'DoFDistribution::permute_indices': New index set"
				" must have same cardinality for swap indices.\n");
		return false;
	}

	permute_indices<Vertex>(vIndNew);
	permute_indices<Edge>(vIndNew);
	//\todo: more

//	we're done
	return true;
}

template <typename TElem, typename TBaseElem>
bool DoFDistribution::
add_connections_for_adjacent(std::vector<std::vector<size_t> >& vvConnection,
                                     TElem* elem,
                                     std::vector<TBaseElem*> vConnElem)
{
//	get index
	const size_t index = first_index(elem);

//	Get connected indices
	for(size_t i = 0; i < vConnElem.size(); ++i)
	{
	//	Get Vertices of adjacent edges
		TBaseElem* adjElem = vConnElem[i];

	//	get index
		const int si1 = m_pISubsetHandler->get_subset_index(adjElem);

	//	skip iff in no subset
	//  This can happen, when no subset is set from the beginning or
	//  even when a surface grid is considered with vertical
	//	copy nodes.
		if(si1 < 0) continue;

	//	get adjacent index
		const size_t adjInd = first_index(adjElem);

		std::vector<size_t>::iterator it;

	//	add connection (index->adjInd)
		it = find(vvConnection[index].begin(), vvConnection[index].end(), adjInd);
		if(it == vvConnection[index].end()) vvConnection[index].push_back(adjInd);

	//	add the opposite direction (adjInd -> index)
		it = find(vvConnection[adjInd].begin(), vvConnection[adjInd].end(), index);
		if(it == vvConnection[adjInd].end()) vvConnection[adjInd].push_back(index);
	}

	return true;
}


template <typename TBaseElem>
bool DoFDistribution::
get_connections(std::vector<std::vector<size_t> >& vvConnection)
{
//	Adjacent geometric objects
	std::vector<VertexBase*> vVrts;
	std::vector<EdgeBase*> vEdges;
	std::vector<Face*> vFaces;
	std::vector<Volume*> vVols;

// 	Grid
	Grid* grid = m_pStorageManager->get_assigned_grid();

// 	Iterators
	typename geometry_traits<TBaseElem>::iterator iter, iterBegin, iterEnd;

//	success flag
	bool bRet = true;

//	Loop vertices
	for(int si = 0; si < num_subsets(); ++si)
	{
	//	iterators
		iterBegin = this->template begin<TBaseElem>(si);
		iterEnd =  this->template end<TBaseElem>(si);

	//	loop vertices
		for(iter = iterBegin; iter != iterEnd; ++iter)
		{
		// 	Get vertex
			TBaseElem* elem = *iter;

		//	skip shadows
			if(m_pSurfaceView)
				if(m_pSurfaceView->is_shadow(elem))
					continue;

		//	Get connected elements
			if(m_vMaxDoFsInDim[VERTEX] > 0) CollectVertices(vVrts, *grid, elem);
			if(m_vMaxDoFsInDim[EDGE] > 0)	CollectEdges(vEdges, *grid, elem);
			if(m_vMaxDoFsInDim[FACE] > 0)	CollectFaces(vFaces, *grid, elem);
			if(m_vMaxDoFsInDim[VOLUME] > 0) CollectVolumes(vVols, *grid, elem);

		//	Get connections via shadow
			if(m_pSurfaceView)
			{
			//	skip if not shadowing
				if(m_pSurfaceView->shadows(elem))
				{
				//	get parent
					TBaseElem* elemParent =
						dynamic_cast<TBaseElem*>(m_pSurfaceView->get_parent(elem));

				//	Get connected elements
					if(elemParent != NULL){
						if(m_vMaxDoFsInDim[VERTEX] > 0) CollectVertices(vVrts, *grid, elemParent, false);
						if(m_vMaxDoFsInDim[EDGE] > 0) CollectEdges(vEdges, *grid, elemParent, false);
						if(m_vMaxDoFsInDim[FACE] > 0) CollectFaces(vFaces, *grid, elemParent, false);
						if(m_vMaxDoFsInDim[VOLUME] > 0) CollectVolumes(vVols, *grid, elemParent, false);
					}
				}
			}

		//	add indices of connected elems to the connection list
			if(m_vMaxDoFsInDim[VERTEX] > 0)
				bRet &= add_connections_for_adjacent<TBaseElem, VertexBase>(vvConnection, elem, vVrts);
			if(m_vMaxDoFsInDim[EDGE] > 0)
				bRet &= add_connections_for_adjacent<TBaseElem, EdgeBase>(vvConnection, elem, vEdges);
			if(m_vMaxDoFsInDim[FACE] > 0)
				bRet &= add_connections_for_adjacent<TBaseElem, Face>(vvConnection, elem, vFaces);
			if(m_vMaxDoFsInDim[VOLUME] > 0)
				bRet &= add_connections_for_adjacent<TBaseElem, Volume>(vvConnection, elem, vVols);
		}
	}

	return true;
}

bool DoFDistribution::get_connections(std::vector<std::vector<size_t> >& vvConnection)
{
//	storage manager required
	if(m_pStorageManager == NULL)
	{
		UG_LOG("ERROR in 'DoFDistribution::get_connections':"
				" No Storage Manager");
		return false;
	}

//	subset handler required
	if(m_pISubsetHandler == NULL)
	{
		UG_LOG("ERROR in 'DoFDistribution::get_connections':"
				" No Subset Handler");
		return false;
	}

//	check that in all subsets same number of functions and at least one
//	if this is not the case for non-grouped DoFs, we cannot allow reordering
	if(!m_bGrouped)
	{
		size_t numDoFs = 0;
		for(int si = 0; si < num_subsets(); ++si)
			for(int roid = 0; roid < NUM_REFERENCE_OBJECTS; ++roid)
		{
			if(m_vvNumDoFsOnROID[roid][si] == 0) continue;
			if(numDoFs == 0) {numDoFs = m_vvNumDoFsOnROID[roid][si]; continue;}
			if(m_vvNumDoFsOnROID[roid][si] != numDoFs)
			{
				UG_LOG("ERROR in 'DoFDistribution::get_connections':"
						" Currently only implemented iff same number of DoFs on"
						" all geometric objects in all subsets.\n");
				return false;
			}
		}
	}

//	clear neighbors
	vvConnection.clear(); vvConnection.resize(m_numIndex);

//	if no subset given, we're done
	if(num_subsets() == 0) return true;

//	success flag
	bool bRet = true;

	if(m_vMaxDoFsInDim[VERTEX] > 0)
		bRet &= get_connections<VertexBase>(vvConnection);
	if(m_vMaxDoFsInDim[EDGE] > 0)
		bRet &= get_connections<EdgeBase>(vvConnection);
	if(m_vMaxDoFsInDim[FACE] > 0)
		bRet &= get_connections<Face>(vvConnection);
	if(m_vMaxDoFsInDim[VOLUME] > 0)
		bRet &= get_connections<Volume>(vvConnection);

//	we're done
	return bRet;
}


template <typename TElem>
bool DoFDistribution::defragment(std::vector<std::pair<size_t, size_t> >& vReplaced)
{
//	get reference element type
	typedef typename reference_element_traits<TElem>::reference_element_type
			reference_element_type;

//	get reference object id
	ReferenceObjectID roid = reference_element_type::REFERENCE_OBJECT_ID;

//	loop subsets
	for(int si = 0; si < num_subsets(); ++si)
	{
	// 	skip if no DoFs to be distributed
		if(m_vvNumDoFsOnROID[roid][si]) continue;

		typename geometry_traits<TElem>::iterator iter, iterBegin, iterEnd;
		iterBegin = this->begin<TElem>(si);
		iterEnd =  this->end<TElem>(si);

	// 	loop Vertices
		for(iter = iterBegin; iter != iterEnd; ++iter)
		{
		// 	get vertex
			TElem* elem = *iter;

		//	get old (current) index
			const size_t oldIndex = first_index(elem);

		// 	check if index must be replaced by lower one
			if(first_index(elem) < m_numIndex) continue;

		//	get free index
			const size_t newIndex = get_free_index(si, roid);

		//	remember replacement
			vReplaced.push_back(std::pair<size_t,size_t>(oldIndex, newIndex));

		//	overwrite index
			first_index(elem) = newIndex;

		//	adjust counters, since this was an replacement
			if(!m_bGrouped)
			{
				m_numIndex -= m_vvNumDoFsOnROID[roid][si];
				m_sizeIndexSet -= m_vvNumDoFsOnROID[roid][si];
				m_vNumIndex[si] -= m_vvNumDoFsOnROID[roid][si];
			}
			else
			{
				--m_numIndex;
				--m_sizeIndexSet;
				--(m_vNumIndex[si]);
			}
		}
	}

//	done
	return true;
}

bool DoFDistribution::defragment()
{
//	we loop all indices and those with highest degree are replaced by a free
//	index. This operation is performed in O(number of Indices)
//	All indices with an index >= m_numDoFs are replaced.

//	check, if holes exist. If not, we're done
	if(m_vFreeIndex.empty())
	{
	//	the index set might have been increased. Thus resize grid functions
		num_indices_changed(m_numIndex);

	//	we're done
		return true;
	}

//	pairs of replaced indices
	std::vector<std::pair<size_t, size_t> > vReplaced;

//	success flag
	bool bSuc = true;

//	replace for all element type
	bSuc &= defragment<Vertex>(vReplaced);
	//\todo: more ...

//	check success
	if(!bSuc)
	{
		UG_LOG("ERROR in 'DoFDistribution::defragment': Something wrong during"
				" defragment.\n");
		return false;
	}

//	check that all holes have been removed
	if(m_numIndex != m_sizeIndexSet)
	{
		UG_LOG("ERROR in 'DoFDistribution::defragment': Still holes in index "
				"set after compression. Check implementation.\n");
		return false;
	}

//	copy values (from back into holes) for managed grid functions
	if(!indices_swaped(vReplaced, true)) return false;

//	cut of unused tail of managed grid functions
	num_indices_changed(m_numIndex);

//	we're done
	return true;
}

size_t DoFDistribution::
inner_multi_indices(multi_index_vector_type& ind,
                    const size_t firstIndex, const int si, const size_t fct,
                    const ReferenceObjectID roid) const
{
//	check that function is def on subset
	if(!is_def_in_subset(fct, si)) return ind.size();

//	only in case of vertex, we have DoFs
	if(m_vvNumDoFsOnROID[roid][si] > 0)
	{
	//	get local shape function id
		LFEID lsfsID = local_finite_element_id(fct);

	//	get trial space
		const ILocalDoFSet& lsfs = LocalDoFSetProvider::get(lsfsID, roid);

	//	get number of DoFs in this sub-geometric object
		const size_t numDoFsOnSub = lsfs.num_dof(roid);

		if(!m_bGrouped)
		{
		//	compute index
			const size_t index = firstIndex + m_vvvOffsets[roid][si][fct];

		//	add dof to local indices
			for(size_t j = 0; j < numDoFsOnSub; ++j)
				ind.push_back(index_type(index+j,0));
		}
		else
		{
		//	compute index
			const size_t comp = m_vvvOffsets[roid][si][fct];

		//	add dof to local indices
			for(size_t j = 0; j < numDoFsOnSub; ++j)
				ind.push_back(index_type(firstIndex, comp+j));
		}
	}

//	return number of indices
	return ind.size();
}

size_t
DoFDistribution::inner_algebra_indices(algebra_index_vector_type& ind,
                                       const size_t firstIndex, const int si,
                                       const ReferenceObjectID roid) const
{
//	only in case of vertex, we have DoFs
	if(m_vvNumDoFsOnROID[roid][si] > 0)
	{
		if(!m_bGrouped)
		{
			for(size_t fct = 0; fct < num_fct(); ++fct)
			{
			//	check that function is def on subset
				if(!is_def_in_subset(fct, si)) continue;

			//	get local shape function id
				LFEID lsfsID = local_finite_element_id(fct);

			//	get trial space
				const ILocalDoFSet& lsfs = LocalDoFSetProvider::get(lsfsID, roid);

			//	get number of DoFs in this sub-geometric object
				const size_t numDoFsOnSub = lsfs.num_dof(roid);

			//	compute index
				const size_t index = firstIndex + m_vvvOffsets[roid][si][fct];

			//	add dof to local indices
				for(size_t j = 0; j < numDoFsOnSub; ++j)
					ind.push_back(index+j);
			}
		}
		else
		{
		//	add dof to local indices
			ind.push_back(firstIndex);
		}
	}

//	return number of indices
	return ind.size();
}


template <typename TBaseElem>
void DoFDistribution::grid_obj_added(TBaseElem* elem)
{
//	subset handler required
	if(m_pISubsetHandler == NULL)
	{
		UG_LOG("ERROR in 'DoFDistribution::grid_obj_added':"
				" No Subset Handler");
		throw(UGFatalError("Subset Handler missng."));
	}

//	Only for the surface dof manager:
//	If the added elem is a shadow, use index of child (shadowing) elem
	if(m_pSurfaceView && m_pSurfaceView->is_shadow(elem))
	{
	//	get child
		TBaseElem* elemChild = m_pSurfaceView->get_shadow_child<TBaseElem>(elem);

	//	get indices of child
		const size_t indexChild = first_index(elemChild);

	//	set index of shadow to index of child
		first_index(elem) = indexChild;
	}
//	normally, we can set a new free index
	else
	{
	//	get subset index
		const int si = m_pISubsetHandler->get_subset_index(elem);

	//	get type
		ReferenceObjectID roid = (ReferenceObjectID) elem->reference_object_id();

	// 	write next free index
		first_index(elem) = get_free_index(si, roid);
	}
}

template <typename TBaseElem>
void DoFDistribution::grid_obj_to_be_removed(TBaseElem* elem)
{
//	get subset index
	const int si = m_pISubsetHandler->get_subset_index(elem);

//	get type
	ReferenceObjectID roid = (ReferenceObjectID) elem->reference_object_id();

// 	remember free index
	push_free_index(first_index(elem), si, roid);
}

template <typename TBaseElem>
void DoFDistribution::
grid_obj_replaced(TBaseElem* elemNew, TBaseElem* elemOld)
{
//	check subset index
	UG_ASSERT(m_pISubsetHandler->get_subset_index(elemNew) ==
			  m_pISubsetHandler->get_subset_index(elemOld),
	          "New elem does not have same subset as replaced on.");

//	copy index
	first_index(elemNew) = first_index(elemOld);
}


void DoFDistribution::grid_obj_added(VertexBase* vrt)
{
	grid_obj_added<VertexBase>(vrt);
}
void DoFDistribution::grid_obj_added(EdgeBase* edge)
{
	grid_obj_added<EdgeBase>(edge);
}
void DoFDistribution::grid_obj_added(Face* face)
{
	grid_obj_added<Face>(face);
}
void DoFDistribution::grid_obj_added(Volume* vol)
{
	grid_obj_added<Volume>(vol);
}

void DoFDistribution::grid_obj_to_be_removed(VertexBase* vrt)
{
	grid_obj_to_be_removed<VertexBase>(vrt);
}
void DoFDistribution::grid_obj_to_be_removed(EdgeBase* edge)
{
	grid_obj_to_be_removed<EdgeBase>(edge);
}
void DoFDistribution::grid_obj_to_be_removed(Face* face)
{
	grid_obj_to_be_removed<Face>(face);
}
void DoFDistribution::grid_obj_to_be_removed(Volume* vol)
{
	grid_obj_to_be_removed<Volume>(vol);
}

void DoFDistribution::grid_obj_replaced(VertexBase* vrtNew, VertexBase* vrtOld)
{
	grid_obj_replaced<VertexBase>(vrtNew, vrtOld);
}
void DoFDistribution::grid_obj_replaced(EdgeBase* edgeNew, EdgeBase* edgeOld)
{
	grid_obj_replaced<EdgeBase>(edgeNew, edgeOld);
}
void DoFDistribution::grid_obj_replaced(Face* faceNew, Face* faceOld)
{
	grid_obj_replaced<Face>(faceNew, faceOld);
}
void DoFDistribution::grid_obj_replaced(Volume* volNew, Volume* volOld)
{
	grid_obj_replaced<Volume>(volNew, volOld);
}


} // end namespace ug

