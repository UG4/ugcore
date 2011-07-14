
//	standard included
#include <vector>
#include <queue>
#include <algorithm>

//	header
#include "conform.h"

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
  m_numIndex(0), m_sizeIndexSet(0)
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
  m_numIndex(0), m_sizeIndexSet(0)
{
	m_vNumIndex.clear();
	m_vNumIndex.resize(this->num_subsets(), 0);

// 	Attach indices
	if(!m_pStorageManager->update_attachments())
		throw(UGFatalError("Attachment missing in DoF Storage Manager."));

// 	create offsets
	create_offsets();
}

bool DoFDistribution::has_dofs_on(int dim) const
{
//	return if at least on dof on that type
	return (m_vMaxDoFsInDim[dim] > 0);
}

template <typename TRefElem>
void DoFDistribution::create_offsets_of_type()
{
	ReferenceObjectID type = TRefElem::REFERENCE_OBJECT_ID;

//	clear offsets
	m_vvvOffsets[type].clear();
	m_vvNumDoFsOnType[type].clear();

//	resize for all subsets
	m_vvvOffsets[type].resize(num_subsets());
	m_vvNumDoFsOnType[type].resize(num_subsets(), 0);

// 	loop subsets
	for(int si = 0; si < num_subsets(); ++si)
	{
	// 	counter
		size_t count = 0;

	//	resize for each function
		m_vvvOffsets[type][si].resize(num_fct());

	//	loop functions
		for(size_t fct = 0; fct < num_fct(); ++fct)
		{
		//	if function is not defined, we set the offset as invalid.
			if(!is_def_in_subset(fct, si))
				m_vvvOffsets[type][si][fct] = NOT_DEF_ON_SUBSET;

		//	get local shape function id
			LSFSID lsfsID = local_shape_function_set_id(fct);

		//	get trial space
			const LocalShapeFunctionSet<TRefElem>& lsfs =
						LocalShapeFunctionSetProvider::get<TRefElem>(lsfsID);

		//	get number of DoFs on the reference element need for the space
			const size_t numDoF = lsfs.num_sh(type);

		//	set offset for each function defined in the subset
			m_vvvOffsets[type][si][fct] = count;

		//	increase number of DoFs per type on this subset
			m_vvNumDoFsOnType[type][si] += numDoF;

		//	increase the count
			count += numDoF;
		}
	}

// 	lets find out the maximum number for a type on all subsets
	m_vMaxDoFsOnType[type] = 0;
	for(int si = 0; si < num_subsets(); ++si)
	{
		m_vMaxDoFsOnType[type] = std::max(m_vMaxDoFsOnType[type],
		                                  m_vvNumDoFsOnType[type][si]);
	}

//	lets find out the maximum number of DoFs for objects in dimension
	const int d = TRefElem::dim;
	m_vMaxDoFsInDim[d] = std::max(m_vMaxDoFsInDim[d],
	                              m_vMaxDoFsOnType[type]);
}

void DoFDistribution::create_offsets()
{
//	reset dimension maximum
	for(size_t d=0; d <= 3; ++d) m_vMaxDoFsInDim[d] = 0;

	create_offsets_of_type<ReferenceVertex>();
	create_offsets_of_type<ReferenceEdge>();
	create_offsets_of_type<ReferenceTriangle>();
	create_offsets_of_type<ReferenceQuadrilateral>();
	create_offsets_of_type<ReferenceTetrahedron>();
	create_offsets_of_type<ReferencePrism>();
	create_offsets_of_type<ReferencePyramid>();
	create_offsets_of_type<ReferenceHexahedron>();
}

size_t DoFDistribution::get_free_index(size_t si, ReferenceObjectID type)
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
		if(!m_bGrouped) m_sizeIndexSet += m_vvNumDoFsOnType[type][si];
		else ++m_sizeIndexSet;
	}

//	adjust counters
	if(!m_bGrouped)
	{
		m_numIndex += m_vvNumDoFsOnType[type][si];
		m_vNumIndex[si] += m_vvNumDoFsOnType[type][si];
	}
	else
	{
		++ m_numIndex;
		++ (m_vNumIndex[si]);
	}

//	return new index
	return freeIndex;
}

void DoFDistribution::push_free_index(size_t freeIndex, size_t si, ReferenceObjectID type)
{
//	remember index
	m_vFreeIndex.push_back(freeIndex);

//	decrease number of distributed indices
	if(!m_bGrouped)
	{
		m_numIndex -= m_vvNumDoFsOnType[type][si];
		m_vNumIndex[si] -= m_vvNumDoFsOnType[type][si];
	}
	else
	{
		-- m_numIndex;
		-- (m_vNumIndex[si]);
	}
}

template <typename TElem>
bool DoFDistribution::distribute_dofs()
{
//	iterator type
	typedef typename geometry_traits<TElem>::iterator iterator;

//	get reference element type
	typedef typename reference_element_traits<TElem>::reference_element_type
			reference_element_type;

//	get reference object id
	ReferenceObjectID type = reference_element_type::REFERENCE_OBJECT_ID;

// 	loop subsets
	for(int si = 0; si < num_subsets(); ++si)
	{
	// 	skip if no dofs to be distributed
		if(m_vvNumDoFsOnType[type][si] == 0) continue;

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
				m_numIndex += m_vvNumDoFsOnType[type][si];
				m_vNumIndex[si] += m_vvNumDoFsOnType[type][si];
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

bool DoFDistribution::distribute_dofs()
{
//	storage manage required
	if(m_pStorageManager == NULL)
	{
		UG_LOG("In 'DoFDistribution::distribute_dofs:"
				"Storage Manager not set. Aborting.\n");
		return false;
	}

//	function pattern required
	if(this->m_pFuncPattern == NULL)
	{
		UG_LOG("In 'DoFDistribution::distribute_dofs:"
				"Function Pattern not set. Aborting.\n");
		return false;
	}

// 	Attach indices
	if(!m_pStorageManager->update_attachments()) return false;

// 	create offsets
	create_offsets();

// 	reset counter for all dofs
	m_numIndex = 0;

// 	reset number of dofs
	m_vNumIndex.clear(); m_vNumIndex.resize(num_subsets(), 0);

//	add dofs on all elems
	bool bSuccess = true;
	bSuccess &= distribute_dofs<Vertex>();
	bSuccess &= distribute_dofs<Edge>();
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
	ReferenceObjectID type = reference_element_type::REFERENCE_OBJECT_ID;

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
			first_index(elem) = vIndNew[oldIndex];

		//	set correct permutation in permuting vector
			if(!m_bGrouped)
				for(size_t fct = 0; fct < m_vvNumDoFsOnType[type][si]; ++fct)
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
	if(this->num_dofs() != vIndNew.size())
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
			CollectVertices(vVrts, *grid, elem);
			CollectEdges(vEdges, *grid, elem);
			CollectFaces(vFaces, *grid, elem);
			CollectVolumes(vVols, *grid, elem);

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
						CollectVertices(vVrts, *grid, elemParent, false);
						CollectEdges(vEdges, *grid, elemParent, false);
						CollectFaces(vFaces, *grid, elemParent, false);
						CollectVolumes(vVols, *grid, elemParent, false);
					}
				}
			}

		//	add indices of connected elems to the connection list
			bRet &= add_connections_for_adjacent<TBaseElem, VertexBase>(vvConnection, elem, vVrts);
			bRet &= add_connections_for_adjacent<TBaseElem, EdgeBase>(vvConnection, elem, vEdges);
			bRet &= add_connections_for_adjacent<TBaseElem, Face>(vvConnection, elem, vFaces);
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

//	check that in all subsets same number of functions and at least one
//	if this is not the case for ungrouped DoFs, we cannot allow reordering
	if(!m_bGrouped)
	{
		size_t numDoFs = 0;
		for(int si = 0; si < num_subsets(); ++si)
			for(int type = 0; type < NUM_REFERENCE_OBJECTS; ++type)
		{
			if(m_vvNumDoFsOnType[type][si] == 0) continue;
			if(numDoFs == 0) {numDoFs = m_vvNumDoFsOnType[type][si]; continue;}
			if(m_vvNumDoFsOnType[type][si] != numDoFs)
			{
				UG_LOG("ERROR in 'P1DoFDistribution::get_connections':"
						" Currently only implemented iff same number of DoFs on"
						" all geometric objects in all subsets.\n");
				return false;
			}
		}
	}

//	clear neighbours
	vvConnection.clear(); vvConnection.resize(m_numIndex);

//	if no subset given, we're done
	if(num_subsets() == 0) return true;

//	success flag
	bool bRet = true;

	bRet &= get_connections<VertexBase>(vvConnection);
	bRet &= get_connections<EdgeBase>(vvConnection);
	bRet &= get_connections<Face>(vvConnection);
	bRet &= get_connections<Volume>(vvConnection);

//	we're done
	return bRet;
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
		ReferenceObjectID type = elem->reference_object_id();

	// 	write next free index
		first_index(elem) = get_free_index(si, type);
	}
}

template <typename TBaseElem>
void DoFDistribution::grid_obj_to_be_removed(TBaseElem* elem)
{
//	get subset index
	const int si = m_pISubsetHandler->get_subset_index(elem);

//	get type
	ReferenceObjectID type = elem->reference_object_id();

// 	remember free index
	push_free_index(first_index(elem), si, type);
}

template <typename TBaseElem>
void DoFDistribution::
grid_obj_replaced(TBaseElem* elemNew, TBaseElem* elemOld)
{
//	get subset index
	const int si = m_pISubsetHandler->get_subset_index(elemOld);

	UG_ASSERT(m_pISubsetHandler->get_subset_index(elemNew) == si,
	          "New elem does not have same subset as replaced on.");

//	copy index
	first_index(elemNew) = first_index(elemOld);
}

template <typename TElem>
bool DoFDistribution::defragment(std::vector<std::pair<size_t, size_t> >& vReplaced)
{
//	get reference element type
	typedef typename reference_element_traits<TElem>::reference_element_type
			reference_element_type;

//	get reference object id
	ReferenceObjectID type = reference_element_type::REFERENCE_OBJECT_ID;

//	loop subsets
	for(int si = 0; si < num_subsets(); ++si)
	{
	// 	skip if no DoFs to be distributed
		if(m_vvNumDoFsOnType[type][si]) continue;

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
			const size_t newIndex = get_free_index(si, type);

		//	remember replacement
			vReplaced.push_back(std::pair<size_t,size_t>(oldIndex, newIndex));

		//	overwrite index
			first_index(elem) = newIndex;

		//	adjust counters, since this was an replacement
			if(!m_bGrouped)
			{
				m_numIndex -= m_vvNumDoFsOnType[type][si];
				m_sizeIndexSet -= m_vvNumDoFsOnType[type][si];
				m_vNumIndex[si] -= m_vvNumDoFsOnType[type][si];
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
                    const ReferenceObjectID type) const
{
//	check that function is def on subset
	if(!is_def_in_subset(fct, si)) return ind.size();

//	only in case of vertex, we have DoFs
	if(m_vvNumDoFsOnType[type][si] > 0)
	{
	//	get local shape function id
		LSFSID lsfsID = local_shape_function_set_id(fct);

	//	get trial space
		const LocalShapeFunctionSetBase& lsfs =
						LocalShapeFunctionSetProvider::get(lsfsID, type);

	//	get number of DoFs in this sub-geometric object
		const size_t numDoFsOnSub = lsfs.num_sh(type);

		if(!m_bGrouped)
		{
		//	compute index
			const size_t index = firstIndex + m_vvvOffsets[type][si][fct];

		//	add dof to local indices
			for(size_t j = 0; j < numDoFsOnSub; ++j)
				ind.push_back(index_type(index+j,0));
		}
		else
		{
		//	compute index
			const size_t comp = m_vvvOffsets[type][si][fct];

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
                                       const ReferenceObjectID type) const
{
//	only in case of vertex, we have DoFs
	if(m_vvNumDoFsOnType[type][si] > 0)
	{
		if(!m_bGrouped)
		{
			for(size_t fct = 0; fct < num_fct(); ++fct)
			{
			//	check that function is def on subset
				if(!is_def_in_subset(fct, si)) continue;

			//	get local shape function id
				LSFSID lsfsID = local_shape_function_set_id(fct);

			//	get trial space
				const LocalShapeFunctionSetBase& lsfs =
							LocalShapeFunctionSetProvider::get(lsfsID, type);

			//	get number of DoFs in this sub-geometric object
				const size_t numDoFsOnSub = lsfs.num_sh(type);

			//	compute index
				const size_t index = firstIndex + m_vvvOffsets[type][si][fct];

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

} // end namespace ug

