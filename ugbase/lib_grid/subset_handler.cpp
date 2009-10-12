//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m11 d24

#include <cassert>
#include "subset_handler.h"

using namespace std;

namespace ug
{
////////////////////////////////////////////////////////////////////////
//	SubsetInfo implementation
SubsetInfo::SubsetInfo()
{
	name = "";
	materialIndex = -1;
}

////////////////////////////////////////////////////////////////////////
//	SubsetHandler implementation
SubsetHandler::SubsetHandler() : m_aSubsetIndex(false), m_aIterator(false)
{
	m_pGrid = NULL;
	m_defaultSubsetIndex = -1;
	m_bSubsetInheritanceEnabled = true;
}

SubsetHandler::SubsetHandler(Grid& grid) : m_aSubsetIndex(false), m_aIterator(false)
{
	m_pGrid = NULL;
	m_defaultSubsetIndex = -1;
	m_bSubsetInheritanceEnabled = true;
	assign_grid(grid);
}

SubsetHandler::
SubsetHandler(const SubsetHandler& sh) : m_aSubsetIndex(false), m_aIterator(false)
{
	assert(!"WARNING in SubsetHandler::SubsetHandler(const SubsetHandler& sel): Copy-Constructor not yet implemented!");
	LOG("WARNING in SubsetHandler::SubsetHandler(const SubsetHandler& sel): Copy-Constructor not yet implemented! Expect unexpected behaviour!");
}

SubsetHandler::~SubsetHandler()
{
	if(m_pGrid != NULL)
		m_pGrid->unregister_observer(this);
}

void SubsetHandler::assign_grid(Grid& grid)
{
	grid.register_observer(this, OT_GRID_OBSERVER | OT_VERTEX_OBSERVER | OT_EDGE_OBSERVER |
									OT_FACE_OBSERVER | OT_VOLUME_OBSERVER);
}

Grid* SubsetHandler::get_assigned_grid()
{
	return m_pGrid;
}

void SubsetHandler::set_default_subset_index(int subsetIndex)
{
	m_defaultSubsetIndex = subsetIndex;
	if(subsetIndex < 0)
		m_defaultSubsetIndex = -1;
}

void SubsetHandler::enable_subset_inheritance(bool bEnable)
{
	m_bSubsetInheritanceEnabled = bEnable;
}

void SubsetHandler::clear()
{
	for(uint i = 0; i < num_subsets(); ++i)
		delete m_subsets[i];

	m_subsets.clear();
	m_subsetInfos.clear();
	if(m_pGrid)
	{
		reset_subset_indices<VertexBase>(m_pGrid->vertices_begin(), m_pGrid->vertices_end());
		reset_subset_indices<EdgeBase>(m_pGrid->edges_begin(), m_pGrid->edges_end());
		reset_subset_indices<Face>(m_pGrid->faces_begin(), m_pGrid->faces_end());
		reset_subset_indices<Volume>(m_pGrid->volumes_begin(), m_pGrid->volumes_end());
	}
}

void SubsetHandler::clear_subset(int subsetIndex)
{
	assert((subsetIndex >= 0) && (subsetIndex < (int)num_subsets()) &&
			"ERROR in SubsetHandler::clear_subset(): bad subset index.");

	m_subsetInfos[subsetIndex] = SubsetInfo();

	if(m_pGrid)
	{
		reset_subset_indices<VertexBase>(begin<VertexBase>(subsetIndex),
											end<VertexBase>(subsetIndex));
		reset_subset_indices<EdgeBase>(begin<EdgeBase>(subsetIndex),
											end<EdgeBase>(subsetIndex));
		reset_subset_indices<Face>(begin<Face>(subsetIndex),
											end<Face>(subsetIndex));
		reset_subset_indices<Volume>(begin<Volume>(subsetIndex),
											end<Volume>(subsetIndex));

		for(int i = 0; i < NUM_GEOMETRIC_BASE_OBJECTS; ++i)
			m_subsets[subsetIndex]->m_elements[i].clear();
	}
}

void SubsetHandler::clear_subsets()
{
	for(int i = 0; i < (int)num_subsets(); ++i)
		clear_subset(i);
}

void SubsetHandler::resize_subset_vec(int newSize)
{
	while((int)num_subsets() < newSize)
		m_subsets.push_back(new Subset);

	m_subsetInfos.resize(newSize);
}

void SubsetHandler::set_subset_info(int subsetIndex, const SubsetInfo& subsetInfo)
{
	if(subsetIndex >= (int)num_subsets())
		resize_subset_vec(subsetIndex + 1);

	m_subsetInfos[subsetIndex] = subsetInfo;
}

SubsetInfo& SubsetHandler::subset_info(int subsetIndex)
{
	if(subsetIndex >= (int)num_subsets())
		resize_subset_vec(subsetIndex + 1);

	return m_subsetInfos[subsetIndex];
}

const SubsetInfo& SubsetHandler::subset_info(int subsetIndex) const
{
	assert(((subsetIndex >= 0) && (subsetIndex < (int)num_subsets())) && "ERROR in SubsetHandler::subset_info(..) const: bad subset index. Use non-const version to avoid this Problem.");

	return m_subsetInfos[subsetIndex];
}

void SubsetHandler::assign_subset(VertexBase* elem, int subsetIndex)
{
	assert((m_pGrid != NULL) && "ERROR in SubsetHandler::assign_subset(vrt): No grid assigned to SubsetHandler.");

	if(subsetIndex >= (int)num_subsets())
		resize_subset_vec(subsetIndex + 1);

//	check if we have to remove elem from a subset.
	if(m_aaSubsetIndexVRT[elem] != -1)
		m_subsets[m_aaSubsetIndexVRT[elem]]->m_elements[VERTEX].erase(m_aaIteratorVRT[elem], elem->shared_pipe_section());

//	assign the new index
	m_aaSubsetIndexVRT[elem] = subsetIndex;

//	add the element to the subset.
	if(subsetIndex != -1)
		m_aaIteratorVRT[elem] = m_subsets[subsetIndex]->m_elements[VERTEX].insert(elem, elem->shared_pipe_section());
}

void SubsetHandler::assign_subset(EdgeBase* elem, int subsetIndex)
{
	assert((m_pGrid != NULL) && "ERROR in SubsetHandler::assign_subset(edge): No grid assigned to SubsetHandler.");

	if(subsetIndex >= (int)num_subsets())
		resize_subset_vec(subsetIndex + 1);

//	check if we have to remove elem from a subset.
	if(m_aaSubsetIndexEDGE[elem] != -1)
		m_subsets[m_aaSubsetIndexEDGE[elem]]->m_elements[EDGE].erase(m_aaIteratorEDGE[elem], elem->shared_pipe_section());

//	assign the new index
	m_aaSubsetIndexEDGE[elem] = subsetIndex;

//	add the element to the subset.
	if(subsetIndex != -1)
		m_aaIteratorEDGE[elem] = m_subsets[subsetIndex]->m_elements[EDGE].insert(elem, elem->shared_pipe_section());
}

void SubsetHandler::assign_subset(Face* elem, int subsetIndex)
{
	assert((m_pGrid != NULL) && "ERROR in SubsetHandler::assign_subset(face): No grid assigned to SubsetHandler.");

	if(subsetIndex >= (int)num_subsets())
		resize_subset_vec(subsetIndex + 1);

//	check if we have to remove elem from a subset.
	if(m_aaSubsetIndexFACE[elem] != -1)
		m_subsets[m_aaSubsetIndexFACE[elem]]->m_elements[FACE].erase(m_aaIteratorFACE[elem], elem->shared_pipe_section());

//	assign the new index
	m_aaSubsetIndexFACE[elem] = subsetIndex;

//	add the element to the subset.
	if(subsetIndex != -1)
		m_aaIteratorFACE[elem] = m_subsets[subsetIndex]->m_elements[FACE].insert(elem, elem->shared_pipe_section());
}

void SubsetHandler::assign_subset(Volume* elem, int subsetIndex)
{
	assert((m_pGrid != NULL) && "ERROR in SubsetHandler::assign_subset(vol): No grid assigned to SubsetHandler.");

	if(subsetIndex >= (int)num_subsets())
		resize_subset_vec(subsetIndex + 1);

//	check if we have to remove elem from a subset.
	if(m_aaSubsetIndexVOL[elem] != -1)
		m_subsets[m_aaSubsetIndexVOL[elem]]->m_elements[VOLUME].erase(m_aaIteratorVOL[elem], elem->shared_pipe_section());

//	assign the new index
	m_aaSubsetIndexVOL[elem] = subsetIndex;

//	add the element to the subset.
	if(subsetIndex != -1)
		m_aaIteratorVOL[elem] = m_subsets[subsetIndex]->m_elements[VOLUME].insert(elem, elem->shared_pipe_section());
}

int SubsetHandler::get_subset_index(GeometricObject* elem)
{
	uint type = elem->base_object_type_id();
	switch(type)
	{
		case VERTEX:
		{
			VertexBase* pVrt = dynamic_cast<VertexBase*>(elem);
			assert((pVrt != NULL) && "ERROR in SubsetHandler::get_subset_index(GeometricObject* elem): elem type and base_type_id do not match!");
			return get_subset_index(pVrt);
		}
		case EDGE:
		{
			EdgeBase* pEdge = dynamic_cast<EdgeBase*>(elem);
			assert((pEdge != NULL) && "ERROR in SubsetHandler::get_subset_index(GeometricObject* elem): elem type and base_type_id do not match!");
			return get_subset_index(pEdge);
		}
		case FACE:
		{
			Face* pFace = dynamic_cast<Face*>(elem);
			assert((pFace != NULL) && "ERROR in SubsetHandler::get_subset_index(GeometricObject* elem): elem type and base_type_id do not match!");
			return get_subset_index(pFace);
		}
		case VOLUME:
		{
			Volume* pVol = dynamic_cast<Volume*>(elem);
			assert((pVol != NULL) && "ERROR in SubsetHandler::get_subset_index(GeometricObject* elem): elem type and base_type_id do not match!");
			return get_subset_index(pVol);
		}
	}

//	we should never arrive at this point
	assert(!"ERROR in SubsetHandler::get_subset_index(GeometricObject* elem): Program should never reach this point!");
	return -1;
}

int SubsetHandler::get_subset_index(VertexBase* elem)
{
	assert((m_pGrid != NULL) && "ERROR in SubsetHandler::get_subset_index(vrt): No grid assigned to SubsetHandler.");
	return m_aaSubsetIndexVRT[elem];
}

int SubsetHandler::get_subset_index(EdgeBase* elem)
{
	assert((m_pGrid != NULL) && "ERROR in SubsetHandler::get_subset_index(edge): No grid assigned to SubsetHandler.");
	return m_aaSubsetIndexEDGE[elem];
}

int SubsetHandler::get_subset_index(Face* elem)
{
	assert((m_pGrid != NULL) && "ERROR in SubsetHandler::get_subset_index(face): No grid assigned to SubsetHandler.");
	return m_aaSubsetIndexFACE[elem];
}

int SubsetHandler::get_subset_index(Volume* elem)
{
	assert((m_pGrid != NULL) && "ERROR in SubsetHandler::get_subset_index(vol): No grid assigned to SubsetHandler.");
	return m_aaSubsetIndexVOL[elem];
}

GeometricObjectCollection
SubsetHandler::get_geometric_object_collection(int subsetIndex)
{
	assert((subsetIndex >= 0) && (subsetIndex < (int)num_subsets()) && "invalid subset index!");

	return GeometricObjectCollection(&m_subsets[subsetIndex]->m_elements[VERTEX],
									 &m_subsets[subsetIndex]->m_elements[EDGE],
									 &m_subsets[subsetIndex]->m_elements[FACE],
									 &m_subsets[subsetIndex]->m_elements[VOLUME]);
}

//	grid callbacks
void SubsetHandler::registered_at_grid(Grid* grid)
{
	if(m_pGrid != NULL)
		m_pGrid->unregister_observer(this);

	m_pGrid = grid;
	if(m_pGrid != NULL)
	{
	//	attach data
		grid->attach_to_vertices(m_aSubsetIndex);
		grid->attach_to_edges(m_aSubsetIndex);
		grid->attach_to_faces(m_aSubsetIndex);
		grid->attach_to_volumes(m_aSubsetIndex);

		grid->attach_to_vertices(m_aIterator);
		grid->attach_to_edges(m_aIterator);
		grid->attach_to_faces(m_aIterator);
		grid->attach_to_volumes(m_aIterator);

		m_aaSubsetIndexVRT.access(*m_pGrid, m_aSubsetIndex);
		m_aaSubsetIndexEDGE.access(*m_pGrid, m_aSubsetIndex);
		m_aaSubsetIndexFACE.access(*m_pGrid, m_aSubsetIndex);
		m_aaSubsetIndexVOL.access(*m_pGrid, m_aSubsetIndex);

		m_aaIteratorVRT.access(*m_pGrid, m_aIterator);
		m_aaIteratorEDGE.access(*m_pGrid, m_aIterator);
		m_aaIteratorFACE.access(*m_pGrid, m_aIterator);
		m_aaIteratorVOL.access(*m_pGrid, m_aIterator);

		reset_subset_indices<VertexBase>(m_pGrid->vertices_begin(), m_pGrid->vertices_end());
		reset_subset_indices<EdgeBase>(m_pGrid->edges_begin(), m_pGrid->edges_end());
		reset_subset_indices<Face>(m_pGrid->faces_begin(), m_pGrid->faces_end());
		reset_subset_indices<Volume>(m_pGrid->volumes_begin(), m_pGrid->volumes_end());
	}
}

void SubsetHandler::unregistered_from_grid(Grid* grid)
{
	assert((m_pGrid != NULL) && "ERROR in SubsetHandler::unregistered_from_grid(...): No grid assigned to SubsetHandler.");
	assert((m_pGrid == grid) && "ERROR in SubsetHandler::unregistered_from_grid(...): Grids do not match.");

//	clear the subsets
	for(uint i = 0; i < num_subsets(); ++i)
		delete m_subsets[i];

	m_subsets.clear();

	grid->detach_from_vertices(m_aSubsetIndex);
	grid->detach_from_edges(m_aSubsetIndex);
	grid->detach_from_faces(m_aSubsetIndex);
	grid->detach_from_volumes(m_aSubsetIndex);

	grid->detach_from_vertices(m_aIterator);
	grid->detach_from_edges(m_aIterator);
	grid->detach_from_faces(m_aIterator);
	grid->detach_from_volumes(m_aIterator);

	m_pGrid = NULL;
}

void SubsetHandler::elements_to_be_cleared(Grid* grid)
{
	for(uint si = 0; si < num_subsets(); si++)
	{
		for(int i = 0; i < NUM_GEOMETRIC_BASE_OBJECTS; ++i)
			m_subsets[si]->m_elements[i].clear();
	}
}

void SubsetHandler::insert_subset(int subsetIndex)
{
//	append required subsets
	if(subsetIndex >= 0)
	{
		if(subsetIndex >= (int)num_subsets())
		{
			resize_subset_vec(subsetIndex + 1);
		//	nothing else to do...
		}
		else
		{
			resize_subset_vec(num_subsets() + 1);
			move_subset(num_subsets() - 1, subsetIndex);
		}
	}
}

void SubsetHandler::erase_subset(int subsetIndex)
{
//	assign all elements of this subset to -1
//	delete the subset
//	move all subsets with higher index one entry up
//	and correct indices of assigned elements.
	if((subsetIndex >= 0) && (subsetIndex < (int)num_subsets()))
	{
		reset_subset_indices<VertexBase>(begin<VertexBase>(subsetIndex), end<VertexBase>(subsetIndex));
		reset_subset_indices<EdgeBase>(begin<EdgeBase>(subsetIndex), end<EdgeBase>(subsetIndex));
		reset_subset_indices<Face>(begin<Face>(subsetIndex), end<Face>(subsetIndex));
		reset_subset_indices<Volume>(begin<Volume>(subsetIndex), end<Volume>(subsetIndex));

		delete m_subsets[subsetIndex];

		for(uint i = subsetIndex + 1; i < num_subsets(); ++i)
		{
			int iPrev = i-1;
		//	alter indices
			set_subset_indices<VertexBase>(begin<VertexBase>(i), end<VertexBase>(i), iPrev);
			set_subset_indices<EdgeBase>(begin<EdgeBase>(i), end<EdgeBase>(i), iPrev);
			set_subset_indices<Face>(begin<Face>(i), end<Face>(i), iPrev);
			set_subset_indices<Volume>(begin<Volume>(i), end<Volume>(i), iPrev);
		//	move the subset
			m_subsets[i-1] = m_subsets[i];
			m_subsetInfos[i-1] = m_subsetInfos[i];
		}

	//	resize the subset vector
		uint numNewSubsets = num_subsets() - 1;
		m_subsets.resize(numNewSubsets);
		m_subsetInfos.resize(numNewSubsets);
	}
	else
		LOG("WARNING in SubsetHandler::erase_subset(...): bad subset index: " << subsetIndex << endl);
}

void SubsetHandler::swap_subsets(int subsetIndex1, int subsetIndex2)
{
	if((subsetIndex1 >= 0) && (subsetIndex1 < (int)num_subsets())
		&& (subsetIndex2 >= 0) && (subsetIndex2 < (int)num_subsets()))
	{
	//	set all indices of subset1 to subsetIndex2 and vice versa.
	//	swap pointers of subsets afterwards
		set_subset_indices<VertexBase>(begin<VertexBase>(subsetIndex1), end<VertexBase>(subsetIndex1), subsetIndex2);
		set_subset_indices<EdgeBase>(begin<EdgeBase>(subsetIndex1), end<EdgeBase>(subsetIndex1), subsetIndex2);
		set_subset_indices<Face>(begin<Face>(subsetIndex1), end<Face>(subsetIndex1), subsetIndex2);
		set_subset_indices<Volume>(begin<Volume>(subsetIndex1), end<Volume>(subsetIndex1), subsetIndex2);

		set_subset_indices<VertexBase>(begin<VertexBase>(subsetIndex2), end<VertexBase>(subsetIndex2), subsetIndex1);
		set_subset_indices<EdgeBase>(begin<EdgeBase>(subsetIndex2), end<EdgeBase>(subsetIndex2), subsetIndex1);
		set_subset_indices<Face>(begin<Face>(subsetIndex2), end<Face>(subsetIndex2), subsetIndex1);
		set_subset_indices<Volume>(begin<Volume>(subsetIndex2), end<Volume>(subsetIndex2), subsetIndex1);

	//	swap pointers
		Subset* pTmp = m_subsets[subsetIndex1];
		m_subsets[subsetIndex1] = m_subsets[subsetIndex2];
		m_subsets[subsetIndex2] = pTmp;

		SubsetInfo tmpSI = m_subsetInfos[subsetIndex1];
		m_subsetInfos[subsetIndex1] = m_subsetInfos[subsetIndex2];
		m_subsetInfos[subsetIndex2] = tmpSI;
	}
	else
		LOG("WARNING in SubsetHandler::swap_subsets(...): bad indices: " << subsetIndex1 << ", " << subsetIndex2 << endl);
}

void SubsetHandler::move_subset(int indexFrom, int indexTo)
{
	if((indexFrom >= 0) && (indexFrom < (int)num_subsets()) && (indexTo >= 0))
	{
		int moveDir = 0;
	//	we have to distinguish two cases
		if(indexFrom < indexTo)
			moveDir = 1;
		else if(indexTo < indexFrom)
			moveDir = -1;

		if(moveDir != 0)
		{
		//	check if we have to append subsets
			if(indexTo > (int)num_subsets())
				resize_subset_vec(indexTo + 1);

		//	store pointer to the from-subset
			Subset* pFrom = m_subsets[indexFrom];

		//	store the from-subset-info
			SubsetInfo siFrom = m_subsetInfos[indexFrom];

		//	assign new indices to elements in from subset (no iterators are changed)
			set_subset_indices<VertexBase>(begin<VertexBase>(indexFrom), end<VertexBase>(indexFrom), indexTo);
			set_subset_indices<EdgeBase>(begin<EdgeBase>(indexFrom), end<EdgeBase>(indexFrom), indexTo);
			set_subset_indices<Face>(begin<Face>(indexFrom), end<Face>(indexFrom), indexTo);
			set_subset_indices<Volume>(begin<Volume>(indexFrom), end<Volume>(indexFrom), indexTo);

		//	move the subsets between indexFrom and indexTo
			for(int i = indexFrom; i != indexTo; i+= moveDir)
			{
				int iNext = i+moveDir;
			//	assign new indices
				set_subset_indices<VertexBase>(begin<VertexBase>(iNext), end<VertexBase>(iNext), i);
				set_subset_indices<EdgeBase>(begin<EdgeBase>(iNext), end<EdgeBase>(iNext), i);
				set_subset_indices<Face>(begin<Face>(iNext), end<Face>(iNext), i);
				set_subset_indices<Volume>(begin<Volume>(iNext), end<Volume>(iNext), i);

			//	move pointer
				m_subsets[i] = m_subsets[iNext];
				m_subsetInfos[i] = m_subsetInfos[iNext];
			}

		//	assign stored pointer
			m_subsets[indexTo] = pFrom;
			m_subsetInfos[indexTo] = siFrom;
		}
	}
	else
		LOG("WARNING in SubsetHandler::move_subset(...): bad indices: " << indexFrom << ", " << indexTo << endl);
}


//	vertex callbacks
void SubsetHandler::vertex_created(Grid* grid, VertexBase* vrt, GeometricObject* pParent)
{
	assert((m_pGrid == grid) && "ERROR in SubsetHandler::vertex_created(...): Grids do not match.");
	m_aaSubsetIndexVRT[vrt] = -1;

	if((pParent != NULL) && m_bSubsetInheritanceEnabled)
		assign_subset(vrt, get_subset_index(pParent));
	else if(m_defaultSubsetIndex != -1)
		assign_subset(vrt, m_defaultSubsetIndex);
}

void SubsetHandler::vertex_to_be_erased(Grid* grid, VertexBase* vrt)
{
	assert((m_pGrid == grid) && "ERROR in SubsetHandler::vertex_to_be_erased(...): Grids do not match.");
	if(m_aaSubsetIndexVRT[vrt] != -1)
		assign_subset(vrt, -1);

}

//	edge callbacks
void SubsetHandler::edge_created(Grid* grid, EdgeBase* edge, GeometricObject* pParent)
{
	assert((m_pGrid == grid) && "ERROR in SubsetHandler::edge_created(...): Grids do not match.");
	m_aaSubsetIndexEDGE[edge] = -1;

	if((pParent != NULL) && m_bSubsetInheritanceEnabled)
		assign_subset(edge, get_subset_index(pParent));
	else if(m_defaultSubsetIndex != -1)
		assign_subset(edge, m_defaultSubsetIndex);
}

void SubsetHandler::edge_to_be_erased(Grid* grid, EdgeBase* edge)
{
	assert((m_pGrid == grid) && "ERROR in SubsetHandler::edge_to_be_erased(...): Grids do not match.");
	if(m_aaSubsetIndexEDGE[edge] != -1)
		assign_subset(edge, -1);

}

//	face callbacks
void SubsetHandler::face_created(Grid* grid, Face* face, GeometricObject* pParent)
{
	assert((m_pGrid == grid) && "ERROR in SubsetHandler::face_created(...): Grids do not match.");
	m_aaSubsetIndexFACE[face] = -1;

	if((pParent != NULL) && m_bSubsetInheritanceEnabled)
		assign_subset(face, get_subset_index(pParent));
	else if(m_defaultSubsetIndex != -1)
		assign_subset(face, m_defaultSubsetIndex);
}

void SubsetHandler::face_to_be_erased(Grid* grid, Face* face)
{
	assert((m_pGrid == grid) && "ERROR in SubsetHandler::face_to_be_erased(...): Grids do not match.");
	if(m_aaSubsetIndexFACE[face] != -1)
		assign_subset(face, -1);
}

//	volume callbacks
void SubsetHandler::volume_created(Grid* grid, Volume* vol, GeometricObject* pParent)
{
	assert((m_pGrid == grid) && "ERROR in SubsetHandler::volume_created(...): Grids do not match.");
	m_aaSubsetIndexVOL[vol] = -1;

	if((pParent != NULL) && m_bSubsetInheritanceEnabled)
		assign_subset(vol, get_subset_index(pParent));
	else if(m_defaultSubsetIndex != -1)
		assign_subset(vol, m_defaultSubsetIndex);
}

void SubsetHandler::volume_to_be_erased(Grid* grid, Volume* vol)
{
	assert((m_pGrid == grid) && "ERROR in SubsetHandler::volume_to_be_erased(...): Grids do not match.");
	if(m_aaSubsetIndexVOL[vol] != -1)
		assign_subset(vol, -1);
}

}//	end of namespace
