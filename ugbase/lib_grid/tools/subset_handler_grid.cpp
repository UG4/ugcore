//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m11 d24

#include <cassert>
#include "subset_handler_grid.h"

using namespace std;

namespace ug
{
////////////////////////////////////////////////////////////////////////
//	GridSubsetHandler implementation
GridSubsetHandler::
GridSubsetHandler(uint supportedElements) :
	ISubsetHandler(supportedElements),
	m_aSharedEntry("MGSubsetHandler_SharedListEntry", false)
{
}

GridSubsetHandler::
GridSubsetHandler(Grid& grid, uint supportedElements) :
	ISubsetHandler(supportedElements),
	m_aSharedEntry("MGSubsetHandler_SharedListEntry", false)
{
	assign_grid(grid);
}

GridSubsetHandler::GridSubsetHandler(const GridSubsetHandler& sh) :
	ISubsetHandler(sh.m_supportedElements),
	m_aSharedEntry("MGSubsetHandler_SharedListEntry", false)
{
	Grid* pGrid = sh.get_assigned_grid();

	if(pGrid){
//TODO: remove virtual function calls from constructor
		assign_grid(*pGrid);		
		assign_subset_handler(sh);
	}
}

GridSubsetHandler::~GridSubsetHandler()
{
	if(m_pGrid != NULL){
		cleanup();
		m_pGrid->unregister_observer(this);
	}
}

void GridSubsetHandler::grid_to_be_destroyed(Grid* grid)
{
	cleanup();
	ISubsetHandler::grid_to_be_destroyed(grid);
}

void GridSubsetHandler::cleanup()
{
	for(size_t i = 0; i < m_subsets.size(); ++i){
		delete m_subsets[i];
	}

	if(elements_are_supported(SHE_VERTEX))
		m_pGrid->detach_from_vertices(m_aSharedEntry);
	if(elements_are_supported(SHE_EDGE))
		m_pGrid->detach_from_edges(m_aSharedEntry);
	if(elements_are_supported(SHE_FACE))
		m_pGrid->detach_from_faces(m_aSharedEntry);
	if(elements_are_supported(SHE_VOLUME))
		m_pGrid->detach_from_volumes(m_aSharedEntry);
}

void GridSubsetHandler::assign_grid(Grid& grid)
{
	if(m_pGrid)
		cleanup();

	ISubsetHandler::set_grid(&grid);

//	attach shared entries
	if(elements_are_supported(SHE_VERTEX))
		m_pGrid->attach_to_vertices(m_aSharedEntry);
	if(elements_are_supported(SHE_EDGE))
		m_pGrid->attach_to_edges(m_aSharedEntry);
	if(elements_are_supported(SHE_FACE))
		m_pGrid->attach_to_faces(m_aSharedEntry);
	if(elements_are_supported(SHE_VOLUME))
		m_pGrid->attach_to_volumes(m_aSharedEntry);
}

GridSubsetHandler& GridSubsetHandler::operator = (const GridSubsetHandler& sh)
{
	ISubsetHandler::operator =(sh);
	return *this;
}

void GridSubsetHandler::erase_subset_lists()
{
	for(uint i = 0; i < m_subsets.size(); ++i)
		delete m_subsets[i];

	m_subsets.clear();
}


void GridSubsetHandler::clear_subset_lists(int index)
{
	if(m_pGrid)
	{
		for(int i = 0; i < NUM_GEOMETRIC_BASE_OBJECTS; ++i)
			m_subsets[index]->m_elements[i].clear();
	}
}

template<class TElemPtr>
void
GridSubsetHandler::
assign_subset(TElemPtr elem, int subsetIndex, int elemType)
{
	assert((m_pGrid != NULL) && "ERROR in SubsetHandler::assign_subset(): No grid assigned to SubsetHandler.");
	subset_required(subsetIndex);

//	check if we have to remove elem from a subset.
	int oldIndex = get_subset_index(elem);

	if(oldIndex != -1)
		m_subsets[oldIndex]->m_elements[elemType].erase(get_list_iterator(elem), elem->shared_pipe_section());

//	add the element to the subset.
	if(subsetIndex != -1)
	{
		m_subsets[subsetIndex]->m_elements[elemType].insert(elem, elem->shared_pipe_section());
		subset_assigned(elem, subsetIndex);
	}
	else{
//TODO:	iterator is useless!
		subset_assigned(elem, -1);
	}

}

void GridSubsetHandler::assign_subset(VertexBase* elem, int subsetIndex)
{
	assign_subset(elem, subsetIndex, VERTEX);
}

void GridSubsetHandler::assign_subset(EdgeBase* elem, int subsetIndex)
{
	assign_subset(elem, subsetIndex, EDGE);
}

void GridSubsetHandler::assign_subset(Face* elem, int subsetIndex)
{
	assign_subset(elem, subsetIndex, FACE);
}

void GridSubsetHandler::assign_subset(Volume* elem, int subsetIndex)
{
	assign_subset(elem, subsetIndex, VOLUME);
}

void GridSubsetHandler::
change_subset_indices(int indOld, int indNew)
{
	if(m_pGrid)
	{
		if(elements_are_supported(SHE_VERTEX))
			change_elem_subset_indices<VertexBase>(indOld, indNew);
		if(elements_are_supported(SHE_EDGE))
			change_elem_subset_indices<EdgeBase>(indOld, indNew);
		if(elements_are_supported(SHE_FACE))
			change_elem_subset_indices<Face>(indOld, indNew);
		if(elements_are_supported(SHE_VOLUME)){
			change_elem_subset_indices<Volume>(indOld, indNew);
		}
	}
}				

void GridSubsetHandler::add_required_subset_lists(int maxIndex)
{
	while((int)m_subsets.size() <= maxIndex){
	//	initialize section containers
		Subset* sub = new Subset;
		if(elements_are_supported(SHE_VERTEX))
			sub->m_elements[VERTEX].get_container().set_pipe(
					&m_pGrid->get_attachment_pipe<VertexBase>(), m_aSharedEntry);
		if(elements_are_supported(SHE_EDGE))
			sub->m_elements[EDGE].get_container().set_pipe(
					&m_pGrid->get_attachment_pipe<EdgeBase>(), m_aSharedEntry);
		if(elements_are_supported(SHE_FACE))
			sub->m_elements[FACE].get_container().set_pipe(
					&m_pGrid->get_attachment_pipe<Face>(), m_aSharedEntry);
		if(elements_are_supported(SHE_VOLUME))
			sub->m_elements[VOLUME].get_container().set_pipe(
					&m_pGrid->get_attachment_pipe<Volume>(), m_aSharedEntry);
		m_subsets.push_back(sub);
	}
}

void GridSubsetHandler::erase_subset_lists(int index)
{
	delete m_subsets[index];
	for(uint i = index + 1; i < num_subsets_in_list(); ++i)
		m_subsets[i-1] = m_subsets[i];
	m_subsets.resize(m_subsets.size() - 1);
}

void GridSubsetHandler::swap_subset_lists(int ind1, int ind2)
{
	Subset* pTmp = m_subsets[ind1];
	m_subsets[ind1] = m_subsets[ind2];
	m_subsets[ind2] = pTmp;
}

void GridSubsetHandler::move_subset_lists(int indexFrom, int indexTo)
{
	int moveDir = 0;
//	we have to distinguish two cases
	if(indexFrom < indexTo)
		moveDir = 1;
	else if(indexTo < indexFrom)
		moveDir = -1;

	if(moveDir != 0)
	{
	//	store pointer to the from-subset
		Subset* pFrom = m_subsets[indexFrom];

		for(int i = indexFrom; i != indexTo; i+= moveDir)
		{
			int iNext = i+moveDir;

		//	move pointer
			m_subsets[i] = m_subsets[iNext];
		}

	//	assign stored pointer
		m_subsets[indexTo] = pFrom;
	}
}

void GridSubsetHandler::
register_subset_elements_at_pipe()
{
	for(size_t i = 0; i < num_subsets_in_list(); ++i)
	{
	//	register vertices
		for(VertexBaseIterator iter = begin<VertexBase>(i);
			iter != end<VertexBase>(i); ++iter)
			register_at_pipe(*iter);

	//	register edges
		for(EdgeBaseIterator iter = begin<EdgeBase>(i);
			iter != end<EdgeBase>(i); ++iter)
			register_at_pipe(*iter);

	//	register faces
		for(FaceIterator iter = begin<Face>(i);
			iter != end<Face>(i); ++iter)
			register_at_pipe(*iter);

	//	register volumes
		for(VolumeIterator iter = begin<Volume>(i);
			iter != end<Volume>(i); ++iter)
			register_at_pipe(*iter);
	}
}


GeometricObjectCollection
GridSubsetHandler::
get_geometric_object_collection(int subsetIndex)
{
	assert((subsetIndex >= 0) && (subsetIndex < (int)num_subsets_in_list()) && "invalid subset index!");

	return GeometricObjectCollection(&m_subsets[subsetIndex]->m_elements[VERTEX],
									 &m_subsets[subsetIndex]->m_elements[EDGE],
									 &m_subsets[subsetIndex]->m_elements[FACE],
									 &m_subsets[subsetIndex]->m_elements[VOLUME]);
}

GeometricObjectCollection
GridSubsetHandler::
get_geometric_object_collection()
{
	uint numSubsets = num_subsets_in_list();
	GeometricObjectCollection goc(numSubsets);
	for(uint i = 0; i < numSubsets; ++i)
	{
		goc.add_level(	&m_subsets[i]->m_elements[VERTEX],
						&m_subsets[i]->m_elements[EDGE],
						&m_subsets[i]->m_elements[FACE],
						&m_subsets[i]->m_elements[VOLUME]);
	}
	
	return goc;
}

size_t GridSubsetHandler::
collect_subset_elements(std::vector<VertexBase*>& vrtsOut, int subsetIndex) const
{
	return collect_subset_elements_impl(vrtsOut, subsetIndex);
}

size_t GridSubsetHandler::
collect_subset_elements(std::vector<EdgeBase*>& edgesOut, int subsetIndex) const
{
	return collect_subset_elements_impl(edgesOut, subsetIndex);
}

size_t GridSubsetHandler::
collect_subset_elements(std::vector<Face*>& facesOut, int subsetIndex) const
{
	return collect_subset_elements_impl(facesOut, subsetIndex);
}

size_t GridSubsetHandler::
collect_subset_elements(std::vector<Volume*>& volsOut, int subsetIndex) const
{
	return collect_subset_elements_impl(volsOut, subsetIndex);
}

template <class TElem>
size_t GridSubsetHandler::
collect_subset_elements_impl(std::vector<TElem*>& elemsOut, int subsetIndex) const
{
	typedef typename geometry_traits<TElem>::iterator ElemIter;
	typedef typename geometry_traits<TElem>::const_iterator ConstElemIter;

	elemsOut.clear();
	
	if(!m_pGrid){
		return 0;
	}
		
	if(subsetIndex < 0){
	//	iterate over all elements of the underlying grid and compare indices
		Grid& grid = *m_pGrid;
		
		for(ElemIter iter = grid.begin<TElem>(); iter != grid.end<TElem>(); ++iter)
		{
			if(get_subset_index(*iter) == -1)
				elemsOut.push_back(*iter);
		}
	}
	else{
		elemsOut.reserve(num<TElem>(subsetIndex));
		for(ConstElemIter iter = begin<TElem>(subsetIndex);
			iter != end<TElem>(subsetIndex); ++iter)
		{
			elemsOut.push_back(*iter);
		}
	}
	
	return elemsOut.size();
}

}//	end of namespace
