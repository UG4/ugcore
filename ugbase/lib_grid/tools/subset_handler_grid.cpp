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
	m_aSharedEntryVRT("SubsetHandler_SharedListEntryVRT", false),
	m_aSharedEntryEDGE("SubsetHandler_SharedListEntryEDGE", false),
	m_aSharedEntryFACE("SubsetHandler_SharedListEntryFACE", false),
	m_aSharedEntryVOL("SubsetHandler_SharedListEntryVOL", false)
{
}

GridSubsetHandler::
GridSubsetHandler(Grid& grid, uint supportedElements) :
	ISubsetHandler(supportedElements),
	m_aSharedEntryVRT("SubsetHandler_SharedListEntryVRT", false),
	m_aSharedEntryEDGE("SubsetHandler_SharedListEntryEDGE", false),
	m_aSharedEntryFACE("SubsetHandler_SharedListEntryFACE", false),
	m_aSharedEntryVOL("SubsetHandler_SharedListEntryVOL", false)
{
	assign_grid(grid);
}

GridSubsetHandler::GridSubsetHandler(const GridSubsetHandler& sh) :
	ISubsetHandler(sh.m_supportedElements),
	m_aSharedEntryVRT("SubsetHandler_SharedListEntryVRT", false),
	m_aSharedEntryEDGE("SubsetHandler_SharedListEntryEDGE", false),
	m_aSharedEntryFACE("SubsetHandler_SharedListEntryFACE", false),
	m_aSharedEntryVOL("SubsetHandler_SharedListEntryVOL", false)
{
	Grid* pGrid = sh.grid();

	if(pGrid){
//TODO: remove virtual function calls from constructor
		assign_grid(*pGrid);		
		assign_subset_handler(sh);
	}
}

GridSubsetHandler::~GridSubsetHandler()
{
	if(m_pGrid != NULL){
		erase_subset_lists_impl();
		detach_data();
	}
}

void GridSubsetHandler::grid_to_be_destroyed(Grid* grid)
{
	assert((m_pGrid == grid) && "ERROR in GridSubsetHandler::grid_to_be_destroyed(...): Grids do not match.");
	cleanup();
	ISubsetHandler::grid_to_be_destroyed(grid);
}

void GridSubsetHandler::cleanup()
{
	erase_subset_lists_impl();
	detach_data();
	//ISubsetHandler::set_grid(NULL);
}

void GridSubsetHandler::assign_grid(Grid* grid)
{
	if(m_pGrid == grid)
		return;

	if(m_pGrid)
		cleanup();

	ISubsetHandler::set_grid(grid);

	//	attach shared entries
	if(m_pGrid){
		if(elements_are_supported(SHE_VERTEX))
			m_pGrid->attach_to_vertices(m_aSharedEntryVRT);
		if(elements_are_supported(SHE_EDGE))
			m_pGrid->attach_to_edges(m_aSharedEntryEDGE);
		if(elements_are_supported(SHE_FACE))
			m_pGrid->attach_to_faces(m_aSharedEntryFACE);
		if(elements_are_supported(SHE_VOLUME))
			m_pGrid->attach_to_volumes(m_aSharedEntryVOL);
	}
}

void GridSubsetHandler::assign_grid(Grid& grid)
{
	assign_grid(&grid);
}

void GridSubsetHandler::detach_data()
{
	if(elements_are_supported(SHE_VERTEX))
		m_pGrid->detach_from_vertices(m_aSharedEntryVRT);
	if(elements_are_supported(SHE_EDGE))
		m_pGrid->detach_from_edges(m_aSharedEntryEDGE);
	if(elements_are_supported(SHE_FACE))
		m_pGrid->detach_from_faces(m_aSharedEntryFACE);
	if(elements_are_supported(SHE_VOLUME))
		m_pGrid->detach_from_volumes(m_aSharedEntryVOL);
}

GridSubsetHandler& GridSubsetHandler::operator = (const GridSubsetHandler& sh)
{
	ISubsetHandler::operator =(sh);
	return *this;
}

GridSubsetHandler& GridSubsetHandler::operator = (const ISubsetHandler& sh)
{
	ISubsetHandler::operator =(sh);
	return *this;
}

void GridSubsetHandler::erase_subset_lists()
{
	erase_subset_lists_impl();
}

void GridSubsetHandler::erase_subset_lists_impl()
{
	for(uint i = 0; i < m_subsets.size(); ++i)
		delete m_subsets[i];

	m_subsets.clear();
}

void GridSubsetHandler::clear_subset_lists(int index)
{
	if(m_pGrid)
	{
		section_container<Vertex>(index).clear();
		section_container<EdgeBase>(index).clear();
		section_container<Face>(index).clear();
		section_container<Volume>(index).clear();
	}
}

template<class TElem>
void
GridSubsetHandler::
assign_subset_impl(TElem* elem, int subsetIndex)
{
	assert((m_pGrid != NULL) && "ERROR in SubsetHandler::assign_subset(): No grid assigned to SubsetHandler.");

//	check if we have to remove elem from a subset.
	int oldIndex = get_subset_index(elem);

	if(oldIndex != -1)
		section_container<TElem>(oldIndex).erase(get_list_iterator(elem), elem->container_section());

//	add the element to the subset.
	if(subsetIndex != -1)
	{
		subset_required(subsetIndex);
		section_container<TElem>(subsetIndex).insert(elem, elem->container_section());
		subset_assigned(elem, subsetIndex);
	}
	else{
//TODO:	iterator is useless!
		subset_assigned(elem, -1);
	}

}

void GridSubsetHandler::assign_subset(Vertex* elem, int subsetIndex)
{
	assign_subset_impl(elem, subsetIndex);
}

void GridSubsetHandler::assign_subset(EdgeBase* elem, int subsetIndex)
{
	assign_subset_impl(elem, subsetIndex);
}

void GridSubsetHandler::assign_subset(Face* elem, int subsetIndex)
{
	assign_subset_impl(elem, subsetIndex);
}

void GridSubsetHandler::assign_subset(Volume* elem, int subsetIndex)
{
	assign_subset_impl(elem, subsetIndex);
}

void GridSubsetHandler::
change_subset_indices(int indOld, int indNew)
{
	if(m_pGrid)
	{
		if(elements_are_supported(SHE_VERTEX))
			change_elem_subset_indices<Vertex>(indOld, indNew);
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
			sub->m_vertices.get_container().set_pipe(
					&m_pGrid->get_attachment_pipe<Vertex>(), m_aSharedEntryVRT);
		if(elements_are_supported(SHE_EDGE))
			sub->m_edges.get_container().set_pipe(
					&m_pGrid->get_attachment_pipe<EdgeBase>(), m_aSharedEntryEDGE);
		if(elements_are_supported(SHE_FACE))
			sub->m_faces.get_container().set_pipe(
					&m_pGrid->get_attachment_pipe<Face>(), m_aSharedEntryFACE);
		if(elements_are_supported(SHE_VOLUME))
			sub->m_volumes.get_container().set_pipe(
					&m_pGrid->get_attachment_pipe<Volume>(), m_aSharedEntryVOL);
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

///	join the subset-lists but do not touch the subset-indices.
void GridSubsetHandler::join_subset_lists(int target, int src1, int src2)
{
	Subset& t = *m_subsets[target];
	Subset& s1 = *m_subsets[src1];
	Subset& s2 = *m_subsets[src2];

	if(target != src1){
		t.m_vertices.transfer_elements(s1.m_vertices);
		t.m_edges.transfer_elements(s1.m_edges);
		t.m_faces.transfer_elements(s1.m_faces);
		t.m_volumes.transfer_elements(s1.m_volumes);
	}
	if(target != src2){
		t.m_vertices.transfer_elements(s2.m_vertices);
		t.m_edges.transfer_elements(s2.m_edges);
		t.m_faces.transfer_elements(s2.m_faces);
		t.m_volumes.transfer_elements(s2.m_volumes);
	}
}

/*
void GridSubsetHandler::
register_subset_elements_at_pipe()
{
	for(size_t i = 0; i < num_subsets_in_list(); ++i)
	{
	//	register vertices
		for(VertexIterator iter = begin<Vertex>(i);
			iter != end<Vertex>(i); ++iter)
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
}*/


GridObjectCollection
GridSubsetHandler::
get_grid_objects_in_subset(int subsetIndex) const
{
//todo: replace with throw
	assert((subsetIndex >= 0) && (subsetIndex < (int)num_subsets_in_list()) && "invalid subset index!");

	subset_required(subsetIndex);
	return GridObjectCollection(&m_subsets[subsetIndex]->m_vertices,
									 &m_subsets[subsetIndex]->m_edges,
									 &m_subsets[subsetIndex]->m_faces,
									 &m_subsets[subsetIndex]->m_volumes);
}

GridObjectCollection
GridSubsetHandler::
get_grid_objects() const
{
	uint numSubsets = num_subsets_in_list();
	GridObjectCollection goc(numSubsets);
	for(uint i = 0; i < numSubsets; ++i)
	{
		goc.add_level(	&m_subsets[i]->m_vertices,
						&m_subsets[i]->m_edges,
						&m_subsets[i]->m_faces,
						&m_subsets[i]->m_volumes);
	}

	return goc;
}
/*
size_t GridSubsetHandler::
collect_subset_elements(std::vector<Vertex*>& vrtsOut, int subsetIndex) const
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
*/
}//	end of namespace
