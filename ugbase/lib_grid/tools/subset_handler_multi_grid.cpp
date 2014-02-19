//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m11 d24

#include <cassert>
#include "subset_handler_multi_grid.h"

using namespace std;

namespace ug
{
////////////////////////////////////////////////////////////////////////
//	MultiGridSubsetHandler implementation
MultiGridSubsetHandler::
MultiGridSubsetHandler(uint supportedElements) :
	ISubsetHandler(supportedElements),
	m_aSharedEntryVRT("MGSubsetHandler_SharedListEntryVRT", false),
	m_aSharedEntryEDGE("MGSubsetHandler_SharedListEntryEDGE", false),
	m_aSharedEntryFACE("MGSubsetHandler_SharedListEntryFACE", false),
	m_aSharedEntryVOL("MGSubsetHandler_SharedListEntryVOL", false)
{
	m_numSubsets = 0;
	m_pMG = NULL;
}

MultiGridSubsetHandler::
MultiGridSubsetHandler(MultiGrid& mg, uint supportedElements) :
	ISubsetHandler(supportedElements),
	m_aSharedEntryVRT("MGSubsetHandler_SharedListEntryVRT", false),
	m_aSharedEntryEDGE("MGSubsetHandler_SharedListEntryEDGE", false),
	m_aSharedEntryFACE("MGSubsetHandler_SharedListEntryFACE", false),
	m_aSharedEntryVOL("MGSubsetHandler_SharedListEntryVOL", false)
{
	m_numSubsets = 0;
	m_pMG = NULL;
	assign_grid(mg);
}

MultiGridSubsetHandler::
MultiGridSubsetHandler(const MultiGridSubsetHandler& sh) :
	ISubsetHandler(sh.m_supportedElements),
	m_aSharedEntryVRT("MGSubsetHandler_SharedListEntryVRT", false),
	m_aSharedEntryEDGE("MGSubsetHandler_SharedListEntryEDGE", false),
	m_aSharedEntryFACE("MGSubsetHandler_SharedListEntryFACE", false),
	m_aSharedEntryVOL("MGSubsetHandler_SharedListEntryVOL", false)
{
	MultiGrid* pGrid = sh.m_pMG;

	if(pGrid){
//TODO: remove virtual function calls from constructor
		assign_grid(*pGrid);		
		assign_subset_handler(sh);
	}
}

MultiGridSubsetHandler::~MultiGridSubsetHandler()
{
	if(m_pGrid != NULL){
		erase_subset_lists_impl();
		detach_data();
	}
}

void MultiGridSubsetHandler::grid_to_be_destroyed(Grid* grid)
{
	assert((m_pGrid == grid) && "ERROR in MultiGridSubsetHandler::grid_to_be_destroyed(...): Grids do not match.");
	cleanup();
	ISubsetHandler::grid_to_be_destroyed(grid);
}

void MultiGridSubsetHandler::cleanup()
{
	erase_subset_lists_impl();
	m_levels.clear();

	if(m_pGrid){
		detach_data();

	//	unregister the previously registered callback
		m_callbackId = MessageHub::SPCallbackId(NULL);

		m_pMG = NULL;
	}

	//ISubsetHandler::set_grid(NULL);
}

void MultiGridSubsetHandler::assign_grid(MultiGrid& mg)
{
	if(m_pMG == &mg)
		return;

	if(m_pMG)
		cleanup();

	m_pMG = &mg;

	ISubsetHandler::set_grid(&mg);

//	attach shared entries
	if(elements_are_supported(SHE_VERTEX))
		m_pGrid->attach_to_vertices(m_aSharedEntryVRT);
	if(elements_are_supported(SHE_EDGE))
		m_pGrid->attach_to_edges(m_aSharedEntryEDGE);
	if(elements_are_supported(SHE_FACE))
		m_pGrid->attach_to_faces(m_aSharedEntryFACE);
	if(elements_are_supported(SHE_VOLUME))
		m_pGrid->attach_to_volumes(m_aSharedEntryVOL);

//	register the callback
	m_callbackId = m_pMG->message_hub()->register_class_callback(
						this, &ug::MultiGridSubsetHandler::multigrid_changed);

	level_required(m_pMG->num_levels());
}

void MultiGridSubsetHandler::detach_data()
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

void MultiGridSubsetHandler::erase_subset_lists()
{
	erase_subset_lists_impl();
}

void MultiGridSubsetHandler::erase_subset_lists_impl()
{
	for(size_t level = 0; level < m_levels.size(); ++level)
	{
		for(size_t i = 0; i < m_levels[level].size(); ++i)
			delete m_levels[level][i];

		m_levels[level].clear();
	}

	m_numSubsets = 0;
}

void MultiGridSubsetHandler::clear_subset_lists(int index)
{
	if(m_pGrid)
	{
		for(size_t level = 0; level < m_levels.size(); ++level){
			section_container<Vertex>(index, level).clear();
			section_container<EdgeBase>(index, level).clear();
			section_container<Face>(index, level).clear();
			section_container<Volume>(index, level).clear();
		}
	}
}

template<class TElem>
void
MultiGridSubsetHandler::
assign_subset_impl(TElem* elem, int subsetIndex)
{
	assert((m_pGrid != NULL) && "ERROR in SubsetHandler::assign_subset(): No grid assigned to SubsetHandler.");

	int level = m_pMG->get_level(elem);

//	check if we have to remove elem from a subset.
	int oldIndex = get_subset_index(elem);
	
	if(oldIndex != -1)
		section_container<TElem>(oldIndex, level).erase(get_list_iterator(elem), elem->container_section());

//	add the element to the subset.
	if(subsetIndex != -1)
	{
		subset_required(subsetIndex);
		level_required(level);
		section_container<TElem>(subsetIndex, level).insert(elem, elem->container_section());
		subset_assigned(elem, subsetIndex);
	}
	else {
//TODO:	iterator is useless!
		subset_assigned(elem, -1);
	}
}

void MultiGridSubsetHandler::assign_subset(Vertex* elem, int subsetIndex)
{
	assign_subset_impl(elem, subsetIndex);
}

void MultiGridSubsetHandler::assign_subset(EdgeBase* elem, int subsetIndex)
{
	assign_subset_impl(elem, subsetIndex);
}

void MultiGridSubsetHandler::assign_subset(Face* elem, int subsetIndex)
{
	assign_subset_impl(elem, subsetIndex);
}

void MultiGridSubsetHandler::assign_subset(Volume* elem, int subsetIndex)
{
	assign_subset_impl(elem, subsetIndex);
}

void MultiGridSubsetHandler::
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
		if(elements_are_supported(SHE_VOLUME))
			change_elem_subset_indices<Volume>(indOld, indNew);
	}
}

void MultiGridSubsetHandler::add_required_subset_lists(int maxIndex)
{
	while((int)num_subsets_in_list() <= maxIndex)
		add_subset_to_all_levels();
}

void MultiGridSubsetHandler::erase_subset_lists(int index)
{
	for(size_t level = 0; level < m_levels.size(); ++level)
	{
		SubsetVec& subsets = m_levels[level];
		delete subsets[index];
		for(uint i = index + 1; i < num_subsets_in_list(); ++i)
			subsets[i-1] = subsets[i];
		subsets.resize(subsets.size() - 1);
	}
	m_numSubsets--;
}

void MultiGridSubsetHandler::swap_subset_lists(int ind1, int ind2)
{
	for(size_t level = 0; level < m_levels.size(); ++level)
	{
		SubsetVec& subsets = m_levels[level];
		Subset* pTmp = subsets[ind1];
		subsets[ind1] = subsets[ind2];
		subsets[ind2] = pTmp;
	}
}

void MultiGridSubsetHandler::move_subset_lists(int indexFrom, int indexTo)
{
	int moveDir = 0;
//	we have to distinguish two cases
	if(indexFrom < indexTo)
		moveDir = 1;
	else if(indexTo < indexFrom)
		moveDir = -1;

	if(moveDir != 0)
	{
		for(size_t level = 0; level < m_levels.size(); ++level)
		{
			SubsetVec& subsets = m_levels[level];

		//	store pointer to the from-subset
			Subset* pFrom = subsets[indexFrom];

			for(int i = indexFrom; i != indexTo; i+= moveDir)
			{
				int iNext = i+moveDir;

			//	move pointer
				subsets[i] = subsets[iNext];
			}

		//	assign stored pointer
			subsets[indexTo] = pFrom;
		}
	}
}

void MultiGridSubsetHandler::join_subset_lists(int target, int src1, int src2)
{
	for(size_t level = 0; level < m_levels.size(); ++level){
		SubsetVec& subsets = m_levels[level];

	//	store pointer to the from-subset
		Subset& t = *subsets[target];
		Subset& s1 = *subsets[src1];
		Subset& s2 = *subsets[src2];

		if(target != src1){
			t.m_vertices.append(s1.m_vertices);
			t.m_edges.append(s1.m_edges);
			t.m_faces.append(s1.m_faces);
			t.m_volumes.append(s1.m_volumes);
		}
		if(target != src2){
			t.m_vertices.append(s2.m_vertices);
			t.m_edges.append(s2.m_edges);
			t.m_faces.append(s2.m_faces);
			t.m_volumes.append(s2.m_volumes);
		}
	}

	if(target != src1)
		clear_subset_lists(src1);
	if(target != src2)
		clear_subset_lists(src2);
}

/*
void MultiGridSubsetHandler::
register_subset_elements_at_pipe()
{
	for(size_t l = 0; l < num_levels(); ++l)
	{
		for(size_t i = 0; i < num_subsets_in_list(); ++i)
		{
		//	register vertices
			for(VertexIterator iter = begin<Vertex>(i, l);
				iter != end<Vertex>(i, l); ++iter)
				register_at_pipe(*iter);

		//	register edges
			for(EdgeBaseIterator iter = begin<EdgeBase>(i, l);
				iter != end<EdgeBase>(i, l); ++iter)
				register_at_pipe(*iter);

		//	register faces
			for(FaceIterator iter = begin<Face>(i, l);
				iter != end<Face>(i, l); ++iter)
				register_at_pipe(*iter);

		//	register volumes
			for(VolumeIterator iter = begin<Volume>(i, l);
				iter != end<Volume>(i, l); ++iter)
				register_at_pipe(*iter);
		}
	}
}
*/

GridObjectCollection
MultiGridSubsetHandler::
get_grid_objects(int subsetIndex, int level) const
{
	subset_required(subsetIndex);
	level_required(level);

	return GridObjectCollection(&m_levels[level][subsetIndex]->m_vertices,
									 &m_levels[level][subsetIndex]->m_edges,
									 &m_levels[level][subsetIndex]->m_faces,
									 &m_levels[level][subsetIndex]->m_volumes);
}

GridObjectCollection
MultiGridSubsetHandler::
get_grid_objects_in_subset(int subsetIndex) const
{
	subset_required(subsetIndex);
	GridObjectCollection goc(m_levels.size());
	for(size_t i = 0; i < m_levels.size(); ++i)
	{
		goc.add_level(	&m_levels[i][subsetIndex]->m_vertices,
						&m_levels[i][subsetIndex]->m_edges,
						&m_levels[i][subsetIndex]->m_faces,
						&m_levels[i][subsetIndex]->m_volumes);
	}
	
	return goc;
}

GridObjectCollection
MultiGridSubsetHandler::
get_grid_objects_in_level(int level) const
{
	level_required(level);
	uint numSubsets = num_subsets_in_list();
	GridObjectCollection goc(numSubsets);
	for(uint i = 0; i < numSubsets; ++i)
	{
		goc.add_level(	&m_levels[level][i]->m_vertices,
						&m_levels[level][i]->m_edges,
						&m_levels[level][i]->m_faces,
						&m_levels[level][i]->m_volumes);
	}
	
	return goc;
}

MultiGridSubsetHandler::Subset* MultiGridSubsetHandler::new_subset()
{

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
	return sub;
}

void MultiGridSubsetHandler::add_level()
{
	m_levels.push_back(SubsetVec());
	int topLevel = m_levels.size() - 1;
	for(uint i = 0; i < num_subsets_in_list(); ++i)
		m_levels[topLevel].push_back(new_subset());
}

void MultiGridSubsetHandler::add_subset_to_all_levels()
{
	for(size_t i = 0; i < m_levels.size(); ++i)
		m_levels[i].push_back(new_subset());

	m_numSubsets++;
}

void MultiGridSubsetHandler::
multigrid_changed(const GridMessage_MultiGridChanged& gm)
{
	if(gm.message_type() == GMMGCT_LEVEL_ADDED)
		level_required(gm.num_levels_in_grid() - 1);
}

/*
size_t MultiGridSubsetHandler::
collect_subset_elements(std::vector<Vertex*>& vrtsOut, int subsetIndex) const
{
	return collect_subset_elements_impl(vrtsOut, subsetIndex);
}

size_t MultiGridSubsetHandler::
collect_subset_elements(std::vector<EdgeBase*>& edgesOut, int subsetIndex) const
{
	return collect_subset_elements_impl(edgesOut, subsetIndex);
}

size_t MultiGridSubsetHandler::
collect_subset_elements(std::vector<Face*>& facesOut, int subsetIndex) const
{
	return collect_subset_elements_impl(facesOut, subsetIndex);
}

size_t MultiGridSubsetHandler::
collect_subset_elements(std::vector<Volume*>& volsOut, int subsetIndex) const
{
	return collect_subset_elements_impl(volsOut, subsetIndex);
}

template <class TElem>
size_t MultiGridSubsetHandler::
collect_subset_elements_impl(std::vector<TElem*>& elemsOut, int subsetIndex) const
{
	typedef typename geometry_traits<TElem>::iterator ElemIter;
	typedef typename geometry_traits<TElem>::const_iterator ConstElemIter;
	
	elemsOut.clear();
	if(!m_pGrid)
		return 0;
		
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
		for(size_t i = 0; i < num_levels(); ++i){
			for(ConstElemIter iter = begin<TElem>(subsetIndex, i);
				iter != end<TElem>(subsetIndex, i); ++iter)
			{
				elemsOut.push_back(*iter);
			}
		}
	}
	
	return elemsOut.size();
}
*/
}//	end of namespace
